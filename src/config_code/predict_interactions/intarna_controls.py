"""Run and resume sequential IntaRNA controls from a prepared control plan."""

import argparse
import contextlib
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import signal
import statistics
import subprocess
import sys
import time
from datetime import datetime, timezone
import uuid


CSV_FIELDS = "id1,start1,end1,id2,start2,end2,E,ED1,ED2".split(",")
SUCCESS_STATES = {"completed", "completed_no_interaction"}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path, value):
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def unique_directory(parent, prefix):
    parent.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    path = parent / f"{prefix}_{stamp}_{uuid.uuid4().hex[:8]}"
    path.mkdir()
    return path


def read_fasta(path):
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    if not lines or not lines[0].startswith(">") or sum(line.startswith(">") for line in lines) != 1:
        raise ValueError(f"Expected one FASTA record in {path}")
    identifier = lines[0][1:].split()[0]
    sequence = "".join(line.strip() for line in lines[1:])
    if not sequence or set(sequence) - set("ACGU"):
        raise ValueError(f"Expected unambiguous uppercase RNA in {path}")
    return identifier, sequence


def validate_result(path, human_id, human_length, viral_id, viral_length):
    with Path(path).open(newline="", encoding="utf-8-sig") as source:
        reader = csv.DictReader(source, delimiter=";")
        if reader.fieldnames != CSV_FIELDS:
            raise ValueError("Prediction CSV header differs from the requested nine columns")
        rows = list(reader)
    if len(rows) > 1:
        raise ValueError("Expected at most one interaction")
    if not rows:
        return {"status": "completed_no_interaction", "prediction_rows": 0, "energy": None, "fields": None}
    row = rows[0]
    if set(row) != set(CSV_FIELDS) or any(value is None for value in row.values()):
        raise ValueError("Prediction CSV contains missing or extra fields")
    if row["id1"] != human_id or row["id2"] != viral_id:
        raise ValueError("Prediction CSV IDs differ from the requested input sequences")
    start1, end1, start2, end2 = [int(row[key]) for key in ("start1", "end1", "start2", "end2")]
    if not (1 <= start1 <= end1 <= human_length and 1 <= start2 <= end2 <= viral_length):
        raise ValueError("Prediction coordinates lie outside the full input sequences")
    energies = [float(row[key]) for key in ("E", "ED1", "ED2")]
    if not all(math.isfinite(value) for value in energies):
        raise ValueError("Prediction energies must be finite")
    return {"status": "completed", "prediction_rows": 1, "energy": energies[0], "fields": row}


def parse_versions(text):
    intarna = re.search(r"IntaRNA\b[^\r\n0-9]*([0-9]+\.[0-9]+\.[0-9]+)", text, re.IGNORECASE)
    vienna = re.search(r"Vienna\s*RNA\s+package\s+([0-9]+\.[0-9]+\.[0-9]+)", text, re.IGNORECASE)
    if intarna is None or vienna is None:
        raise ValueError("Version output must identify IntaRNA and its ViennaRNA package")
    return {"intarna": intarna.group(1), "viennarna": vienna.group(1)}


def bundle_path(bundle, relative):
    path = Path(relative)
    if path.is_absolute():
        raise ValueError("Control plan paths must be relative to the bundle")
    resolved = (bundle / path).resolve()
    try:
        resolved.relative_to(bundle)
    except ValueError:
        raise ValueError("Control plan path leaves the bundle") from None
    return resolved


def verify_sequence_entry(bundle, entry, identifier):
    path = bundle_path(bundle, entry["fasta"])
    if sha256(path) != entry["file_sha256"]:
        raise ValueError(f"Input FASTA checksum differs for {identifier}")
    observed_id, sequence = read_fasta(path)
    if observed_id != identifier or len(sequence) != entry["length_nt"]:
        raise ValueError(f"Input FASTA ID or length differs for {identifier}")
    if hashlib.sha256(sequence.encode("ascii")).hexdigest() != entry["sequence_sha256"]:
        raise ValueError(f"Input RNA checksum differs for {identifier}")


def load_plan(bundle):
    bundle = Path(bundle).resolve()
    path = bundle / "control_plan.json"
    plan = json.loads(path.read_text(encoding="utf-8"))
    if plan["schema_version"] != 1 or plan["num_shuffles"] != 100:
        raise ValueError("Expected a version-one plan with exactly 100 shuffled controls")
    if plan["expected_versions"] != {"intarna": "3.4.1", "viennarna": "2.7.0"}:
        raise ValueError("The plan must match IntaRNA 3.4.1 and ViennaRNA 2.7.0")
    required = ["--model=X", "--mode=H", "--acc=C", "--tAccW=150", "--tAccL=100",
                "--qAccW=150", "--qAccL=100", "--seedBP=7", "--outNumber=1", "--outMode=C",
                "--outCsvCols=" + ",".join(CSV_FIELDS)]
    if plan["settings_args"] != required:
        raise ValueError("Control settings differ from the original local prediction settings")
    humans = plan["humans"]
    human_ids = [human["transcript_id"] for human in humans]
    if not human_ids or len(set(human_ids)) != len(human_ids):
        raise ValueError("Human transcript IDs must be unique and nonempty")
    if any(re.fullmatch(r"[A-Za-z0-9_.-]+", identifier) is None for identifier in human_ids):
        raise ValueError("Human transcript IDs must be safe as filenames")
    shuffles = plan["shuffles"]
    if [query["shuffle_id"] for query in shuffles] != list(range(1, 101)):
        raise ValueError("The shuffled inputs must contain IDs 1 through 100 exactly once")
    if len({query["id"] for query in shuffles}) != 100:
        raise ValueError("Native shuffled FASTA IDs must be unique")
    original_viral = plan["original_viral"]
    verify_sequence_entry(bundle, original_viral, original_viral["id"])
    original_log = plan["original_version_log"]
    log_path = bundle_path(bundle, original_log["path"])
    if sha256(log_path) != original_log["file_sha256"]:
        raise ValueError("Original version log checksum differs")
    if parse_versions(log_path.read_text(encoding="utf-8")) != plan["expected_versions"]:
        raise ValueError("Original prediction versions differ from the plan")
    for human in humans:
        verify_sequence_entry(bundle, human, human["transcript_id"])
        original_csv = bundle_path(bundle, human["original_csv"])
        original_command = bundle_path(bundle, human["original_command"])
        if sha256(original_csv) != human["original_csv_sha256"] or sha256(original_command) != human["original_command_sha256"]:
            raise ValueError("Original prediction provenance checksum differs")
        result = validate_result(original_csv, human["transcript_id"], human["length_nt"],
                                 original_viral["id"], original_viral["length_nt"])
        if result["status"] != "completed" or result["energy"] != float(human["original_energy"]):
            raise ValueError("Original energy differs from its verified prediction CSV")
        arguments = shlex.split(original_command.read_text(encoding="utf-8"))
        settings = [argument for argument in arguments[1:]
                    if not argument.startswith(("--target=", "--query=", "--out="))]
        if settings != plan["settings_args"]:
            raise ValueError("Original command settings differ from the control plan")
    for query in shuffles:
        verify_sequence_entry(bundle, query, query["id"])
    plan["_bundle"] = str(bundle)
    plan["_plan_sha256"] = sha256(path)
    return plan


def inspect_versions(executable, invocation_dir):
    result = subprocess.run([str(executable), "--version"], capture_output=True, text=True, check=False)
    (invocation_dir / "version.stdout.log").write_text(result.stdout, encoding="utf-8")
    (invocation_dir / "version.stderr.log").write_text(result.stderr, encoding="utf-8")
    version_log = invocation_dir / "version.log"
    version_log.write_text(result.stdout + result.stderr, encoding="utf-8")
    if result.returncode != 0:
        raise ValueError("IntaRNA version probe failed")
    versions = parse_versions(result.stdout + result.stderr)
    versions.update({"version_log_sha256": sha256(version_log), "executable_sha256": sha256(executable)})
    return versions


def job_fingerprint(plan, human, query, versions):
    return {
        "plan_sha256": plan["_plan_sha256"], "settings_args": plan["settings_args"],
        "versions": {key: versions[key] for key in ("intarna", "viennarna")},
        "human": {key: human[key] for key in ("transcript_id", "fasta", "file_sha256", "sequence_sha256", "length_nt")},
        "shuffle": {key: query[key] for key in ("shuffle_id", "id", "fasta", "file_sha256", "sequence_sha256", "length_nt")},
    }


def verified_attempt(attempt_dir, expected):
    marker = attempt_dir / "completed.json"
    record = json.loads(marker.read_text(encoding="utf-8"))
    if record["status"] not in SUCCESS_STATES or record["return_code"] != 0 or record["fingerprint"] != expected:
        raise ValueError("Attempt metadata differs from the requested inputs, settings, or versions")
    required = {"result.csv", "stdout.log", "stderr.log", "command.json", "command.sh", "request.json", "version.log"}
    if set(record["artifact_sha256"]) != required:
        raise ValueError("Attempt provenance is incomplete")
    for name, digest in record["artifact_sha256"].items():
        if sha256(attempt_dir / name) != digest:
            raise ValueError(f"Attempt artifact checksum differs: {name}")
    if json.loads((attempt_dir / "request.json").read_text()) != expected:
        raise ValueError("Attempt input request differs")
    if parse_versions((attempt_dir / "version.log").read_text()) != expected["versions"]:
        raise ValueError("Attempt version log differs")
    if record["versions"]["version_log_sha256"] != sha256(attempt_dir / "version.log"):
        raise ValueError("Attempt version metadata checksum differs")
    if any(record["versions"][key] != value for key, value in expected["versions"].items()):
        raise ValueError("Attempt recorded versions differ")
    if re.fullmatch(r"[0-9a-f]{64}", record["versions"]["executable_sha256"]) is None:
        raise ValueError("Attempt executable checksum is missing or invalid")
    command = json.loads((attempt_dir / "command.json").read_text())
    settings = [argument for argument in command[1:]
                if not argument.startswith(("--target=", "--query=", "--out="))]
    if settings != expected["settings_args"] or record["command"] != command:
        raise ValueError("Attempt command settings differ")
    human, query = expected["human"], expected["shuffle"]
    for prefix, suffix in (("--target=", human["fasta"]), ("--query=", query["fasta"]),
                           ("--out=", attempt_dir.name + "/result.csv.partial")):
        values = [argument[len(prefix):].replace("\\", "/") for argument in command if argument.startswith(prefix)]
        if len(values) != 1 or not values[0].endswith("/" + suffix.replace("\\", "/")):
            raise ValueError("Attempt command input or output paths differ")
    result = validate_result(attempt_dir / "result.csv", human["transcript_id"], human["length_nt"], query["id"], query["length_nt"])
    if any(record[key] != result[key] for key in ("status", "prediction_rows", "energy", "fields")):
        raise ValueError("Attempt outcome differs from the saved result CSV")
    record["attempt_dir"] = str(attempt_dir)
    return record


def find_verified_attempt(job_dir, expected):
    invalid = []
    if not job_dir.exists():
        return None, invalid
    for attempt_dir in sorted(job_dir.glob("attempt_*"), reverse=True):
        try:
            return verified_attempt(attempt_dir, expected), invalid
        except (OSError, ValueError, KeyError, TypeError) as error:
            invalid.append({"attempt_dir": str(attempt_dir), "reason": f"{type(error).__name__}: {error}"})
    return None, invalid


def run_process(command, stdout_path, stderr_path):
    with stdout_path.open("xb") as stdout, stderr_path.open("xb") as stderr:
        process = subprocess.Popen(command, stdout=stdout, stderr=stderr)
        try:
            return process.wait()
        except BaseException:
            if process.poll() is None:
                process.terminate()
                try:
                    process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
            raise


def run_job(plan, human, query, versions, executable, invocation_dir):
    bundle = Path(plan["_bundle"])
    job_dir = bundle / "results/jobs" / human["transcript_id"] / f"shuffle_{query['shuffle_id']:03d}"
    attempt_dir = unique_directory(job_dir, "attempt")
    expected = job_fingerprint(plan, human, query, versions)
    partial = attempt_dir / "result.csv.partial"
    command = [str(executable), f"--target={bundle_path(bundle, human['fasta'])}",
               f"--query={bundle_path(bundle, query['fasta'])}", *plan["settings_args"], f"--out={partial}"]
    (attempt_dir / "command.json").write_text(json.dumps(command, indent=2) + "\n", encoding="utf-8")
    (attempt_dir / "command.sh").write_text(shlex.join(command) + "\n", encoding="utf-8")
    (attempt_dir / "request.json").write_text(json.dumps(expected, indent=2) + "\n", encoding="utf-8")
    shutil.copyfile(invocation_dir / "version.log", attempt_dir / "version.log")
    started = time.monotonic()
    record = {"status": "running", "started_utc": utc_now(), "fingerprint": expected,
              "command": command, "versions": versions, "attempt_dir": str(attempt_dir)}
    atomic_json(attempt_dir / "state.json", record)
    try:
        verify_sequence_entry(bundle, human, human["transcript_id"])
        verify_sequence_entry(bundle, query, query["id"])
        if sha256(executable) != versions["executable_sha256"]:
            raise ValueError("IntaRNA executable changed after its version check")
        if sha256(attempt_dir / "version.log") != versions["version_log_sha256"]:
            raise ValueError("Version log changed before this attempt")
        code = run_process(command, attempt_dir / "stdout.log", attempt_dir / "stderr.log")
        record["return_code"] = code
        if code != 0:
            raise RuntimeError(f"IntaRNA exited with code {code}")
        verify_sequence_entry(bundle, human, human["transcript_id"])
        verify_sequence_entry(bundle, query, query["id"])
        if sha256(executable) != versions["executable_sha256"]:
            raise ValueError("IntaRNA executable changed during this attempt")
        result = validate_result(partial, human["transcript_id"], human["length_nt"], query["id"], query["length_nt"])
        partial.rename(attempt_dir / "result.csv")
        record.update(result)
        artifacts = ("result.csv", "stdout.log", "stderr.log", "command.json", "command.sh", "request.json", "version.log")
        record["artifact_sha256"] = {name: sha256(attempt_dir / name) for name in artifacts}
        record["finished_utc"] = utc_now()
        record["elapsed_seconds"] = round(time.monotonic() - started, 3)
        atomic_json(attempt_dir / "state.json", record)
        atomic_json(attempt_dir / "completed.json", record)
        return record
    except BaseException as error:
        record.update({"status": "interrupted" if isinstance(error, KeyboardInterrupt) else "failed",
                       "error": f"{type(error).__name__}: {error}", "finished_utc": utc_now(),
                       "elapsed_seconds": round(time.monotonic() - started, 3)})
        atomic_json(attempt_dir / "state.json", record)
        if isinstance(error, KeyboardInterrupt):
            raise
        return record


def summarize_results(plan, outcomes):
    summaries = []
    for human in plan["humans"]:
        reported = []
        no_interaction = failed = missing = unverified = 0
        for query in plan["shuffles"]:
            record = outcomes.get((human["transcript_id"], query["shuffle_id"]), {})
            state = record.get("status", "missing")
            if state == "completed":
                reported.append(float(record["energy"]))
            elif state == "completed_no_interaction":
                no_interaction += 1
            elif state == "missing":
                missing += 1
            elif state == "unverified":
                unverified += 1
            else:
                failed += 1
        verified = len(reported) + no_interaction
        complete = verified == plan["num_shuffles"]
        extended = reported + [math.inf] * no_interaction
        mean = statistics.mean(extended) if extended else None
        median = statistics.median(extended) if extended else None
        controls_le = sum(energy <= float(human["original_energy"]) for energy in reported)
        summaries.append({
            "transcript_id": human["transcript_id"], "related_group": human.get("related_group", ""),
            "original_energy": float(human["original_energy"]), "status": "complete" if complete else "pending",
            "verified_controls": verified, "completed_with_interaction": len(reported),
            "completed_no_interaction": no_interaction, "failed": failed, "unverified": unverified, "missing": missing,
            "shuffled_mean": "inf" if mean == math.inf else mean,
            "shuffled_median": "inf" if median == math.inf else median,
            "reported_energy_mean": statistics.mean(reported) if reported else None,
            "reported_energy_median": statistics.median(reported) if reported else None,
            "controls_le_original": controls_le if complete else None,
            "empirical_tail_fraction": (controls_le + 1) / (plan["num_shuffles"] + 1) if complete else None,
            "observed_controls_le_original_so_far": controls_le,
            "evidence_label": "exploratory",
        })
    return summaries


def latest_unverified_state(job_dir, invalid):
    if not invalid:
        return {"status": "missing"}
    for attempt_dir in sorted(job_dir.glob("attempt_*"), reverse=True):
        try:
            state = json.loads((attempt_dir / "state.json").read_text())
            if state.get("status") in ("failed", "interrupted", "running"):
                return {"status": "interrupted" if state["status"] == "running" else state["status"],
                        "attempt_dir": str(attempt_dir), "error": state.get("error", state["status"])}
        except (OSError, ValueError):
            pass
    return {"status": "unverified"}


def write_csv(path, rows):
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as destination:
        writer = csv.DictWriter(destination, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def save_snapshot(plan, outcomes, invocation_dir, state, started, invalid_attempts):
    summaries = summarize_results(plan, outcomes)
    jobs = []
    for human in plan["humans"]:
        for query in plan["shuffles"]:
            record = outcomes.get((human["transcript_id"], query["shuffle_id"]), {"status": "missing"})
            jobs.append({"transcript_id": human["transcript_id"], "shuffle_id": query["shuffle_id"],
                         "status": record["status"], "energy": record.get("energy"),
                         "attempt_dir": record.get("attempt_dir", ""), "error": record.get("error", "")})
    write_csv(invocation_dir / "job_status.csv", jobs)
    write_csv(invocation_dir / "summary.csv", summaries)
    atomic_json(invocation_dir / "summary.json", summaries)
    report = {"status": state, "updated_utc": utc_now(), "plan_sha256": plan["_plan_sha256"],
              "requested_jobs": len(jobs), "verified_jobs": sum(row["verified_controls"] for row in summaries),
              "failed_jobs": sum(row["failed"] for row in summaries),
              "unverified_jobs": sum(row["unverified"] for row in summaries),
              "missing_jobs": sum(row["missing"] for row in summaries),
              "elapsed_seconds": round(time.monotonic() - started, 3),
              "invalid_or_incomplete_attempts": invalid_attempts, "evidence_label": "exploratory"}
    atomic_json(invocation_dir / "run_status.json", report)
    lines = ["Each transcript is summarized separately using its full RNA and the same 100 existing shuffled viral sequences.", "",
             "The empirical energy tail fraction is (1 + controls with reported E <= original E) / 101. "
             "The official fraction is available when all 100 controls have verified terminal outcomes. "
             "Failed, missing and unverified jobs are counted explicitly. A successful header-only CSV is a no-interaction outcome.", "",
             "For the primary shuffled mean and median, a verified no-interaction outcome has analysis value +infinity. "
             "This convention contributes zero to the lower-energy tail event and preserves all 100 outcomes. "
             "The reported_energy_mean and reported_energy_median columns summarize only finite reported energies, "
             "with completed_with_interaction giving their sample size. Pending summaries describe verified outcomes so far.", "",
             "| Transcript | Original E | Verified | No interaction | Failed | Missing | Shuffled mean | Shuffled median | E <= original | Tail fraction |",
             "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
    for row in summaries:
        display = lambda value: "pending" if value is None else str(value)
        lines.append(f"| {row['transcript_id']} | {row['original_energy']} | {row['verified_controls']} | "
                     f"{row['completed_no_interaction']} | {row['failed'] + row['unverified']} | {row['missing']} | "
                     f"{display(row['shuffled_mean'])} | {display(row['shuffled_median'])} | "
                     f"{display(row['controls_le_original'])} | {display(row['empirical_tail_fraction'])} |")
    lines.extend(["", "The four transcripts tagged shared_626nt share sequence and remain related observations. "
                  "Every candidate remains exploratory. These controls preserve the Step 4 multiple-testing result and do not establish biological binding."])
    (invocation_dir / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report


@contextlib.contextmanager
def bundle_lock(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a+b") as handle:
        handle.seek(0, os.SEEK_END)
        if handle.tell() == 0:
            handle.write(b"0")
            handle.flush()
        handle.seek(0)
        if os.name == "nt":
            import msvcrt
            try:
                msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
            except OSError:
                raise RuntimeError("Another control runner is using this bundle") from None
            try:
                yield
            finally:
                handle.seek(0)
                msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
        else:
            import fcntl
            try:
                fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            except OSError:
                raise RuntimeError("Another control runner is using this bundle") from None
            try:
                yield
            finally:
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--intarna", default=os.environ.get("INTARNA_BIN", "IntaRNA"))
    parser.add_argument("--summarize-only", action="store_true")
    args = parser.parse_args(argv)
    plan = load_plan(args.bundle)
    bundle = Path(plan["_bundle"])
    with bundle_lock(bundle / "results/.runner.lock"):
        invocation_dir = unique_directory(bundle / "results/invocations", "run")
        started = time.monotonic()
        outcomes = {}
        invalid_attempts = []
        executable = None
        versions = plan["expected_versions"]
        active_key = None
        active_job_dir = None
        previous_term_handler = signal.getsignal(signal.SIGTERM)
        def interrupt(signum, frame):
            raise KeyboardInterrupt
        signal.signal(signal.SIGTERM, interrupt)
        try:
            if not args.summarize_only:
                executable = shutil.which(args.intarna)
                if executable is None:
                    raise ValueError(f"IntaRNA executable unavailable: {args.intarna}")
                versions = inspect_versions(Path(executable), invocation_dir)
                if {key: versions[key] for key in plan["expected_versions"]} != plan["expected_versions"]:
                    raise ValueError("Runtime versions differ from the original prediction versions")
            for human in plan["humans"]:
                for query in plan["shuffles"]:
                    job_dir = bundle / "results/jobs" / human["transcript_id"] / f"shuffle_{query['shuffle_id']:03d}"
                    expected = job_fingerprint(plan, human, query, versions)
                    record, invalid = find_verified_attempt(job_dir, expected)
                    invalid_attempts.extend(invalid)
                    outcomes[(human["transcript_id"], query["shuffle_id"])] = record or latest_unverified_state(job_dir, invalid)
            save_snapshot(plan, outcomes, invocation_dir, "summary_only" if args.summarize_only else "running", started, invalid_attempts)
            if args.summarize_only:
                print(f"Summary snapshot saved to {invocation_dir}", flush=True)
                return 0
            total = len(plan["humans"]) * plan["num_shuffles"]
            progress = 0
            for human in plan["humans"]:
                for query in plan["shuffles"]:
                    progress += 1
                    key = (human["transcript_id"], query["shuffle_id"])
                    record = outcomes[key]
                    resumed = record["status"] in SUCCESS_STATES
                    verify_sequence_entry(bundle, human, human["transcript_id"])
                    verify_sequence_entry(bundle, query, query["id"])
                    if resumed:
                        try:
                            outcomes[key] = verified_attempt(Path(record["attempt_dir"]),
                                                             job_fingerprint(plan, human, query, versions))
                        except (OSError, ValueError, KeyError, TypeError) as error:
                            invalid_attempts.append({"attempt_dir": record["attempt_dir"],
                                                     "reason": f"{type(error).__name__}: {error}"})
                            resumed = False
                    print(f"[{progress}/{total}] {human['transcript_id']} shuffle {query['shuffle_id']:03d}/100 "
                          f"{'verified resume' if resumed else 'starting'} elapsed={time.monotonic() - started:.1f}s", flush=True)
                    if not resumed:
                        active_key = key
                        active_job_dir = bundle / "results/jobs" / human["transcript_id"] / f"shuffle_{query['shuffle_id']:03d}"
                        outcomes[key] = run_job(plan, human, query, versions, executable, invocation_dir)
                        active_key = None
                        active_job_dir = None
                    print(f"[{progress}/{total}] {human['transcript_id']} shuffle {query['shuffle_id']:03d}/100 "
                          f"{outcomes[key]['status']} job_elapsed={outcomes[key].get('elapsed_seconds', 0):.1f}s "
                          f"elapsed={time.monotonic() - started:.1f}s", flush=True)
                    save_snapshot(plan, outcomes, invocation_dir, "running", started, invalid_attempts)
            complete = all(record["status"] in SUCCESS_STATES for record in outcomes.values())
            final = save_snapshot(plan, outcomes, invocation_dir, "complete" if complete else "completed_with_failures", started, invalid_attempts)
            print(f"{final['status']}: {final['verified_jobs']}/{total} verified controls. Summary: {invocation_dir}", flush=True)
            return 0 if complete else 1
        except BaseException as error:
            # A killed attempt stays uncommitted and is retried in a fresh directory.
            state = "interrupted" if isinstance(error, KeyboardInterrupt) else "failed"
            if active_key is not None:
                outcomes[active_key] = latest_unverified_state(active_job_dir, [{"reason": "Active attempt interrupted"}])
            report = save_snapshot(plan, outcomes, invocation_dir, state, started, invalid_attempts)
            report["error"] = f"{type(error).__name__}: {error}"
            atomic_json(invocation_dir / "run_status.json", report)
            print(f"{state}: {error}. Progress saved to {invocation_dir}", file=sys.stderr, flush=True)
            if isinstance(error, KeyboardInterrupt):
                return 130
            raise
        finally:
            signal.signal(signal.SIGTERM, previous_term_handler)


if __name__ == "__main__":
    raise SystemExit(main())
