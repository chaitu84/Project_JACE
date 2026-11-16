"""Checks for resumable controls using synthetic inputs and a fake executable."""

import contextlib
import csv
import hashlib
import io
import json
import os
import shlex
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


SOURCE_DIR = Path(__file__).resolve().parents[1] / "src/config_code/predict_interactions"
sys.path.insert(0, str(SOURCE_DIR))

import intarna_controls as controls  # noqa: E402
import prepare_exploratory_controls as preparation  # noqa: E402


SETTINGS = [
    "--model=X", "--mode=H", "--acc=C", "--tAccW=150", "--tAccL=100",
    "--qAccW=150", "--qAccL=100", "--seedBP=7", "--outNumber=1", "--outMode=C",
    "--outCsvCols=id1,start1,end1,id2,start2,end2,E,ED1,ED2",
]
CSV_FIELDS = "id1,start1,end1,id2,start2,end2,E,ED1,ED2".split(",")


def file_hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sequence_hash(sequence):
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def write_prediction(path, human_id="H1.1", viral_id="V_shuffle_001", energy=-20.0, no_interaction=False):
    with path.open("w", newline="", encoding="utf-8") as target:
        writer = csv.writer(target, delimiter=";")
        writer.writerow(CSV_FIELDS)
        if not no_interaction:
            writer.writerow([human_id, 1, 10, viral_id, 2, 11, energy, 0.4, 0.8])


class ResumableControlTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name)
        self.bundle = self.directory / "control bundle"
        self.bundle.mkdir()

    def make_plan(self, human_count=1):
        (self.bundle / "human").mkdir(exist_ok=True)
        (self.bundle / "viral_shuffles").mkdir(exist_ok=True)
        (self.bundle / "original_predictions").mkdir(exist_ok=True)
        original_sequence = "UGCA" * 6
        original_fasta = self.bundle / "original_viral.rna.fasta"
        original_fasta.write_text(f">V\n{original_sequence}\n", encoding="utf-8")
        version_log = self.bundle / "original_predictions/intarna_version.log"
        version_log.write_text("IntaRNA 3.4.1\n using Vienna RNA package 2.7.0 and boost 1.85.0\n", encoding="utf-8")
        humans = []
        for number in range(1, human_count + 1):
            identifier = f"H{number}.1"
            sequence = "ACGU" * (5 + number)
            fasta = self.bundle / "human" / f"{identifier}.rna.fasta"
            fasta.write_text(f">{identifier} synthetic full transcript\n{sequence}\n", encoding="utf-8")
            original_csv = self.bundle / "original_predictions" / f"{identifier}.csv"
            write_prediction(original_csv, human_id=identifier, viral_id="V", energy=-20.0)
            original_command = original_csv.with_suffix(".command.sh")
            command = ["IntaRNA", f"--target={fasta}", f"--query={original_fasta}", *SETTINGS, f"--out={original_csv}"]
            original_command.write_text(shlex.join(command) + "\n", encoding="utf-8")
            humans.append({
                "transcript_id": identifier,
                "fasta": fasta.relative_to(self.bundle).as_posix(),
                "file_sha256": file_hash(fasta),
                "sequence_sha256": sequence_hash(sequence),
                "length_nt": len(sequence),
                "original_energy": -20.0,
                "related_group": "shared_region" if number > 1 else "distinct_region",
                "original_csv": original_csv.relative_to(self.bundle).as_posix(),
                "original_csv_sha256": file_hash(original_csv),
                "original_command": original_command.relative_to(self.bundle).as_posix(),
                "original_command_sha256": file_hash(original_command),
            })
        shuffles = []
        for number in range(1, 101):
            identifier = f"V_shuffle_{number:03d}"
            sequence = "UGCA" * 6
            fasta = self.bundle / "viral_shuffles" / f"viral_shuffle_{number:03d}.fasta"
            fasta.write_text(f">{identifier}\n{sequence}\n", encoding="utf-8")
            shuffles.append({
                "shuffle_id": number, "id": identifier,
                "fasta": fasta.relative_to(self.bundle).as_posix(),
                "file_sha256": file_hash(fasta), "sequence_sha256": sequence_hash(sequence),
                "length_nt": len(sequence),
            })
        plan = {
            "schema_version": 1,
            "expected_versions": {"intarna": "3.4.1", "viennarna": "2.7.0"},
            "settings_args": SETTINGS.copy(), "num_shuffles": 100,
            "humans": humans, "shuffles": shuffles,
            "original_viral": {
                "id": "V", "fasta": original_fasta.relative_to(self.bundle).as_posix(),
                "file_sha256": file_hash(original_fasta), "sequence_sha256": sequence_hash(original_sequence),
                "length_nt": len(original_sequence),
            },
            "original_version_log": {
                "path": version_log.relative_to(self.bundle).as_posix(),
                "file_sha256": file_hash(version_log),
            },
        }
        (self.bundle / "control_plan.json").write_text(json.dumps(plan, indent=2) + "\n", encoding="utf-8")
        return plan

    def fake_process(self, command, stdout_path, stderr_path):
        options = {argument.split("=", 1)[0]: argument.split("=", 1)[1]
                   for argument in command[1:] if "=" in argument}
        human_id = Path(options["--target"]).read_text().splitlines()[0][1:].split()[0]
        viral_id = Path(options["--query"]).read_text().splitlines()[0][1:].split()[0]
        output = Path(options["--out"])
        write_prediction(output, human_id, viral_id, energy=-20.0)
        Path(stdout_path).write_text("fake control completed\n", encoding="utf-8")
        Path(stderr_path).write_text("", encoding="utf-8")
        return 0

    def test_native_csv_validation_includes_ids_coordinates_and_finite_energies(self):
        result = self.directory / "result.csv"
        write_prediction(result)
        validated = controls.validate_result(result, "H1.1", 24, "V_shuffle_001", 24)
        self.assertEqual(validated["status"], "completed")
        self.assertEqual(validated["prediction_rows"], 1)
        self.assertEqual(validated["energy"], -20.0)

        for bad_id in ("wrong", "H1"):
            write_prediction(result, human_id=bad_id)
            with self.assertRaises(ValueError):
                controls.validate_result(result, "H1.1", 24, "V_shuffle_001", 24)
        for bad_energy in ("nan", "inf", "-inf", "unreadable"):
            write_prediction(result, energy=bad_energy)
            with self.assertRaises(ValueError):
                controls.validate_result(result, "H1.1", 24, "V_shuffle_001", 24)
        write_prediction(result)
        result.write_text(result.read_text().replace(";1;10;", ";0;10;"), encoding="utf-8")
        with self.assertRaises(ValueError):
            controls.validate_result(result, "H1.1", 24, "V_shuffle_001", 24)

    def test_header_only_is_a_verified_no_interaction_outcome(self):
        result = self.directory / "none.csv"
        write_prediction(result, no_interaction=True)
        validated = controls.validate_result(result, "H1.1", 24, "V_shuffle_001", 24)
        self.assertEqual(validated["status"], "completed_no_interaction")
        self.assertEqual(validated["prediction_rows"], 0)
        self.assertIsNone(validated["energy"])

    def test_missing_header_extra_rows_and_missing_values_are_rejected(self):
        result = self.directory / "invalid.csv"
        invalid_contents = [
            "",
            "wrong;header\n",
            ";".join(CSV_FIELDS) + "\nH1.1;1;10;V_shuffle_001;2;11;-20;0.4\n",
            ";".join(CSV_FIELDS) + "\nH1.1;1;10;V_shuffle_001;2;11;-20;0.4;0.8\n" * 2,
        ]
        for content in invalid_contents:
            with self.subTest(content=content):
                result.write_text(content, encoding="utf-8")
                with self.assertRaises(ValueError):
                    controls.validate_result(result, "H1.1", 24, "V_shuffle_001", 24)

    def test_plan_rejects_changed_human_shuffle_or_original_prediction(self):
        plan = self.make_plan()
        controls.load_plan(self.bundle)
        for key, entry, path_key in (
            ("human", plan["humans"][0], "fasta"),
            ("shuffle", plan["shuffles"][0], "fasta"),
            ("original", plan["humans"][0], "original_csv"),
            ("command", plan["humans"][0], "original_command"),
        ):
            with self.subTest(source=key):
                path = self.bundle / entry[path_key]
                before = path.read_bytes()
                path.write_bytes(before + b"\n")
                with self.assertRaises(ValueError):
                    controls.load_plan(self.bundle)
                path.write_bytes(before)

    def prepare_runtime(self):
        self.make_plan()
        plan = controls.load_plan(self.bundle)
        invocation = self.directory / "runtime invocation"
        invocation.mkdir()
        version_log = invocation / "version.log"
        version_log.write_text("IntaRNA 3.4.1\n using Vienna RNA package 2.7.0 and boost 1.85.0\n", encoding="utf-8")
        executable = self.directory / "fake runtime.exe"
        executable.write_bytes(b"synthetic executable boundary only")
        versions = {"intarna": "3.4.1", "viennarna": "2.7.0",
                    "version_log_sha256": file_hash(version_log), "executable_sha256": file_hash(executable)}
        return plan, plan["humans"][0], plan["shuffles"][0], versions, executable, invocation

    def test_completed_attempt_requires_matching_request_settings_versions_and_artifact_hashes(self):
        plan, human, query, versions, executable, invocation = self.prepare_runtime()
        with mock.patch.object(controls, "run_process", side_effect=self.fake_process):
            record = controls.run_job(plan, human, query, versions, executable, invocation)
        self.assertEqual(record["status"], "completed")
        attempt = Path(record["attempt_dir"])
        expected = controls.job_fingerprint(plan, human, query, versions)
        verified, invalid = controls.find_verified_attempt(attempt.parent, expected)
        self.assertEqual(verified["energy"], -20.0)
        self.assertEqual(invalid, [])

        for name in ("stderr.log", "request.json", "command.json", "version.log"):
            with self.subTest(artifact=name):
                artifact = attempt / name
                original = artifact.read_bytes()
                artifact.write_bytes(original + b"\n")
                verified, invalid = controls.find_verified_attempt(attempt.parent, expected)
                self.assertIsNone(verified)
                self.assertEqual(len(invalid), 1)
                artifact.write_bytes(original)
        changed_versions = dict(versions, intarna="3.4.2")
        changed = controls.job_fingerprint(plan, human, query, changed_versions)
        self.assertIsNone(controls.find_verified_attempt(attempt.parent, changed)[0])
        changed = json.loads(json.dumps(expected))
        changed["human"]["sequence_sha256"] = "0" * 64
        self.assertIsNone(controls.find_verified_attempt(attempt.parent, changed)[0])
        changed = json.loads(json.dumps(expected))
        changed["settings_args"] = [*changed["settings_args"], "--seedBP=8"]
        self.assertIsNone(controls.find_verified_attempt(attempt.parent, changed)[0])

    def test_failed_and_interrupted_attempts_are_preserved_and_retried_in_fresh_directories(self):
        plan, human, query, versions, executable, invocation = self.prepare_runtime()

        def fail_process(command, stdout_path, stderr_path):
            self.fake_process(command, stdout_path, stderr_path)
            return 9

        with mock.patch.object(controls, "run_process", side_effect=fail_process):
            failed = controls.run_job(plan, human, query, versions, executable, invocation)
        self.assertEqual(failed["status"], "failed")
        failed_attempt = Path(failed["attempt_dir"])
        self.assertTrue((failed_attempt / "result.csv.partial").exists())
        self.assertFalse((failed_attempt / "completed.json").exists())
        failed_bytes = {path.name: path.read_bytes() for path in failed_attempt.iterdir()}

        def interrupt_process(command, stdout_path, stderr_path):
            self.fake_process(command, stdout_path, stderr_path)
            raise KeyboardInterrupt("synthetic interruption")

        with mock.patch.object(controls, "run_process", side_effect=interrupt_process):
            with self.assertRaises(KeyboardInterrupt):
                controls.run_job(plan, human, query, versions, executable, invocation)
        attempts = list(failed_attempt.parent.glob("attempt_*"))
        self.assertEqual(len(attempts), 2)
        interrupted = next(path for path in attempts if path != failed_attempt)
        self.assertEqual(json.loads((interrupted / "state.json").read_text())["status"], "interrupted")
        self.assertFalse((interrupted / "completed.json").exists())

        with mock.patch.object(controls, "run_process", side_effect=self.fake_process):
            completed = controls.run_job(plan, human, query, versions, executable, invocation)
        self.assertEqual(completed["status"], "completed")
        self.assertEqual(len(list(failed_attempt.parent.glob("attempt_*"))), 3)
        self.assertEqual({path.name: path.read_bytes() for path in failed_attempt.iterdir()}, failed_bytes)

    def test_input_change_during_a_process_cannot_commit_a_completed_result(self):
        plan, human, query, versions, executable, invocation = self.prepare_runtime()
        human_path = self.bundle / human["fasta"]
        original = human_path.read_bytes()

        def mutate_input(command, stdout_path, stderr_path):
            self.fake_process(command, stdout_path, stderr_path)
            human_path.write_bytes(original + b"\n")
            return 0

        with mock.patch.object(controls, "run_process", side_effect=mutate_input):
            result = controls.run_job(plan, human, query, versions, executable, invocation)
        self.assertEqual(result["status"], "failed")
        self.assertFalse((Path(result["attempt_dir"]) / "completed.json").exists())
        human_path.write_bytes(original)

    def test_resume_rejects_inconsistent_version_metadata_even_when_result_is_valid(self):
        plan, human, query, versions, executable, invocation = self.prepare_runtime()
        with mock.patch.object(controls, "run_process", side_effect=self.fake_process):
            record = controls.run_job(plan, human, query, versions, executable, invocation)
        attempt = Path(record["attempt_dir"])
        marker = attempt / "completed.json"
        saved = json.loads(marker.read_text())
        saved["versions"]["version_log_sha256"] = "0" * 64
        marker.write_text(json.dumps(saved), encoding="utf-8")
        found, invalid = controls.find_verified_attempt(
            attempt.parent, controls.job_fingerprint(plan, human, query, versions)
        )
        self.assertIsNone(found)
        self.assertEqual(len(invalid), 1)

    def test_changed_executable_cannot_commit_under_an_older_version_probe(self):
        plan, human, query, versions, executable, invocation = self.prepare_runtime()

        def change_executable(command, stdout_path, stderr_path):
            self.fake_process(command, stdout_path, stderr_path)
            executable.write_bytes(b"changed synthetic executable")
            return 0

        with mock.patch.object(controls, "run_process", side_effect=change_executable):
            record = controls.run_job(plan, human, query, versions, executable, invocation)
        self.assertEqual(record["status"], "failed")
        self.assertFalse((Path(record["attempt_dir"]) / "completed.json").exists())

    def test_lock_rejects_another_writer_and_is_released_afterward(self):
        lock = self.directory / "results/.runner.lock"
        with controls.bundle_lock(lock):
            with self.assertRaises(RuntimeError):
                with controls.bundle_lock(lock):
                    self.fail("A second writer acquired the same bundle")
        with controls.bundle_lock(lock):
            pass

    def test_summarize_only_never_probes_or_starts_intarna(self):
        self.make_plan()
        with mock.patch.object(controls, "inspect_versions", side_effect=AssertionError("version probe forbidden")):
            with mock.patch.object(controls, "run_process", side_effect=AssertionError("predictions forbidden")):
                with contextlib.redirect_stdout(io.StringIO()):
                    code = controls.main(["--bundle", str(self.bundle), "--summarize-only"])
        self.assertEqual(code, 0)
        self.assertFalse((self.bundle / "results/jobs").exists())
        invocation = next((self.bundle / "results/invocations").iterdir())
        with (invocation / "summary.csv").open(newline="") as source:
            summary = list(csv.DictReader(source))[0]
        self.assertEqual(int(summary["verified_controls"]), 0)
        self.assertEqual(int(summary["missing"]), 100)
        self.assertEqual(summary["empirical_tail_fraction"], "")

    def test_preparer_reads_exact_settings_from_shell_escaped_original_commands(self):
        command = self.directory / "original.command.sh"
        escaped_settings = " ".join(option.replace(",", "\\,") for option in SETTINGS)
        command.write_text(
            "IntaRNA --target='human input.fasta' --query='viral input.fasta' "
            + escaped_settings + " --out='old result.csv'\n",
            encoding="utf-8",
        )
        self.assertEqual(preparation.command_settings(command), SETTINGS)
        command.write_text(command.read_text().replace("--seedBP=7", "--seedBP=8"), encoding="utf-8")
        with self.assertRaises(ValueError):
            preparation.command_settings(command)

    def test_preparer_existing_output_directory_preserves_every_byte(self):
        sentinel = self.bundle / "existing.csv"
        sentinel.write_bytes(b"preserve original output\r\n")
        with self.assertRaises(FileExistsError):
            preparation.prepare(self.directory, self.bundle)
        self.assertEqual({path.name: path.read_bytes() for path in self.bundle.iterdir()},
                         {"existing.csv": b"preserve original output\r\n"})

    def test_tail_counts_ties_and_uses_all_100_controls_per_transcript(self):
        plan = self.make_plan(human_count=5)
        outcomes = {}
        for human in plan["humans"]:
            for shuffle_id in range(1, 101):
                energy = -21.0 if shuffle_id <= 40 else -20.0 if shuffle_id <= 60 else -19.0
                outcomes[(human["transcript_id"], shuffle_id)] = {"status": "completed", "energy": energy}
        summaries = controls.summarize_results(plan, outcomes)
        self.assertEqual(len(summaries), 5)
        self.assertEqual([row["transcript_id"] for row in summaries], [f"H{number}.1" for number in range(1, 6)])
        for row in summaries:
            self.assertEqual(row["verified_controls"], 100)
            self.assertEqual(row["original_energy"], -20.0)
            self.assertAlmostEqual(float(row["shuffled_mean"]), -20.0)
            self.assertAlmostEqual(float(row["shuffled_median"]), -20.0)
            self.assertEqual(row["controls_le_original"], 60)
            self.assertAlmostEqual(row["empirical_tail_fraction"], 61 / 101)
        self.assertEqual([row["related_group"] for row in summaries[1:]], ["shared_region"] * 4)

    def test_missing_or_failed_control_keeps_official_tail_fraction_unreported(self):
        plan = self.make_plan()
        outcomes = {("H1.1", number): {"status": "completed", "energy": -21.0} for number in range(1, 100)}
        missing = controls.summarize_results(plan, outcomes)[0]
        self.assertEqual(missing["verified_controls"], 99)
        self.assertIsNone(missing["empirical_tail_fraction"])
        self.assertEqual(missing["missing"], 1)
        outcomes[("H1.1", 100)] = {"status": "failed", "energy": None}
        failed = controls.summarize_results(plan, outcomes)[0]
        self.assertEqual(failed["verified_controls"], 99)
        self.assertIsNone(failed["empirical_tail_fraction"])
        self.assertEqual(failed["failed"], 1)

    def test_no_interaction_remains_explicit_in_extended_and_reported_energy_summaries(self):
        plan = self.make_plan()
        outcomes = {("H1.1", number): {"status": "completed", "energy": -20.0} for number in range(1, 100)}
        outcomes[("H1.1", 100)] = {"status": "completed_no_interaction", "energy": None}
        summary = controls.summarize_results(plan, outcomes)[0]
        self.assertEqual(summary["verified_controls"], 100)
        self.assertEqual(summary["completed_no_interaction"], 1)
        self.assertTrue(float(summary["shuffled_mean"]) == float("inf"))
        self.assertEqual(float(summary["shuffled_median"]), -20.0)
        self.assertEqual(summary["reported_energy_mean"], -20.0)
        self.assertEqual(summary["reported_energy_median"], -20.0)
        self.assertEqual(summary["controls_le_original"], 99)
        self.assertAlmostEqual(summary["empirical_tail_fraction"], 100 / 101)

    def make_fake_executable(self):
        fake_script = self.directory / "fake_intarna.py"
        fake_script.write_text(
            """import csv
import os
import sys
from pathlib import Path

if '--version' in sys.argv:
    print('IntaRNA 3.4.1')
    print(' using Vienna RNA package 2.7.0 and boost 1.85.0')
    raise SystemExit(0)

options = dict(argument.split('=', 1) for argument in sys.argv[1:])
human_id = Path(options['--target']).read_text().splitlines()[0][1:].split()[0]
viral_id = Path(options['--query']).read_text().splitlines()[0][1:].split()[0]
shuffle_number = int(viral_id.rsplit('_', 1)[1])
events = Path(os.environ['FAKE_INTARNA_EVENTS'])
with events.open('a') as stream:
    stream.write(f'start {human_id} {shuffle_number}\\n')
try:
    with Path(options['--out']).open('w', newline='') as stream:
        writer = csv.writer(stream, delimiter=';')
        writer.writerow('id1,start1,end1,id2,start2,end2,E,ED1,ED2'.split(','))
        if shuffle_number != int(os.environ.get('FAKE_NO_INTERACTION', '0')):
            energy = -21 if shuffle_number <= 40 else -20 if shuffle_number <= 60 else -19
            writer.writerow([human_id, 1, 10, viral_id, 2, 11, energy, .4, .8])
    print('fake prediction completed')
finally:
    with events.open('a') as stream:
        stream.write(f'end {human_id} {shuffle_number}\\n')
""",
            encoding="utf-8",
        )
        if os.name == "nt":
            executable = self.directory / "fake_intarna.cmd"
            executable.write_text(
                f'@echo off\n"{sys.executable}" "%~dp0fake_intarna.py" %*\n',
                encoding="utf-8",
            )
        else:
            executable = self.directory / "fake_intarna"
            executable.write_text(f'#!/bin/sh\nexec {shlex.quote(sys.executable)} {shlex.quote(str(fake_script))} "$@"\n', encoding="utf-8")
            executable.chmod(0o755)
        return executable

    def run_fake_cli(self, executable, events):
        environment = os.environ.copy()
        environment.update(FAKE_INTARNA_EVENTS=str(events), FAKE_NO_INTERACTION="100")
        return subprocess.run(
            [sys.executable, str(SOURCE_DIR / "intarna_controls.py"), "--bundle", str(self.bundle), "--intarna", str(executable)],
            capture_output=True, text=True, env=environment, timeout=120,
        )

    def test_real_fake_cli_runs_sequentially_resumes_and_retries_only_tampered_result(self):
        self.make_plan()
        executable = self.make_fake_executable()
        events = self.directory / "fake_events.log"
        first = self.run_fake_cli(executable, events)
        self.assertEqual(first.returncode, 0, first.stdout + first.stderr)
        calls = events.read_text().splitlines()
        self.assertEqual(len(calls), 200)
        for number in range(100):
            self.assertEqual(calls[number * 2], f"start H1.1 {number + 1}")
            self.assertEqual(calls[number * 2 + 1], f"end H1.1 {number + 1}")

        completed = sorted((self.bundle / "results/jobs").rglob("completed.json"))
        self.assertEqual(len(completed), 100)
        snapshots = {path: path.read_bytes() for path in completed}
        second = self.run_fake_cli(executable, events)
        self.assertEqual(second.returncode, 0, second.stdout + second.stderr)
        self.assertEqual(events.read_text().splitlines(), calls)
        self.assertEqual({path: path.read_bytes() for path in completed}, snapshots)

        result = next(path.parent / "result.csv" for path in completed if "shuffle_001" in path.parts)
        result.write_text(result.read_text().replace(";-21;", ";-22;"), encoding="utf-8")
        tampered_bytes = result.read_bytes()
        third = self.run_fake_cli(executable, events)
        self.assertEqual(third.returncode, 0, third.stdout + third.stderr)
        self.assertEqual(events.read_text().splitlines()[-2:], ["start H1.1 1", "end H1.1 1"])
        self.assertEqual(len(events.read_text().splitlines()), 202)
        self.assertEqual(result.read_bytes(), tampered_bytes)
        self.assertEqual(len(list((self.bundle / "results/jobs").rglob("completed.json"))), 101)
        invocation_dirs = list((self.bundle / "results/invocations").iterdir())
        self.assertEqual(len(invocation_dirs), 3)
        for invocation in invocation_dirs:
            with (invocation / "summary.csv").open(newline="") as source:
                summary = list(csv.DictReader(source))[0]
            self.assertEqual(int(summary["verified_controls"]), 100)
            self.assertEqual(int(summary["completed_no_interaction"]), 1)
            self.assertEqual(int(summary["controls_le_original"]), 60)
            self.assertAlmostEqual(float(summary["empirical_tail_fraction"]), 61 / 101)


if __name__ == "__main__":
    unittest.main()
