"""Prepare existing viral shuffles and a resumable runner without launching it."""

import argparse
import csv
import hashlib
import json
import math
import re
import shlex
import shutil
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from Bio import SeqIO


PROJECT_PATH = Path(__file__).resolve().parents[3]
ORIGINAL_RUN = "run_20261116T044527Z_4D0EGI"
CSV_COLUMNS = "id1,start1,end1,id2,start2,end2,E,ED1,ED2".split(",")
EXPECTED_SETTINGS = [
    "--model=X", "--mode=H", "--acc=C", "--tAccW=150", "--tAccL=100",
    "--qAccW=150", "--qAccL=100", "--seedBP=7", "--outNumber=1", "--outMode=C",
    "--outCsvCols=id1,start1,end1,id2,start2,end2,E,ED1,ED2",
]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fasta_metadata(path):
    record = SeqIO.read(path, "fasta")
    sequence = str(record.seq)
    require(bool(sequence) and set(sequence) <= set("ACGU"), f"Expected original unambiguous RNA in {path}")
    return {"id": record.id, "file_sha256": file_sha256(path),
            "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
            "length_nt": len(sequence)}, sequence


def command_settings(path):
    arguments = shlex.split(Path(path).read_text(encoding="utf-8"))
    require(bool(arguments), f"Empty command in {path}")
    settings = [argument for argument in arguments[1:]
                if not argument.startswith(("--target=", "--query=", "--out="))]
    require(settings == EXPECTED_SETTINGS, f"Original prediction settings differ in {path}")
    return settings


def prepare(project, output_dir, original_run=None):
    project = Path(project).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=False)
    bundle = project / "results/step5/exploratory_shortlist"
    original_run = Path(original_run).resolve() if original_run else bundle / "predictions" / ORIGINAL_RUN
    old_controls = project / "results/step5/shuffled_controls"
    preparation = {"status": "preparing", "predictions_launched": False,
                   "started_utc": datetime.now(timezone.utc).isoformat(), "evidence_label": "exploratory"}
    status_path = output_dir / "preparation_status.json"
    status_path.write_text(json.dumps(preparation, indent=2) + "\n", encoding="utf-8")
    source_files = {}

    def remember(path):
        path = Path(path)
        source_files[str(path)] = file_sha256(path)

    def copy_input(source, relative):
        source = Path(source)
        remember(source)
        destination = output_dir / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        require(not destination.exists(), f"Prepared file already exists: {destination}")
        shutil.copyfile(source, destination)
        require(file_sha256(destination) == source_files[str(source)], f"Copy differs for {source}")

    try:
        original_status = json.loads((original_run / "run_status.json").read_text(encoding="utf-8"))
        require(original_status["status"] == "complete" and original_status["completed_with_interaction"] == 5
                and original_status["failed_jobs"] == 0, "Original five predictions must be complete")
        prepared = json.loads((bundle / "preparation_manifest.json").read_text(encoding="utf-8"))
        require(prepared["status"] == "prepared", "Exploratory inputs are incomplete")
        for relative, expected in prepared["output_sha256"].items():
            require(file_sha256(bundle / relative) == expected, f"Prepared exploratory input changed: {relative}")
        require((original_run / "prediction_inputs.tsv").read_bytes() == (bundle / "prediction_inputs.tsv").read_bytes(),
                "Original predictions used a different input manifest")
        with (bundle / "candidate_manifest.csv").open(newline="", encoding="utf-8") as source:
            candidates = list(csv.DictReader(source))
        require(len(candidates) == 5 and len({row["transcript_id"] for row in candidates}) == 5,
                "Expected five distinct versioned transcript IDs")
        require(all(row["evidence_label"] == "exploratory" and float(row["bh_q_value"]) > 0.05 for row in candidates),
                "Candidate labels or Step 4 BH results differ")
        blocks = json.loads((bundle / "shared_sequence_blocks.json").read_text(encoding="utf-8"))
        shared = next(block for block in blocks if block["length_nt"] == 626 and len(block["transcript_intervals0"]) == 4)
        related_ids = set(shared["transcript_intervals0"])
        require(related_ids <= {row["transcript_id"] for row in candidates}, "Unknown related candidate")

        version = (original_run / "intarna_version.log").read_text(encoding="utf-8")
        require(re.search(r"IntaRNA\s+3\.4\.1\b", version) and re.search(r"Vienna\s+RNA\s+package\s+2\.7\.0\b", version),
                "Original IntaRNA/ViennaRNA versions differ")
        original_viral_name = "HPV16_E6starI_E7_spliced_region.original.rna.fasta"
        viral_metadata, original_viral = fasta_metadata(bundle / original_viral_name)
        copy_input(bundle / original_viral_name, f"original_predictions/{original_viral_name}")
        viral_metadata["fasta"] = f"original_predictions/{original_viral_name}"
        humans = []
        for candidate in candidates:
            identifier = candidate["transcript_id"]
            require(re.fullmatch(r"ENST\d+\.\d+", identifier), f"Expected a versioned transcript ID: {identifier}")
            human_relative = f"human/{identifier}.full.rna.fasta"
            human, sequence = fasta_metadata(bundle / candidate["human_fasta"])
            require(human["id"] == identifier and human["sequence_sha256"] == candidate["sequence_sha256"]
                    and len(sequence) == int(candidate["length_nt"]), f"Candidate identity differs for {identifier}")
            copy_input(bundle / candidate["human_fasta"], human_relative)
            csv_path = original_run / f"{identifier}.csv"
            with csv_path.open(newline="", encoding="utf-8-sig") as source:
                reader = csv.DictReader(source, delimiter=";")
                require(reader.fieldnames == CSV_COLUMNS, f"Wrong original CSV header in {csv_path}")
                predictions = list(reader)
            require(len(predictions) == 1, f"Expected one original prediction for {identifier}")
            prediction = predictions[0]
            require(prediction["id1"] == identifier and prediction["id2"] == viral_metadata["id"], "Original prediction IDs differ")
            require(all(math.isfinite(float(prediction[key])) for key in ("E", "ED1", "ED2")), "Nonfinite original energy")
            require(1 <= int(prediction["start1"]) <= int(prediction["end1"]) <= len(sequence)
                    and 1 <= int(prediction["start2"]) <= int(prediction["end2"]) <= len(original_viral),
                    "Original prediction coordinates differ")
            command_path = original_run / f"{identifier}.command.sh"
            command_settings(command_path)
            for suffix in ("csv", "command.sh", "stdout.log", "stderr.log", "validation.log"):
                copy_input(original_run / f"{identifier}.{suffix}", f"original_predictions/{identifier}.{suffix}")
            humans.append({
                "transcript_id": identifier, "fasta": human_relative,
                "file_sha256": human["file_sha256"], "sequence_sha256": human["sequence_sha256"],
                "length_nt": human["length_nt"], "original_energy": float(prediction["E"]),
                "original_csv": f"original_predictions/{identifier}.csv", "original_csv_sha256": file_sha256(csv_path),
                "original_command": f"original_predictions/{identifier}.command.sh",
                "original_command_sha256": file_sha256(command_path),
                "related_group": "shared_626nt" if identifier in related_ids else "LINC00869",
                "gene_id": candidate["gene_id"], "step4_bh_q_value": float(candidate["bh_q_value"]),
            })
        require(len({human["sequence_sha256"] for human in humans}) == 5, "Expected five distinct full human sequences")

        # Compare existing saved sequences. No shuffling or normalization is performed.
        batches = []
        for name in ("viral_shuffles_batch_1.fasta", "viral_shuffles_batch_2.fasta"):
            path = old_controls / name
            remember(path)
            batches.extend((record.id, str(record.seq)) for record in SeqIO.parse(path, "fasta"))
        step4_path = project / "results/step4/viral_shuffles.fasta"
        remember(step4_path)
        step4 = [(record.id, str(record.seq)) for record in SeqIO.parse(step4_path, "fasta")]
        require(len(batches) == len(step4) == 100 and batches == step4, "Existing batch and Step 4 shuffles differ")
        shuffles = []
        for number in range(1, 101):
            source = old_controls / f"viral_shuffle_{number:03d}.fasta"
            metadata, sequence = fasta_metadata(source)
            require((metadata["id"], sequence) == batches[number - 1], f"Saved shuffle {number} differs from the existing pool")
            require(metadata["id"] == f"{viral_metadata['id']}_shuffle_{number:03d}", f"Shuffle ID differs for {number}")
            require(Counter(sequence) == Counter(original_viral), f"Viral base composition differs for shuffle {number}")
            relative = f"viral_shuffles/{source.name}"
            copy_input(source, relative)
            shuffles.append({"shuffle_id": number, "fasta": relative, **metadata})
        require(len({shuffle["sequence_sha256"] for shuffle in shuffles}) == 100, "Existing shuffle pool contains duplicate sequences")
        for name in ("run_status.json", "job_status.tsv", "intarna_version.log", "prediction_inputs.tsv"):
            copy_input(original_run / name, f"original_predictions/{name}")
        for name in ("preparation_manifest.json", "candidate_manifest.csv", "shared_sequence_blocks.json"):
            copy_input(bundle / name, f"provenance/{name}")
        remember(old_controls / "intarna_100_controls.zip")
        engine = Path(__file__).with_name("intarna_controls.py")
        copy_input(engine, "intarna_controls.py")
        plan = {
            "schema_version": 1, "evidence_label": "exploratory", "num_shuffles": 100,
            "expected_versions": {"intarna": "3.4.1", "viennarna": "2.7.0"},
            "settings_args": EXPECTED_SETTINGS, "humans": humans, "shuffles": shuffles,
            "original_viral": viral_metadata,
            "original_version_log": {"path": "original_predictions/intarna_version.log",
                                     "file_sha256": file_sha256(original_run / "intarna_version.log")},
            "original_prediction_run": str(original_run), "source_files": source_files,
            "related_candidates": {"shared_626nt": sorted(related_ids)},
            "no_interaction_analysis_convention": "Verified successful header-only output is a +infinity outcome for the energy tail and extended-distribution statistics; reported-only statistics are also provided",
            "statistical_context": "All five are exploratory. None passed Step 4 BH FDR 0.05. These controls do not establish biological binding.",
        }
        (output_dir / "control_plan.json").write_text(json.dumps(plan, indent=2) + "\n", encoding="utf-8")
        launcher = '''#!/usr/bin/env bash
set -euo pipefail
bundle_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
exec "${PYTHON_BIN:-python3}" "$bundle_dir/intarna_controls.py" --bundle "$bundle_dir" "$@"
'''
        (output_dir / "run_controls.sh").write_text(launcher, encoding="utf-8", newline="\n")
        ubuntu_launcher = (output_dir / "run_controls.sh").as_posix()
        if output_dir.drive:
            ubuntu_launcher = f"/mnt/{output_dir.drive[0].lower()}{ubuntu_launcher[2:]}"
        instructions = [
            "Prepared 500 sequential exploratory control jobs, with 100 saved viral shuffles for each of five full human RNAs.", "",
            "The 100 individual shuffled FASTAs were copied byte for byte and verified against the existing two batches and Step 4 pool. "
            "Human FASTAs, original prediction CSVs, commands and logs were also copied unchanged. The earlier completed control results remain preserved.", "",
            "The runner requires IntaRNA 3.4.1 and ViennaRNA 2.7.0. It uses exactly the saved local arguments from all five original predictions. "
            "Boost version and the complete version output are recorded. One IntaRNA pair runs at a time.", "",
            "Run or resume with the same Ubuntu command, using the standard WSL drive mount.", "", "```bash",
            "conda activate binf-intarna", f"bash {ubuntu_launcher}", "```", "",
            "Inspect a saved-results summary without launching predictions.", "", "```bash",
            f"bash {ubuntu_launcher} --summarize-only", "```", "",
            "Every attempt has its own folder under results/jobs/<transcript>/shuffle_<number>/. "
            "CSV output, stdout, stderr, command, input hashes and versions are saved immediately. "
            "Each invocation has a separate folder under results/invocations/ with progress, run status and per-transcript summaries. "
            "Successful resumption requires a completed receipt that agrees with the current input sequences, settings, versions and saved file hashes. "
            "Failed, interrupted, malformed and unverified attempts are retained. A retry creates a new attempt folder. "
            "An OS lock restricts the bundle to one active invocation.", "",
            "A zero-exit-code CSV with the required header and no data row is recorded as completed_no_interaction. "
            "An empty file, malformed result, missing file or nonzero exit is recorded as failed. "
            "A no-interaction outcome contributes zero to count(energy <= original_energy). "
            "For extended-distribution means and medians, it is represented analytically as +infinity. "
            "The summary also gives the mean and median of reported finite energies with their reported-energy count. "
            "The convention assigns no fabricated finite energy to a no-interaction outcome.", "",
            "Each transcript has a separate original energy, shuffled mean/median, control count at or below its original energy, "
            "and exploratory empirical tail fraction (count + 1) / 101. The final fraction requires 100 verified successful outcomes, "
            "including completed_no_interaction cases. Failed and missing jobs have explicit counts and keep the final fraction pending.", "",
            "ENST00000732616.1, ENST00000780340.1, ENST00000732615.1 and ENST00000754042.1 share a 626-nt sequence region. "
            "They retain separate full-transcript controls and carry the related_group value shared_626nt. "
            "These controls do not undo the Step 4 multiple-testing result or establish biological binding. All five remain exploratory.",
        ]
        (output_dir / "README.md").write_text("\n".join(instructions) + "\n", encoding="utf-8")
        for source, expected in source_files.items():
            require(file_sha256(source) == expected, f"Source changed during preparation: {source}")
        preparation.update({"status": "prepared", "finished_utc": datetime.now(timezone.utc).isoformat(),
                            "candidate_count": 5, "shuffles_per_candidate": 100, "planned_jobs": 500,
                            "existing_shuffle_bytes_preserved": True, "all_source_hashes_preserved": True,
                            "control_plan_sha256": file_sha256(output_dir / "control_plan.json"),
                            "runner_sha256": file_sha256(output_dir / "intarna_controls.py"),
                            "ubuntu_command": f"conda activate binf-intarna\nbash {ubuntu_launcher}"})
        status_path.write_text(json.dumps(preparation, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(preparation, indent=2))
        return preparation
    except BaseException as error:
        preparation.update({"status": "failed", "error": f"{type(error).__name__}: {error}"})
        status_path.write_text(json.dumps(preparation, indent=2) + "\n", encoding="utf-8")
        raise


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, default=PROJECT_PATH)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--original-run", type=Path)
    args = parser.parse_args(argv)
    output = args.output_dir or args.project / "results/step5/exploratory_controls"
    return prepare(args.project, output, args.original_run)


if __name__ == "__main__":
    main()
