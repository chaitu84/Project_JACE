"""Verify saved scientific outputs and prepare compact, local website assets.

This script reads scientific inputs and writes only src/components_UI/data/.
It never invokes IntaRNA or changes source results.
"""

import csv
import hashlib
import json
import math
import re
import shlex
import statistics
from datetime import datetime, timezone
from decimal import Decimal
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).resolve().parent / "data"
CONTROL_BUNDLE = ROOT / "results/step5/exploratory_controls"
CANDIDATE_BUNDLE = ROOT / "results/step5/exploratory_shortlist"
CONTROL_RUN = CONTROL_BUNDLE / "results/invocations/run_20261116T051139789910Z_0558ca3f"
ORIGINAL_RUN = CANDIDATE_BUNDLE / "predictions/run_20261116T044527Z_4D0EGI"
CSV_FIELDS = "id1,start1,end1,id2,start2,end2,E,ED1,ED2".split(",")
EXPECTED = {
    "ENST00000778197.1": (-21.23, -15.4354, 8),
    "ENST00000732616.1": (-25.83, -20.1227, 16),
    "ENST00000780340.1": (-25.83, -20.1231, 14),
    "ENST00000732615.1": (-25.83, -20.0387, 14),
    "ENST00000754042.1": (-25.83, -20.0384, 14),
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path, delimiter=","):
    with Path(path).open(newline="", encoding="utf-8-sig") as source:
        return list(csv.DictReader(source, delimiter=delimiter))


def read_fasta(path):
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    require(lines and lines[0].startswith(">") and sum(line.startswith(">") for line in lines) == 1,
            f"Expected one FASTA record in {path}")
    sequence = "".join(line.strip() for line in lines[1:])
    require(sequence and set(sequence) <= set("ACGU"), f"Invalid RNA in {path}")
    return lines[0][1:].split()[0], sequence


def versions(text):
    intarna = re.search(r"IntaRNA\s+(\d+\.\d+\.\d+)", text)
    vienna = re.search(r"Vienna\s*RNA\s+package\s+(\d+\.\d+\.\d+)", text)
    boost = re.search(r"boost\s+(\d+\.\d+\.\d+)", text)
    require(intarna and vienna, "Version log omits IntaRNA or ViennaRNA")
    return {"intarna": intarna.group(1), "viennarna": vienna.group(1),
            "boost": boost.group(1) if boost else None}


def parse_prediction(path, human, query):
    with path.open(newline="", encoding="utf-8-sig") as source:
        reader = csv.DictReader(source, delimiter=";")
        require(reader.fieldnames == CSV_FIELDS, f"Wrong prediction columns in {path}")
        rows = list(reader)
    require(len(rows) == 1, f"Expected one reported interaction in {path}")
    row = rows[0]
    require(row["id1"] == human["transcript_id"] and row["id2"] == query["id"], f"Wrong prediction IDs in {path}")
    require(1 <= int(row["start1"]) <= int(row["end1"]) <= human["length_nt"], f"Human coordinates outside {path}")
    require(1 <= int(row["start2"]) <= int(row["end2"]) <= query["length_nt"], f"Viral coordinates outside {path}")
    require(all(math.isfinite(float(row[key])) for key in ("E", "ED1", "ED2")), f"Nonfinite energy in {path}")
    return row


def resolve_attempt(remote, transcript, shuffle):
    normalized = remote.replace("\\", "/").rstrip("/")
    name = normalized.rsplit("/", 1)[-1]
    require(re.fullmatch(r"attempt_[A-Za-z0-9_.-]+", name), "Invalid saved attempt name")
    relative = f"results/jobs/{transcript}/shuffle_{shuffle:03d}/{name}"
    require(normalized.endswith("/results/step5/exploratory_controls/" + relative), "Saved attempt path disagrees with its job")
    return CONTROL_BUNDLE / relative


def same_number(observed, expected, label):
    require(math.isclose(float(observed), float(expected), rel_tol=0, abs_tol=1e-11), f"Numeric discrepancy in {label}")


def write_json(path, data):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def write_csv(path, rows, columns):
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as destination:
        writer = csv.DictWriter(destination, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    protected = {}
    primary_sources = []

    def protect(path, role=None):
        path = Path(path).resolve()
        require(path.is_relative_to(ROOT), f"Source outside the project: {path}")
        if path not in protected:
            protected[path] = sha256(path)
        if role:
            primary_sources.append({"path": path.relative_to(ROOT).as_posix(), "sha256": protected[path], "role": role})
        return path

    try:
        plan_path = protect(CONTROL_BUNDLE / "control_plan.json", "Verified input and settings plan")
        plan = json.loads(plan_path.read_text(encoding="utf-8"))
        require(plan["num_shuffles"] == 100 and len(plan["humans"]) == 5, "Wrong control plan dimensions")
        plan_hash = protected[plan_path]
        state_path = protect(CONTROL_RUN / "run_status.json", "Completed control-run status")
        state = json.loads(state_path.read_text(encoding="utf-8"))
        require(state["status"] == "complete" and state["requested_jobs"] == state["verified_jobs"] == 500, "Control run is incomplete")
        require(all(state[key] == 0 for key in ("failed_jobs", "missing_jobs", "unverified_jobs")), "Control run contains unresolved jobs")
        require(state["plan_sha256"] == plan_hash, "Control run plan checksum differs")
        software_path = protect(CONTROL_RUN / "version.stdout.log", "Software version output")
        software = versions(software_path.read_text())
        require({key: software[key] for key in ("intarna", "viennarna")} == plan["expected_versions"], "Software versions differ")
        protect(CONTROL_RUN / "summary.csv", "Saved per-transcript energy summary")
        protect(CONTROL_RUN / "summary.json", "Saved machine-readable energy summary")
        protect(CONTROL_RUN / "summary.md", "Saved human-readable energy summary")
        saved_summaries = {row["transcript_id"]: row for row in read_csv(CONTROL_RUN / "summary.csv")}
        saved_json = {row["transcript_id"]: row for row in json.loads((CONTROL_RUN / "summary.json").read_text())}
        jobs_path = protect(CONTROL_RUN / "job_status.csv", "All 500 actual control energies and attempts")
        jobs = read_csv(jobs_path)
        ids = [human["transcript_id"] for human in plan["humans"]]
        require(set(ids) == set(EXPECTED) == set(saved_summaries) == set(saved_json), "Candidate sets differ")
        pairs = {(row["transcript_id"], int(row["shuffle_id"])) for row in jobs}
        expected_pairs = {(identifier, number) for identifier in ids for number in range(1, 101)}
        require(len(jobs) == len(pairs) == 500 and pairs == expected_pairs, "Transcript/shuffle pairs are duplicated or incomplete")
        require(all(row["status"] == "completed" and not row["error"] for row in jobs), "Job table contains a non-interaction or failed outcome")
        jobs_by_pair = {(row["transcript_id"], int(row["shuffle_id"])): row for row in jobs}

        for entry in [plan["original_viral"], *plan["humans"], *plan["shuffles"]]:
            path = protect(CONTROL_BUNDLE / entry["fasta"])
            require(protected[path] == entry["file_sha256"], f"Input file checksum differs: {path}")
            identifier, sequence = read_fasta(path)
            expected_id = entry.get("transcript_id", entry.get("id"))
            require(identifier == expected_id and len(sequence) == entry["length_nt"], f"Input identity differs: {path}")
            require(hashlib.sha256(sequence.encode("ascii")).hexdigest() == entry["sequence_sha256"], f"Input sequence checksum differs: {path}")

        details_path = protect(CANDIDATE_BUNDLE / "candidate_manifest.csv", "GENCODE annotations and exact-match geometry")
        details = {row["transcript_id"]: row for row in read_csv(details_path)}
        shortlist_path = protect(ROOT / "results/step4/transcript_comparison/exploratory_shortlist.csv", "Step 4 exploratory candidate statistics")
        shortlist = {row["seq_id"]: row for row in read_csv(shortlist_path)}
        comparison_path = protect(ROOT / "results/step4/transcript_comparison/per_transcript_comparison.csv", "Complete Step 4 testing family")
        with comparison_path.open(newline="", encoding="utf-8") as source:
            complete_rows = {row["seq_id"]: row for row in csv.DictReader(source) if row["seq_id"] in ids}
        require(set(complete_rows) == set(ids), "Candidates missing from full Step 4 comparison")
        for identifier in ids:
            require(all(complete_rows[identifier][key] == shortlist[identifier][key] for key in shortlist[identifier]
                        if key not in ("selection_rank", "evidence_label")), "Step 4 shortlist differs from full comparison")
        step4_path = protect(ROOT / "results/step4/transcript_comparison/analysis_summary.json", "Step 4 full-family BH results")
        step4_source = json.loads(step4_path.read_text())
        regions_path = protect(CANDIDATE_BUNDLE / "match_regions.csv", "Exact-match overlap regions")
        regions = read_csv(regions_path)
        prediction_summary_path = protect(ORIGINAL_RUN / "prediction_summary.csv", "Original full-sequence prediction coordinates")
        originals = {row["transcript_id"]: row for row in read_csv(prediction_summary_path)}
        shared_path = protect(CANDIDATE_BUNDLE / "shared_sequence_blocks.json", "Verified shared sequence blocks")
        shared_blocks = json.loads(shared_path.read_text())
        shared_626 = next(block for block in shared_blocks if block["length_nt"] == 626 and len(block["transcript_intervals0"]) == 4)
        candidates, control_rows, comparison_rows = [], [], []
        for human in plan["humans"]:
            identifier = human["transcript_id"]
            detail, screen = details[identifier], shortlist[identifier]
            original_path = protect(CONTROL_BUNDLE / human["original_csv"])
            require(protected[original_path] == human["original_csv_sha256"], "Original prediction file changed")
            prediction = parse_prediction(original_path, human, plan["original_viral"])
            native_original = protect(ORIGINAL_RUN / f"{identifier}.csv", f"Original full-sequence prediction for {identifier}")
            require(protected[native_original] == protected[original_path], "Native original prediction differs from its verified copy")
            require(parse_prediction(native_original, human, plan["original_viral"]) == prediction, "Original prediction fields differ")
            original_energy = Decimal(prediction["E"])
            same_number(original_energy, human["original_energy"], f"{identifier} original energy")
            require(protected[original_path] == originals[identifier]["prediction_csv_sha256"], "Original prediction summary checksum differs")
            original_command = protect(CONTROL_BUNDLE / human["original_command"])
            require(protected[original_command] == human["original_command_sha256"], "Original command checksum differs")
            original_arguments = shlex.split(original_command.read_text())[1:]
            original_settings = [argument for argument in original_arguments if not argument.startswith(("--target=", "--query=", "--out="))]
            require(original_settings == plan["settings_args"], "Original prediction settings differ from the controls")
            energies = []
            for query in plan["shuffles"]:
                number = query["shuffle_id"]
                job = jobs_by_pair[(identifier, number)]
                attempt = resolve_attempt(job["attempt_dir"], identifier, number)
                marker = protect(attempt / "completed.json")
                receipt = json.loads(marker.read_text())
                require(receipt["status"] == "completed" and receipt["return_code"] == 0 and receipt["prediction_rows"] == 1,
                        f"Control receipt is incomplete for {identifier}/{number}")
                expected_fingerprint = {
                    "plan_sha256": plan_hash, "settings_args": plan["settings_args"], "versions": plan["expected_versions"],
                    "human": {key: human[key] for key in ("transcript_id", "fasta", "file_sha256", "sequence_sha256", "length_nt")},
                    "shuffle": {key: query[key] for key in ("shuffle_id", "id", "fasta", "file_sha256", "sequence_sha256", "length_nt")},
                }
                require(receipt["fingerprint"] == expected_fingerprint, f"Control fingerprint differs for {identifier}/{number}")
                required_artifacts = {"result.csv", "stdout.log", "stderr.log", "command.json", "command.sh", "request.json", "version.log"}
                require(set(receipt["artifact_sha256"]) == required_artifacts, "Control provenance is incomplete")
                for name, expected_hash in receipt["artifact_sha256"].items():
                    path = protect(attempt / name)
                    require(protected[path] == expected_hash, f"Control artifact checksum differs: {path}")
                require(json.loads((attempt / "request.json").read_text()) == expected_fingerprint, "Saved control request differs")
                command = json.loads((attempt / "command.json").read_text())
                require(command == receipt["command"], "Saved control command differs")
                require(shlex.split((attempt / "command.sh").read_text()) == command, "Saved shell command differs from its structured arguments")
                arguments = [argument for argument in command[1:] if not argument.startswith(("--target=", "--query=", "--out="))]
                require(arguments == plan["settings_args"], "Control energy/output settings differ")
                observed_versions = versions((attempt / "version.log").read_text())
                require({key: observed_versions[key] for key in ("intarna", "viennarna")} == plan["expected_versions"], "Control software versions differ")
                require(receipt["versions"]["version_log_sha256"] == protected[(attempt / "version.log").resolve()], "Receipt version-log checksum differs")
                result = parse_prediction(attempt / "result.csv", human, query)
                energy = Decimal(result["E"])
                same_number(energy, job["energy"], "Job-table energy")
                same_number(energy, receipt["energy"], "Receipt energy")
                require(result == receipt["fields"], "Receipt/native CSV fields differ")
                energies.append(energy)
                control_rows.append({"transcript_id": identifier, "shuffle_id": number, "energy_kcal_mol": str(energy), "status": "completed"})
            mean = sum(energies) / Decimal(100)
            median = statistics.median(energies)
            count = sum(energy <= original_energy for energy in energies)
            tail = Decimal(count + 1) / Decimal(101)
            for saved in (saved_summaries[identifier], saved_json[identifier]):
                require(saved["status"] == "complete" and int(saved["verified_controls"]) == int(saved["completed_with_interaction"]) == 100,
                        "Saved per-candidate summary is incomplete")
                require(all(int(saved[key]) == 0 for key in ("completed_no_interaction", "failed", "unverified", "missing")), "Saved summary has missing outcomes")
                for key, value in (("original_energy", original_energy), ("shuffled_mean", mean), ("shuffled_median", median),
                                   ("controls_le_original", count), ("empirical_tail_fraction", tail)):
                    same_number(saved[key], value, f"{identifier} saved {key}")
            expected_original, expected_mean, expected_count = EXPECTED[identifier]
            same_number(original_energy, expected_original, f"{identifier} requested original energy")
            same_number(mean, expected_mean, f"{identifier} requested shuffled mean")
            require(count == expected_count, f"Requested lower-energy control count differs for {identifier}")
            candidate = {
                "id": identifier, "gene_id": detail["gene_id"], "gene_name": detail["gene_name"],
                "rank": int(detail["selection_rank"]), "length_nt": int(detail["length_nt"]),
                "original_energy": float(original_energy), "shuffled_mean": float(mean), "shuffled_median": float(median),
                "controls_le_original": count, "tail_fraction": float(tail), "energies": [float(energy) for energy in energies],
                "predicted_human": [int(prediction["start1"]), int(prediction["end1"])],
                "predicted_viral": [int(prediction["start2"]), int(prediction["end2"])],
                "ED1": float(prediction["ED1"]), "ED2": float(prediction["ED2"]),
                "predicted_human_sequence": originals[identifier]["human_predicted_span_sequence"],
                "predicted_viral_sequence": originals[identifier]["viral_predicted_span_sequence"],
                "step4": {"real_count": int(screen["real_count"]), "shuffled_mean": float(screen["shuffled_mean"]),
                          "shuffled_sd": float(screen["shuffled_sd"]), "shuffles_ge_real": int(screen["shuffles_ge_real"]),
                          "empirical_p": float(screen["empirical_p"]), "bh_q_value": float(screen["bh_q_value"])},
                "related_group": human["related_group"], "overlap_clusters": int(detail["overlap_cluster_count"]),
                "max_nonoverlapping_sites": int(detail["maximum_nonoverlapping_human_sites"]),
                "dominant_repeat": {"length_nt": int(detail["dominant_overlap_cluster_nt"]), "match_rows": int(detail["dominant_overlap_cluster_rows"])},
                "match_regions": [{"start0": int(region["start0"]), "end0": int(region["end0"]), "match_rows": int(region["match_rows"]),
                                   "sequence": region["rna_sequence"]} for region in regions if region["transcript_id"] == identifier],
            }
            candidates.append(candidate)
            comparison_rows.append({
                "transcript_id": identifier, "gene_id": candidate["gene_id"], "gene_name": candidate["gene_name"],
                "length_nt": candidate["length_nt"], "original_energy_kcal_mol": float(original_energy),
                "shuffled_mean_kcal_mol": float(mean), "shuffled_median_kcal_mol": float(median),
                "controls_energy_le_original": count, "total_controls": 100, "exploratory_energy_tail_fraction": float(tail),
                "step4_empirical_p": candidate["step4"]["empirical_p"], "step4_bh_q_value": candidate["step4"]["bh_q_value"],
                "related_group": candidate["related_group"],
            })
        candidates.sort(key=lambda candidate: candidate["rank"])
        all_values = [value for candidate in candidates for value in candidate["energies"]] + [candidate["original_energy"] for candidate in candidates]
        bin_min, bin_max = math.floor(min(all_values) / 2) * 2, math.ceil(max(all_values) / 2) * 2
        bin_edges = list(range(bin_min, bin_max + 1, 2))
        require(len(bin_edges) >= 2, "Histogram bounds are invalid")
        histogram_max = 0
        for candidate in candidates:
            counts = [0] * (len(bin_edges) - 1)
            for energy in candidate["energies"]:
                index = min(int((energy - bin_min) // 2), len(counts) - 1)
                counts[index] += 1
            require(sum(counts) == 100, "Histogram lost a control energy")
            candidate["histogram_counts"] = counts
            histogram_max = max(histogram_max, max(counts))
        for path, expected_hash in protected.items():
            require(sha256(path) == expected_hash, f"Scientific source changed during preparation: {path}")
        index_text = "\n".join(f"{path.relative_to(ROOT).as_posix()} {digest}" for path, digest in sorted(protected.items()))
        data = {
            "metadata": {
                "project": "Human lncRNA–HPV16 complementarity", "evidence_label": "Exploratory computational analysis",
                "verified": True, "prepared_utc": datetime.now(timezone.utc).isoformat(),
                "candidate_count": 5, "control_count": 500, "controls_per_candidate": 100,
                "failed_controls": 0, "missing_controls": 0, "no_interaction_controls": 0,
                "control_run": CONTROL_RUN.relative_to(ROOT).as_posix(), "original_prediction_run": ORIGINAL_RUN.relative_to(ROOT).as_posix(),
                "software": software, "target": {"id": plan["original_viral"]["id"], "length_nt": plan["original_viral"]["length_nt"],
                                                    "label": "HPV16 E6*I/E7 spliced RNA region"},
                "settings_args": plan["settings_args"], "common_bin_edges": bin_edges,
                "common_histogram_max_count": histogram_max, "sources": primary_sources,
            },
            "definitions": {
                "energy_unit": "kcal/mol", "energy_direction": "Lower, more negative reported energies are more favorable within the model.",
                "energy_components": "Reported energy includes hybridization and accessibility/opening costs.",
                "step5_tail_fraction": "(1 + controls with energy <= original energy) / 101",
                "step5_tail_interpretation": "Exploratory empirical energy tail fraction among these shuffled controls. It does not represent the probability that binding is false.",
                "step4_empirical_p": "(1 + shuffles with exact-match count >= real count) / 101",
                "step4_bh_q": "Benjamini–Hochberg adjusted Step 4 p-value across all 197,195 searched transcripts.",
                "coordinates": "Predicted human and viral spans use one-based, inclusive coordinates. Exact-match region start0/end0 use zero-based, end-exclusive coordinates.",
                "shuffle_method": "The same 100 saved single-base permutations of the viral RNA were used for each full human RNA. They preserve nucleotide counts but do not necessarily preserve dinucleotide frequencies or natural RNA structure.",
                "histograms": "Common 2-kcal/mol bins and comparable axes show 100 actual control energies per candidate. Bins include their left edge and exclude their right edge, except the final bin includes its upper edge.",
            },
            "step4": {key: step4_source[key] for key in ("transcripts_tested", "shuffles_per_transcript", "bh_rejections_at_0_05", "minimum_bh_q_value", "minimum_empirical_p", "shortlist_criteria")},
            "relationships": {"all_full_sequences_distinct": len({human["sequence_sha256"] for human in plan["humans"]}) == 5,
                              "shared_block_length_nt": 626, "related_candidates": list(shared_626["transcript_intervals0"]),
                              "shared_block_intervals0": shared_626["transcript_intervals0"],
                              "selection": "The shortlist was selected from the data using exploratory criteria and was not an independently prespecified confirmation set."},
            "candidates": candidates,
        }
        data["step4"].update({"window_size": 10, "window_step": 1, "viral_windows_per_shuffle": 564,
                              "zero_match_combinations_included": True})
        verification = {
            "status": "verified", "discrepancies": [], "unique_transcript_shuffle_pairs": 500,
            "candidates_checked": 5, "controls_per_candidate": 100, "completed_with_interaction": 500,
            "failed": 0, "missing": 0, "no_interaction": 0,
            "all_500_native_csvs_receipts_requests_commands_versions_and_artifact_hashes_verified": True,
            "means_medians_inclusive_tail_counts_and_fractions_independently_recomputed": True,
            "all_expected_user_results_agree": True, "all_scientific_source_hashes_preserved": True,
            "protected_source_file_count": len(protected),
            "protected_source_index_sha256": hashlib.sha256(index_text.encode("utf-8")).hexdigest(),
            "notes": ["The earlier original-prediction summary retains its historical energy-controls-pending field. Current completion is established by the verified 500-job invocation."],
            "source_files": primary_sources,
        }
        write_csv(OUTPUT / "candidate-comparison.csv", comparison_rows, list(comparison_rows[0]))
        write_csv(OUTPUT / "control-energies.csv", control_rows, ["transcript_id", "shuffle_id", "energy_kcal_mol", "status"])
        write_json(OUTPUT / "website-data.json", data)
        write_json(OUTPUT / "verification.json", verification)
        print(json.dumps({"status": "verified", "candidates": 5, "controls": 500, "source_files_checked": len(protected),
                          "bin_edges": bin_edges, "histogram_max_count": histogram_max}, indent=2))
    except Exception as error:
        write_json(OUTPUT / "verification.json", {"status": "failed", "discrepancies": [f"{type(error).__name__}: {error}"]})
        raise


if __name__ == "__main__":
    main()
