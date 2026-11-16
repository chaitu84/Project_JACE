"""Compare every searched lncRNA with the verified 100 shuffled controls."""

import argparse
import hashlib
import importlib.metadata
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_PATH = Path(__file__).resolve().parents[3]
NUM_SHUFFLES = 100
FDR_ALPHA = 0.05


def bh_adjust(p_values):
    """Benjamini-Hochberg adjusted p-values for the entire supplied family."""
    values = np.asarray(p_values, dtype=float)
    if values.ndim != 1 or not np.all(np.isfinite(values)):
        raise ValueError("p-values must be a finite one-dimensional array")
    if np.any((values < 0) | (values > 1)):
        raise ValueError("p-values must lie between zero and one")
    if not len(values):
        return values.copy()
    order = np.argsort(values, kind="stable")
    ranked = values[order] * len(values) / np.arange(1, len(values) + 1)
    adjusted = np.minimum.accumulate(ranked[::-1])[::-1].clip(max=1)
    restored = np.empty_like(adjusted)
    restored[order] = adjusted
    return restored


def count_matches(real_path, shuffled_path, transcript_ids, num_shuffles=100, chunk_size=250_000):
    """Count all rows in chunks, leaving absent transcript/shuffle pairs at zero."""
    ids = pd.Index(transcript_ids)
    if ids.has_duplicates or ids.isna().any() or any(not str(value).strip() for value in ids):
        raise ValueError("Transcript IDs must be unique and nonempty")
    if num_shuffles < 1 or chunk_size < 1:
        raise ValueError("Number of shuffles and chunk size must be positive")
    real_counts = np.zeros(len(ids), dtype=np.uint64)
    shuffle_counts = np.zeros((len(ids), num_shuffles), dtype=np.uint32)
    real_starts = {}

    with pd.read_csv(
        real_path, usecols=["hum_trspt_id", "human_start", "human_end"],
        dtype=str, chunksize=chunk_size,
    ) as chunks:
        for chunk in chunks:
            rows = ids.get_indexer(chunk["hum_trspt_id"])
            if np.any(rows < 0):
                raise ValueError("Real matches contain an unknown or missing transcript ID")
            starts = pd.to_numeric(chunk["human_start"], errors="coerce").to_numpy(dtype=float)
            ends = pd.to_numeric(chunk["human_end"], errors="coerce").to_numpy(dtype=float)
            valid = np.isfinite(starts) & np.isfinite(ends) & (starts >= 0)
            valid &= (starts == np.floor(starts)) & (ends - starts == 10)
            if not np.all(valid):
                raise ValueError("Real coordinates must describe zero-based, end-exclusive 10-nt matches")
            np.add.at(real_counts, rows, 1)
            for row, start in zip(rows, starts.astype(np.int64)):
                real_starts.setdefault(int(row), []).append(int(start))

    with pd.read_csv(
        shuffled_path, usecols=["hum_trspt_id", "shuffle_id"],
        dtype=str, chunksize=chunk_size,
    ) as chunks:
        for number, chunk in enumerate(chunks, start=1):
            rows = ids.get_indexer(chunk["hum_trspt_id"])
            if np.any(rows < 0):
                raise ValueError("Shuffled matches contain an unknown or missing transcript ID")
            shuffle_ids = pd.to_numeric(chunk["shuffle_id"], errors="coerce").to_numpy(dtype=float)
            valid = np.isfinite(shuffle_ids) & (shuffle_ids == np.floor(shuffle_ids))
            valid &= (shuffle_ids >= 1) & (shuffle_ids <= num_shuffles)
            if not np.all(valid):
                raise ValueError(f"shuffle_id must be an integer from 1 to {num_shuffles}")
            # Repeated indices need add.at so every match row increments its count.
            np.add.at(shuffle_counts, (rows, shuffle_ids.astype(np.int64) - 1), 1)
            if number == 1 or number % 10 == 0:
                print(f"Counted shuffled match chunk {number}", flush=True)
    return real_counts, shuffle_counts, real_starts


def compare_counts(transcripts, real_counts, shuffle_counts, real_starts):
    real_counts = np.asarray(real_counts)
    shuffle_counts = np.asarray(shuffle_counts)
    if shuffle_counts.shape != (len(transcripts), NUM_SHUFFLES) or real_counts.shape != (len(transcripts),):
        raise ValueError("Expected one real count and 100 shuffled counts per transcript")
    if np.any(real_counts < 0) or np.any(shuffle_counts < 0):
        raise ValueError("Match counts must be nonnegative")
    lengths = transcripts["length"].to_numpy(dtype=np.int64)
    if np.any(lengths <= 0):
        raise ValueError("Transcript lengths must be positive")
    mean = shuffle_counts.mean(axis=1)
    sd = shuffle_counts.std(axis=1, ddof=1)
    greater_equal = np.count_nonzero(shuffle_counts >= real_counts[:, None], axis=1)
    empirical_p = (greater_equal + 1) / (NUM_SHUFFLES + 1)
    q_values = bh_adjust(empirical_p)
    sites = np.zeros(len(transcripts), dtype=np.int64)
    for row, starts in real_starts.items():
        previous_end = -1
        for start in sorted(set(starts)):
            if start + 10 > lengths[row]:
                raise ValueError("Real match extends beyond its transcript length")
            if start >= previous_end:
                sites[row] += 1
                previous_end = start + 10
    excess = real_counts.astype(float) - mean
    columns = {
        "seq_id": transcripts["seq_id"].to_numpy(),
        "length_nt": lengths,
        "real_count": real_counts,
        "shuffled_mean": mean,
        "shuffled_sd": sd,
        "shuffles_ge_real": greater_equal,
        "empirical_p": empirical_p,
        "bh_q_value": q_values,
        "bh_reject_0_05": q_values <= FDR_ALPHA,
        "excess_matches": excess,
        "excess_matches_per_kb": excess * 1000 / lengths,
        "real_nonoverlapping_human_sites": sites,
    }
    for number in range(NUM_SHUFFLES):
        columns[f"shuffle_{number + 1:03d}_count"] = shuffle_counts[:, number]
    return pd.DataFrame(columns)


def propose_shortlist(comparison, size=5):
    eligible = comparison[
        (comparison["real_count"] >= 3)
        & (comparison["real_nonoverlapping_human_sites"] >= 2)
        & (comparison["real_count"] > comparison["shuffled_mean"])
        & (comparison["empirical_p"] <= 0.05)
    ]
    ranked = eligible.sort_values(
        ["empirical_p", "excess_matches_per_kb", "real_count", "seq_id"],
        ascending=[True, False, False, True], kind="stable",
    ).head(size)
    descriptive = [column for column in comparison.columns if not column.startswith("shuffle_")]
    shortlist = ranked[descriptive].copy()
    shortlist.insert(0, "selection_rank", np.arange(1, len(shortlist) + 1))
    shortlist["evidence_label"] = "exploratory"
    return shortlist.reset_index(drop=True)


def file_sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_status(output_dir, status):
    temporary = output_dir / "run_status.json.tmp"
    temporary.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    temporary.replace(output_dir / "run_status.json")


def verify_inputs(project):
    step4 = project / "results/step4"
    paths = {
        "human_csv": project / "results/lncrna_checks.csv",
        "real_matches": project / "results/exact_complementarity_matches.csv",
        "shuffled_matches": step4 / "shuffled_matches.csv",
        "shuffle_summary": step4 / "shuffle_summary.csv",
        "step4_run_status": step4 / "run_status.json",
        "step4_audit": step4 / "audit_report.json",
    }
    run = json.loads(paths["step4_run_status"].read_text(encoding="utf-8"))
    audit = json.loads(paths["step4_audit"].read_text(encoding="utf-8"))
    if run["status"] != "complete" or audit["status"] != "verified":
        raise ValueError("Step 4 must have completed and verified results")
    if run["requested_shuffles"] != 100 or run["completed_shuffles"] != 100 or run["test_limit"] is not None:
        raise ValueError("Expected 100 completed full-input shuffles")
    if (run["window_size"], run["window_step"]) != (10, 1):
        raise ValueError("Expected matched 10-mer searches at step 1")
    coverage = run["human_input_coverage"]
    if coverage["excluded_rows"] != 0:
        raise ValueError("This input universe must contain exactly the rows searched by Step 4")
    hashes = {name: file_sha256(path) for name, path in paths.items()}
    expected = {
        "human_csv": run["inputs"]["human_csv"]["sha256"],
        "real_matches": run["real_output_replay"]["sha256"],
        "shuffled_matches": audit["corrected_matches_sha256"],
    }
    for name, expected_hash in expected.items():
        if hashes[name] != expected_hash:
            raise ValueError(f"Verified input checksum differs for {name}")
    transcripts = pd.read_csv(paths["human_csv"], usecols=["seq_id", "length"], dtype={"seq_id": str, "length": np.int64})
    if len(transcripts) != run["human_sequences_searched"] or len(transcripts) != coverage["input_rows"]:
        raise ValueError("Transcript universe differs from the completed Step 4 search")
    summary = pd.read_csv(paths["shuffle_summary"])
    expected_ids = np.arange(1, 101)
    if len(summary) != 100 or not np.array_equal(summary["shuffle_id"], expected_ids):
        raise ValueError("Expected summary entries for shuffle IDs 1 through 100")
    for column, expected_value in {
        "human_sequences_searched": len(transcripts), "viral_windows_searched": 564,
        "window_size": 10, "window_step": 1, "random_seed": run["random_seed"],
    }.items():
        if not (summary[column] == expected_value).all():
            raise ValueError(f"Inconsistent shuffle summary column {column}")
    return paths, hashes, run, audit, transcripts, summary


def summary_text(summary, shortlist):
    lines = [
        f"Compared all {summary['transcripts_tested']:,} searched lncRNA transcripts with 100 shuffled controls.",
        "", "Every match row contributes one count, including overlapping human matches and repeated viral windows. "
        "The comparison CSV contains all 100 shuffle counts for every transcript, including explicit zeros. "
        "Shuffled standard deviation uses the sample formula with ddof=1.",
        "", "The empirical upper-tail p-value is (1 + shuffles_ge_real) / 101. Ties count toward the upper tail. "
        "BH correction includes every searched transcript, including transcripts with zero real matches. "
        "Zero-real transcripts have shuffles_ge_real=100 and p=1.",
        "", "With only 100 shuffles, the minimum p-value and spacing are 1/101, approximately 0.009901. "
        "A transcript exceeding all controls still has p=1/101. Finer tail probabilities require more controls. "
        "This coarse resolution limits detection after correcting for many transcripts. BH uses the ordered p-values "
        "and the total number of tested transcripts.",
        "", f"There are {summary['bh_rejections_at_0_05']:,} BH rejections at FDR 0.05. "
        f"The smallest BH-adjusted p-value is {summary['minimum_bh_q_value']:.8g}. "
        f"There are {summary['unadjusted_p_at_most_0_05']:,} transcripts with unadjusted empirical p <= 0.05.",
        "", "The exploratory shortlist requires real_count >= 3, at least two nonoverlapping human 10-nt sites, "
        "real_count above the shuffled mean, and unadjusted empirical p <= 0.05. "
        "Nonoverlapping sites are counted greedily from sorted unique human starts, using end=start+10. "
        "Ranking uses empirical p ascending, excess matches per kb descending, real count descending, "
        "then transcript ID ascending. Excess matches per kb equals 1000*(real_count-shuffled_mean)/length_nt. "
        "Length is already reflected in each transcript's shuffled baseline, and density serves as a ranking criterion. "
        "The shortlist contains up to the requested number, with no extra candidates added when too few qualify.",
        "", "The shortlist is exploratory. Raw match rows can share sequence or overlapping positions, and "
        "shortlist thresholds are descriptive choices. The shuffled viral controls preserve single-base composition. "
        "Results describe exact-match enrichment under this control model and do not establish biological binding. "
        "BH decisions use the standard Benjamini-Hochberg procedure. Transcripts can share sequence, "
        "so their tests can be dependent and this analysis does not establish the dependence assumptions required "
        "for formal BH FDR control. [BH implementation and assumptions](https://www.statsmodels.org/stable/generated/statsmodels.stats.multitest.fdrcorrection.html).",
        "", "| Rank | Transcript | Real | Shuffled mean | Shuffled SD | Shuffles >= real | Empirical p | BH q |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for row in shortlist.itertuples(index=False):
        lines.append(
            f"| {row.selection_rank} | {row.seq_id} | {row.real_count} | {row.shuffled_mean:.4g} | "
            f"{row.shuffled_sd:.4g} | {row.shuffles_ge_real} | {row.empirical_p:.6g} | {row.bh_q_value:.6g} |"
        )
    return "\n".join(lines) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, default=PROJECT_PATH)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--chunk-size", type=int, default=250_000)
    parser.add_argument("--shortlist-size", type=int, default=5)
    args = parser.parse_args(argv)
    if args.chunk_size < 1 or args.shortlist_size < 1:
        parser.error("Chunk size and shortlist size must be positive")
    project = args.project.resolve()
    output_dir = args.output_dir or project / "results/step4/transcript_comparison"
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    status = {
        "status": "running", "started_utc": datetime.now(timezone.utc).isoformat(),
        "chunk_size": args.chunk_size, "shortlist_size": args.shortlist_size,
        "script_sha256": file_sha256(Path(__file__)),
    }
    write_status(output_dir, status)
    try:
        paths, hashes, run, audit, transcripts, shuffle_summary = verify_inputs(project)
        status["inputs"] = {name: {"path": str(path), "sha256": hashes[name]} for name, path in paths.items()}
        write_status(output_dir, status)
        real, shuffled_counts, starts = count_matches(
            paths["real_matches"], paths["shuffled_matches"], transcripts["seq_id"], chunk_size=args.chunk_size,
        )
        if int(real.sum()) != run["real_output_replay"]["match_rows"]:
            raise ValueError("Real match total differs from verified output")
        if not np.array_equal(shuffled_counts.sum(axis=0), shuffle_summary["total_exact_matches"]):
            raise ValueError("Shuffled match totals differ from the completed summary")
        if not np.array_equal(np.count_nonzero(shuffled_counts, axis=0), shuffle_summary["human_sequences_with_matches"]):
            raise ValueError("Shuffled matched-transcript totals differ from the completed summary")
        if int(shuffled_counts.sum()) != audit["total_corrected_match_rows"]:
            raise ValueError("Shuffled match total differs from the independent audit")
        comparison = compare_counts(transcripts, real, shuffled_counts, starts)
        shortlist = propose_shortlist(comparison, args.shortlist_size)
        summary = {
            "transcripts_tested": len(comparison), "shuffles_per_transcript": 100,
            "real_match_rows": int(real.sum()), "shuffled_match_rows": int(shuffled_counts.sum()),
            "transcripts_with_real_matches": int(np.count_nonzero(real)),
            "transcripts_with_zero_real_matches": int(np.count_nonzero(real == 0)),
            "zero_transcript_shuffle_combinations": int(np.count_nonzero(shuffled_counts == 0)),
            "minimum_empirical_p": float(comparison["empirical_p"].min()),
            "minimum_bh_q_value": float(comparison["bh_q_value"].min()),
            "unadjusted_p_at_most_0_05": int((comparison["empirical_p"] <= 0.05).sum()),
            "bh_rejections_at_0_05": int(comparison["bh_reject_0_05"].sum()),
            "shortlisted_transcripts": len(shortlist), "shuffled_sd_ddof": 1,
            "all_100_count_and_matched_transcript_totals_verified": True,
            "p_value_formula": "(1 + count(shuffled_count >= real_count)) / 101",
            "multiple_testing_method": "Benjamini-Hochberg across all searched transcripts",
            "shortlist_criteria": {
                "minimum_real_count": 3, "minimum_nonoverlapping_human_sites": 2,
                "real_above_shuffled_mean": True, "maximum_unadjusted_empirical_p": 0.05,
                "ranking": ["empirical_p ascending", "excess_matches_per_kb descending", "real_count descending", "seq_id ascending"],
                "label": "exploratory",
            },
            "runtime": {"python": sys.version, "executable": sys.executable,
                        "packages": {name: importlib.metadata.version(name) for name in ("numpy", "pandas")}},
        }
        for filename, frame in (("per_transcript_comparison.csv", comparison), ("exploratory_shortlist.csv", shortlist)):
            frame.to_csv(output_dir / f"{filename}.partial", index=False, chunksize=5_000, mode="x")
        (output_dir / "analysis_summary.md.partial").write_text(summary_text(summary, shortlist), encoding="utf-8")
        (output_dir / "analysis_summary.json.partial").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        # Recheck inputs before marking outputs complete.
        for name, path in paths.items():
            if file_sha256(path) != hashes[name]:
                raise ValueError(f"Input changed during analysis: {name}")
        for filename in ("per_transcript_comparison.csv", "exploratory_shortlist.csv", "analysis_summary.md", "analysis_summary.json"):
            (output_dir / f"{filename}.partial").rename(output_dir / filename)
        status.update({"status": "complete", "finished_utc": datetime.now(timezone.utc).isoformat(),
                       "elapsed_seconds": round(time.monotonic() - started, 2), "summary": summary})
        write_status(output_dir, status)
        print(json.dumps(summary, indent=2))
        print(shortlist.to_string(index=False))
        return status
    except BaseException as error:
        status.update({"status": "failed", "error": f"{type(error).__name__}: {error}",
                       "elapsed_seconds": round(time.monotonic() - started, 2)})
        write_status(output_dir, status)
        raise


if __name__ == "__main__":
    main()
