"""Run reproducible, overlapping-window controls against every input lncRNA."""

import argparse
import csv
import hashlib
import json
import random
import time
from datetime import datetime, timezone
from itertools import zip_longest
from pathlib import Path

if __package__:
    from .complementary import (
        build_match_index, iter_indexed_matches, load_human_sequences,
        prepare_windows, read_fasta,
    )
else:
    from complementary import (
        build_match_index, iter_indexed_matches, load_human_sequences,
        prepare_windows, read_fasta,
    )


PROJECT_PATH = Path(__file__).resolve().parents[3]
WINDOW_SIZE = 10
WINDOW_STEP = 1
NUM_SHUFFLES = 100
RANDOM_SEED = 42
TEST_LIMIT = None

MATCH_COLUMNS = [
    "hum_trspt_id", "viral_rna_id", "human_start", "human_end",
    "viral_start", "viral_end", "viral_window", "viral_reverse_complement",
    "human_matched_sequence", "matches", "match_length", "mismatches", "gaps",
]
SUMMARY_COLUMNS = [
    "shuffle_id", "total_exact_matches", "human_sequences_with_matches",
    "human_sequences_searched", "viral_windows_searched", "window_size",
    "window_step", "random_seed",
]


def shuffled(sequence: str, rng: random.Random) -> str:
    bases = list(sequence)
    rng.shuffle(bases)
    return "".join(bases)


def file_sha256(path: Path) -> str:
    with path.open("rb") as source:
        return hashlib.file_digest(source, "sha256").hexdigest()


def write_status(output_dir: Path, status: dict) -> None:
    temporary = output_dir / "run_status.json.tmp"
    temporary.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    temporary.replace(output_dir / "run_status.json")


def match_values(human_ids, viral_id, match):
    human_index, viral_start, viral_end, chunk, comp, human_start = match
    return (
        human_ids[human_index], viral_id, human_start, human_start + WINDOW_SIZE,
        viral_start, viral_end, chunk, comp, comp, WINDOW_SIZE, WINDOW_SIZE, 0, 0,
    )


def verify_real_output(path, human_ids, viral_id, windows, index):
    """Compare the entire existing real CSV, including row order and duplicates."""
    count = 0
    with path.open(newline="", encoding="utf-8") as source:
        reader = csv.reader(source)
        if next(reader, None) != MATCH_COLUMNS:
            raise ValueError(f"Unexpected real-match CSV header in {path}")
        expected = (
            [str(value) for value in match_values(human_ids, viral_id, match)]
            for match in iter_indexed_matches(windows, index)
        )
        for count, (actual, replayed) in enumerate(
            zip_longest(reader, expected), start=1
        ):
            if actual != replayed:
                raise ValueError(f"Real-output replay differs at match row {count}")
    return {"status": "verified", "match_rows": count, "sha256": file_sha256(path)}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=PROJECT_PATH / "results/step4")
    parser.add_argument("--num-shuffles", type=int, default=NUM_SHUFFLES)
    parser.add_argument("--test-limit", type=int, default=TEST_LIMIT)
    args = parser.parse_args(argv)
    if args.num_shuffles < 1:
        parser.error("--num-shuffles must be positive")
    if args.test_limit is not None and args.test_limit < 1:
        parser.error("--test-limit must be positive")

    # Each run owns a fresh directory, so rerunning cannot truncate old results.
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    status = {
        "status": "running",
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "requested_shuffles": args.num_shuffles,
        "completed_shuffles": 0,
        "window_size": WINDOW_SIZE,
        "window_step": WINDOW_STEP,
        "random_seed": RANDOM_SEED,
        "test_limit": args.test_limit,
        "coordinates": "zero-based, end-exclusive",
        "shuffle_method": "single-base permutation of original viral RNA per replicate",
        "matches_file": "shuffled_matches.csv.partial",
        "summary_file": "shuffle_summary.csv",
        "code_sha256": {
            name: file_sha256(Path(__file__).with_name(name))
            for name in ("complementary.py", "shuffled.py")
        },
    }
    write_status(output_dir, status)

    try:
        viral_path = PROJECT_PATH / (
            "src/config_code/hpv_step2/targets/"
            "HPV16_E6starI_E7_spliced_region.rna.fasta"
        )
        human_path = PROJECT_PATH / "results/lncrna_checks.csv"
        viral_id, viral_sequence = read_fasta(viral_path)
        if len(viral_sequence) < WINDOW_SIZE:
            raise ValueError("The viral sequence is shorter than WINDOW_SIZE")
        human_data = load_human_sequences(human_path, args.test_limit)
        if human_data.empty:
            raise ValueError("No valid human sequences were found")
        human_ids = human_data["seq_id"].tolist()
        status.update({
            "viral_rna_id": viral_id,
            "viral_length": len(viral_sequence),
            "human_sequences_searched": len(human_data),
            "human_input_coverage": human_data.attrs,
            "inputs": {
                "human_csv": {"path": str(human_path), "sha256": file_sha256(human_path)},
                "viral_fasta": {"path": str(viral_path), "sha256": file_sha256(viral_path)},
            },
        })
        write_status(output_dir, status)

        rng = random.Random(RANDOM_SEED)
        controls = []
        real_windows = prepare_windows(viral_sequence, WINDOW_SIZE, WINDOW_STEP)
        queries = {window[3] for window in real_windows}
        with (output_dir / "viral_shuffles.fasta").open("x", encoding="utf-8") as fasta:
            for shuffle_id in range(1, args.num_shuffles + 1):
                sequence = shuffled(viral_sequence, rng)
                windows = prepare_windows(sequence, WINDOW_SIZE, WINDOW_STEP)
                controls.append(windows)
                queries.update(window[3] for window in windows)
                fasta.write(f">{viral_id}_shuffle_{shuffle_id:03d}\n{sequence}\n")

        print(
            f"Indexing {len(human_data):,} human transcripts for "
            f"{len(queries):,} distinct query {WINDOW_SIZE}-mers", flush=True,
        )
        index = build_match_index(human_data["seq"], queries, WINDOW_SIZE)
        real_path = PROJECT_PATH / "results/exact_complementarity_matches.csv"
        if args.test_limit is None and real_path.exists():
            status["real_output_replay"] = verify_real_output(
                real_path, human_ids, viral_id, real_windows, index
            )
            print(
                f"Verified all {status['real_output_replay']['match_rows']:,} "
                "existing real-match rows", flush=True,
            )
        write_status(output_dir, status)

        partial_path = output_dir / "shuffled_matches.csv.partial"
        with (
            partial_path.open("x", newline="", encoding="utf-8") as matches_file,
            (output_dir / "shuffle_summary.csv").open("x", newline="", encoding="utf-8") as summary_file,
        ):
            match_writer = csv.writer(matches_file)
            summary_writer = csv.writer(summary_file)
            match_writer.writerow(["shuffle_id", *MATCH_COLUMNS])
            summary_writer.writerow(SUMMARY_COLUMNS)
            summary_file.flush()
            for shuffle_id, windows in enumerate(controls, start=1):
                match_count = 0
                matched_humans = 0
                previous_human = None
                for match in iter_indexed_matches(windows, index):
                    match_writer.writerow((shuffle_id, *match_values(human_ids, viral_id, match)))
                    match_count += 1
                    if match[0] != previous_human:
                        matched_humans += 1
                        previous_human = match[0]
                # A summary row is a completion marker for this entire shuffle.
                matches_file.flush()
                summary_writer.writerow((
                    shuffle_id, match_count, matched_humans, len(human_data),
                    len(windows), WINDOW_SIZE, WINDOW_STEP, RANDOM_SEED,
                ))
                summary_file.flush()
                status["completed_shuffles"] = shuffle_id
                status["elapsed_seconds"] = round(time.monotonic() - started, 2)
                write_status(output_dir, status)
                print(
                    f"Finished shuffle {shuffle_id}/{args.num_shuffles}: "
                    f"{match_count:,} exact matches, "
                    f"{status['elapsed_seconds']:.1f}s elapsed", flush=True,
                )

        partial_path.rename(output_dir / "shuffled_matches.csv")
        status.update({
            "status": "complete",
            "matches_file": "shuffled_matches.csv",
            "finished_utc": datetime.now(timezone.utc).isoformat(),
            "elapsed_seconds": round(time.monotonic() - started, 2),
        })
        write_status(output_dir, status)
        print(f"Complete. Results saved to {output_dir}", flush=True)
        return status
    except BaseException as error:
        status.update({
            "status": "failed",
            "error": f"{type(error).__name__}: {error}",
            "elapsed_seconds": round(time.monotonic() - started, 2),
        })
        write_status(output_dir, status)
        raise


if __name__ == "__main__":
    main()
