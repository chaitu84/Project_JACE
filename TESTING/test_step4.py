"""Regression checks for the optimized exact-complementarity search.

Run with ``python -m unittest discover -s TESTING -p test_step4.py -v``.
"""

import csv
import hashlib
import io
import json
import random
import sys
import tempfile
import unittest
from collections import Counter
from contextlib import redirect_stdout
from pathlib import Path
from unittest import mock


COMPLEMENTS_DIR = (
    Path(__file__).resolve().parents[1] / "src/config_code/complements"
)
sys.path.insert(0, str(COMPLEMENTS_DIR))

from complementary import (  # noqa: E402
    build_match_index,
    find_matches,
    iter_indexed_matches,
    load_human_sequences,
    prepare_windows,
    rev_comp,
)
from shuffled import shuffled  # noqa: E402
import shuffled as runner  # noqa: E402


def naive_matches(human_sequences, windows):
    """Enumerate every position independently of the optimized index."""
    matches = []
    for human_index, sequence in enumerate(human_sequences):
        for viral_start, viral_end, chunk, complement in windows:
            for human_start in range(len(sequence) - len(complement) + 1):
                if sequence[human_start:human_start + len(complement)] == complement:
                    matches.append(
                        (
                            human_index,
                            viral_start,
                            viral_end,
                            chunk,
                            complement,
                            human_start,
                        )
                    )
    return matches


def indexed_matches(human_sequences, windows, window_size=10):
    index = build_match_index(
        human_sequences,
        [window[3] for window in windows],
        window_size=window_size,
    )
    return list(iter_indexed_matches(windows, index))


class ExactComplementarityRegressionTests(unittest.TestCase):
    def test_reference_find_matches_includes_overlaps(self):
        self.assertEqual(find_matches("AAA", "AAAAA"), [0, 1, 2])

    def test_overlap_duplicate_windows_and_half_open_coordinates(self):
        humans = ["A" * 14, "C" * 12, "A" * 10]
        windows = prepare_windows("U" * 12)
        matches = indexed_matches(humans, windows)

        self.assertEqual([window[0] for window in windows], [0, 1, 2])
        self.assertEqual(len(matches), 18)
        self.assertEqual(matches, naive_matches(humans, windows))
        self.assertEqual(
            Counter(match[1] for match in matches),
            {0: 6, 1: 6, 2: 6},
        )
        for human_index, viral_start, viral_end, chunk, complement, start in matches:
            self.assertEqual(viral_end - viral_start, 10)
            self.assertEqual(humans[human_index][start:start + 10], complement)
            self.assertEqual(rev_comp(chunk), complement)
        self.assertIn((0, 0, 10, "U" * 10, "A" * 10, 0), matches)
        self.assertIn((0, 2, 12, "U" * 10, "A" * 10, 4), matches)

    def test_ambiguous_bases_do_not_create_crossing_matches(self):
        humans = ["A" * 10 + "N" + "A" * 10, "A" * 9 + "R", ""]
        windows = prepare_windows("U" * 10 + "N" + "U" * 10)
        matches = indexed_matches(humans, windows)

        self.assertEqual([window[0] for window in windows], [0, 11])
        self.assertEqual(matches, naive_matches(humans, windows))
        self.assertEqual(len(matches), 4)
        self.assertEqual({match[5] for match in matches}, {0, 11})

    def test_repeated_input_rows_retain_separate_occurrences(self):
        humans = ["ACGUACGUAC", "ACGUACGUAC"]
        windows = prepare_windows(rev_comp(humans[0]))
        matches = indexed_matches(humans, windows)
        self.assertEqual([match[0] for match in matches], [0, 1])
        self.assertEqual(matches, naive_matches(humans, windows))

    def test_transcript_and_batch_boundaries(self):
        windows = prepare_windows("U" * 10)
        self.assertEqual(indexed_matches(["A" * 9, "A", ""], windows), [])

        # The first sequence fills a complete 2,000,000-base index batch.
        humans = ["C" * 1_999_989 + "A" * 10, "A" * 10]
        self.assertEqual(
            indexed_matches(humans, windows),
            [
                (0, 0, 10, "U" * 10, "A" * 10, 1_999_989),
                (1, 0, 10, "U" * 10, "A" * 10, 0),
            ],
        )

    def test_no_matches_and_short_sequences(self):
        windows = prepare_windows("U" * 10)
        self.assertEqual(indexed_matches(["A" * 9, "C" * 20, ""], windows), [])
        self.assertEqual(indexed_matches(["A" * 20], []), [])
        self.assertEqual(indexed_matches([], windows), [])
        self.assertEqual(prepare_windows("U" * 9), [])

    def test_normalization_missing_rows_and_limit(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "human.csv"
            with path.open("w", newline="", encoding="utf-8") as output:
                writer = csv.writer(output)
                writer.writerow(["seq_id", "seq", "other_column"])
                writer.writerows(
                    [
                        ["", "AAAA", "missing ID"],
                        ["missing_seq", "", "missing sequence"],
                        ["blank_seq", "   ", "empty after normalization"],
                        ["duplicate", "  acgt  ", "first valid"],
                        ["duplicate", "ttAA", "second valid"],
                        ["last", "GGcc", "third valid"],
                    ]
                )

            humans = load_human_sequences(path)
            self.assertEqual(humans["seq_id"].tolist(), ["duplicate", "duplicate", "last"])
            self.assertEqual(humans["seq"].tolist(), ["ACGU", "UUAA", "GGCC"])
            limited = load_human_sequences(path, test_limit=1)
            self.assertEqual(limited["seq"].tolist(), ["ACGU"])

        windows = prepare_windows("  atgc  ", window_size=3)
        self.assertEqual(windows, [(0, 3, "AUG", "CAU"), (1, 4, "UGC", "GCA")])

    def test_randomized_index_matches_naive_search_in_original_order(self):
        rng = random.Random(7281)
        for window_size in (1, 3, 10):
            for case in range(15):
                viral = "".join(rng.choices("ACGU", k=window_size + 15))
                humans = [
                    "".join(rng.choices("ACGUN", k=rng.randrange(0, 60)))
                    for _ in range(5)
                ]
                # Guarantee some hits, including two occurrences of a query.
                query = rev_comp(viral[2:2 + window_size])
                humans.append(query + "N" + query)
                for step in (1, 10):
                    with self.subTest(window_size=window_size, case=case, step=step):
                        windows = prepare_windows(viral, window_size, step)
                        self.assertEqual(
                            indexed_matches(humans, windows, window_size),
                            naive_matches(humans, windows),
                        )

    def test_seed42_shuffle_stream_and_base_multiplicity(self):
        sequence = "ACGUACGUAAUUCCGG"
        rng = random.Random(42)
        actual = [shuffled(sequence, rng) for _ in range(3)]
        self.assertEqual(
            actual,
            ["UACGGUCACGCGAUAU", "CCCCGAGUGUUGAAAU", "CUAUGACAGGCCUGAU"],
        )
        for result in actual:
            self.assertEqual(Counter(result), Counter(sequence))


class ShuffledRunnerIntegrationTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.project = Path(temporary.name)
        self.output_dir = self.project / "results/step4"
        self.viral_sequence = "ACGUACGUACGUAC"
        self.viral_id = "fixture_viral"
        self.control_sequences = []
        rng = random.Random(42)
        for _ in range(2):
            bases = list(self.viral_sequence)
            rng.shuffle(bases)
            self.control_sequences.append("".join(bases))

        target_dir = self.project / "src/config_code/hpv_step2/targets"
        target_dir.mkdir(parents=True)
        self.viral_path = target_dir / "HPV16_E6starI_E7_spliced_region.rna.fasta"
        self.viral_path.write_text(
            f">{self.viral_id}\n{self.viral_sequence}\n", encoding="utf-8"
        )
        results_dir = self.project / "results"
        results_dir.mkdir()
        self.human_path = results_dir / "lncrna_checks.csv"
        self.human_ids = ["first", "repeat", "repeat"]
        self.human_sequences = [
            rev_comp(self.control_sequences[0]) * 2,
            rev_comp(self.control_sequences[1]),
            "A" * 10,
        ]
        with self.human_path.open("w", newline="", encoding="utf-8") as output:
            writer = csv.writer(output)
            writer.writerow(["seq_id", "seq"])
            writer.writerows(
                [
                    ["", "ACGUACGUACGU"],
                    ["missing_sequence", ""],
                    ["blank_sequence", "   "],
                    [self.human_ids[0], "  " + self.human_sequences[0].replace("U", "T").lower() + "  "],
                    [self.human_ids[1], self.human_sequences[1]],
                    [self.human_ids[2], self.human_sequences[2]],
                ]
            )

    def run_controls(self):
        return runner.main(
            ["--output-dir", str(self.output_dir), "--num-shuffles", "2"]
        )

    def read_csv(self, filename):
        with (self.output_dir / filename).open(newline="", encoding="utf-8") as source:
            return list(csv.DictReader(source))

    def test_full_entrypoint_outputs_match_naive_search(self):
        with mock.patch.object(runner, "PROJECT_PATH", self.project):
            with redirect_stdout(io.StringIO()):
                returned_status = self.run_controls()

        status = json.loads((self.output_dir / "run_status.json").read_text())
        self.assertEqual(status, returned_status)
        self.assertEqual(status["status"], "complete")
        self.assertEqual(status["requested_shuffles"], 2)
        self.assertEqual(status["completed_shuffles"], 2)
        self.assertEqual(status["window_size"], 10)
        self.assertEqual(status["window_step"], 1)
        self.assertEqual(status["random_seed"], 42)
        self.assertEqual(status["coordinates"], "zero-based, end-exclusive")
        self.assertEqual(status["human_sequences_searched"], 3)
        self.assertEqual(
            status["human_input_coverage"],
            {"input_rows": 6, "valid_rows_before_limit": 3, "excluded_rows": 3, "test_limit": None},
        )
        self.assertEqual(status["matches_file"], "shuffled_matches.csv")
        self.assertFalse((self.output_dir / "shuffled_matches.csv.partial").exists())
        self.assertEqual(
            set(status["code_sha256"]), {"complementary.py", "shuffled.py"}
        )
        for input_name, path in (("human_csv", self.human_path), ("viral_fasta", self.viral_path)):
            self.assertEqual(
                status["inputs"][input_name]["sha256"],
                hashlib.sha256(path.read_bytes()).hexdigest(),
            )

        expected_fasta = "".join(
            f">{self.viral_id}_shuffle_{number:03d}\n{sequence}\n"
            for number, sequence in enumerate(self.control_sequences, start=1)
        )
        self.assertEqual(
            (self.output_dir / "viral_shuffles.fasta").read_text(), expected_fasta
        )

        expected_rows = []
        expected_summary = []
        for shuffle_id, sequence in enumerate(self.control_sequences, start=1):
            windows = prepare_windows(sequence)
            matches = naive_matches(self.human_sequences, windows)
            self.assertTrue(matches)
            for human_index, viral_start, viral_end, chunk, complement, start in matches:
                expected_rows.append(
                    {
                        "shuffle_id": str(shuffle_id),
                        "hum_trspt_id": self.human_ids[human_index],
                        "viral_rna_id": self.viral_id,
                        "human_start": str(start),
                        "human_end": str(start + 10),
                        "viral_start": str(viral_start),
                        "viral_end": str(viral_end),
                        "viral_window": chunk,
                        "viral_reverse_complement": complement,
                        "human_matched_sequence": self.human_sequences[human_index][start:start + 10],
                        "matches": "10",
                        "match_length": "10",
                        "mismatches": "0",
                        "gaps": "0",
                    }
                )
            expected_summary.append(
                {
                    "shuffle_id": str(shuffle_id),
                    "total_exact_matches": str(len(matches)),
                    "human_sequences_with_matches": str(len({match[0] for match in matches})),
                    "human_sequences_searched": "3",
                    "viral_windows_searched": "5",
                    "window_size": "10",
                    "window_step": "1",
                    "random_seed": "42",
                }
            )
        self.assertEqual(self.read_csv("shuffled_matches.csv"), expected_rows)
        self.assertEqual(self.read_csv("shuffle_summary.csv"), expected_summary)

    def test_interrupted_shuffle_keeps_partial_rows_and_completed_summary(self):
        original_iterator = runner.iter_indexed_matches
        calls = 0

        def interrupt_second_shuffle(windows, index):
            nonlocal calls
            calls += 1
            matches = original_iterator(windows, index)
            if calls == 1:
                yield from matches
                return
            yield next(matches)
            raise RuntimeError("forced interruption after a match row")

        with mock.patch.object(runner, "PROJECT_PATH", self.project):
            with mock.patch.object(runner, "iter_indexed_matches", interrupt_second_shuffle):
                with redirect_stdout(io.StringIO()):
                    with self.assertRaisesRegex(RuntimeError, "forced interruption"):
                        self.run_controls()

        status = json.loads((self.output_dir / "run_status.json").read_text())
        self.assertEqual(status["status"], "failed")
        self.assertEqual(status["completed_shuffles"], 1)
        self.assertIn("forced interruption", status["error"])
        self.assertEqual(status["matches_file"], "shuffled_matches.csv.partial")
        self.assertFalse((self.output_dir / "shuffled_matches.csv").exists())
        partial_rows = self.read_csv("shuffled_matches.csv.partial")
        self.assertEqual(Counter(row["shuffle_id"] for row in partial_rows)["2"], 1)
        summary_rows = self.read_csv("shuffle_summary.csv")
        self.assertEqual([row["shuffle_id"] for row in summary_rows], ["1"])
        self.assertEqual(
            int(summary_rows[0]["total_exact_matches"]),
            sum(row["shuffle_id"] == "1" for row in partial_rows),
        )

    def test_existing_output_directory_preserves_every_file(self):
        self.output_dir.mkdir()
        sentinels = {
            "shuffled_matches.csv": b"old matches\r\n",
            "shuffle_summary.csv": b"old summary\r\n",
            "run_status.json": b'{"status":"existing"}\n',
            "sentinel.bin": b"\x00\x01preserve these bytes\xff",
        }
        for filename, content in sentinels.items():
            (self.output_dir / filename).write_bytes(content)

        with mock.patch.object(runner, "PROJECT_PATH", self.project):
            with self.assertRaises(FileExistsError):
                self.run_controls()

        self.assertEqual(
            {path.name: path.read_bytes() for path in self.output_dir.iterdir()},
            sentinels,
        )


if __name__ == "__main__":
    unittest.main()
