"""Checks for transcript counts, empirical tails, FDR, and exploratory selection."""

import csv
import hashlib
import io
import json
import math
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest import mock

import numpy as np
import pandas as pd


COMPLEMENTS_DIR = Path(__file__).resolve().parents[1] / "src/config_code/complements"
sys.path.insert(0, str(COMPLEMENTS_DIR))

from compare_transcripts import (  # noqa: E402
    bh_adjust,
    compare_counts,
    count_matches,
    propose_shortlist,
)
import compare_transcripts as runner  # noqa: E402


class TranscriptComparisonTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name)
        self.real_path = self.directory / "real.csv"
        self.shuffle_path = self.directory / "shuffled.csv"
        self.ids = ["A", "B", "C", "zero"]
        self.transcripts = pd.DataFrame({"seq_id": self.ids, "length": [100] * 4})
        self.real_rows = [
            ["A", 0, 10],
            ["C", 0, 10],
            ["B", 3, 13],
            ["C", 1, 11],
            ["A", 0, 10],
            ["C", 10, 20],
            ["C", 20, 30],
            ["C", 21, 31],
        ]
        self.shuffle_rows = []
        for shuffle_id in range(1, 101):
            self.shuffle_rows.append(["B", shuffle_id])
            if shuffle_id <= 4:
                self.shuffle_rows.extend([["A", shuffle_id]] * (shuffle_id - 1))
        self.write_inputs()

    def write_inputs(self):
        with self.real_path.open("w", newline="", encoding="utf-8") as target:
            writer = csv.writer(target)
            writer.writerow(["hum_trspt_id", "human_start", "human_end"])
            writer.writerows(self.real_rows)
        with self.shuffle_path.open("w", newline="", encoding="utf-8") as target:
            writer = csv.writer(target)
            writer.writerow(["hum_trspt_id", "shuffle_id"])
            writer.writerows(self.shuffle_rows)

    def count(self, transcript_ids=None, chunk_size=2):
        with redirect_stdout(io.StringIO()):
            return count_matches(
                self.real_path,
                self.shuffle_path,
                self.ids if transcript_ids is None else transcript_ids,
                num_shuffles=100,
                chunk_size=chunk_size,
            )

    def comparison(self):
        return compare_counts(self.transcripts, *self.count())

    def test_chunked_counts_preserve_rows_and_fill_all_zero_combinations(self):
        real_counts, shuffle_counts, real_starts = self.count()
        np.testing.assert_array_equal(real_counts, [2, 1, 5, 0])
        expected = np.zeros((4, 100), dtype=np.uint32)
        expected[0, :4] = [0, 1, 2, 3]
        expected[1, :] = 1
        np.testing.assert_array_equal(shuffle_counts, expected)
        self.assertEqual(real_starts[0], [0, 0])
        self.assertEqual(real_starts[2], [0, 1, 10, 20, 21])
        self.assertEqual(real_starts.get(3, []), [])

        # Counts must be independent of how repeated IDs fall across chunks.
        larger_real, larger_shuffles, larger_starts = self.count(chunk_size=1000)
        np.testing.assert_array_equal(larger_real, real_counts)
        np.testing.assert_array_equal(larger_shuffles, shuffle_counts)
        self.assertEqual(larger_starts, real_starts)

    def test_empirical_ties_zero_counts_and_sample_standard_deviation(self):
        table = self.comparison().set_index("seq_id")
        self.assertEqual(table.index.tolist(), self.ids)
        self.assertEqual(table.loc["A", "real_count"], 2)
        self.assertAlmostEqual(table.loc["A", "shuffled_mean"], 0.06)
        self.assertAlmostEqual(table.loc["A", "shuffled_sd"], math.sqrt((14 - 0.36) / 99))
        self.assertEqual(table.loc["A", "shuffles_ge_real"], 2)
        self.assertAlmostEqual(table.loc["A", "empirical_p"], 3 / 101)

        # The count equal to the real count belongs to the upper tail.
        self.assertEqual(table.loc["B", "shuffles_ge_real"], 100)
        self.assertEqual(table.loc["B", "empirical_p"], 1)
        self.assertEqual(table.loc["B", "shuffled_sd"], 0)
        self.assertEqual(table.loc["C", "shuffles_ge_real"], 0)
        self.assertAlmostEqual(table.loc["C", "empirical_p"], 1 / 101)
        self.assertEqual(table.loc["zero", "shuffles_ge_real"], 100)
        self.assertEqual(table.loc["zero", "empirical_p"], 1)
        self.assertEqual(table.loc["zero", "shuffled_mean"], 0)
        self.assertEqual(table.loc["zero", "shuffled_sd"], 0)

    def test_comparison_keeps_all_100_counts_and_uses_complete_bh_family(self):
        table = self.comparison().set_index("seq_id")
        count_columns = [f"shuffle_{number:03d}_count" for number in range(1, 101)]
        self.assertTrue(set(count_columns).issubset(table.columns))
        np.testing.assert_array_equal(table[count_columns].to_numpy(), self.count()[1])
        # Including the zero-real transcript makes the testing family size four.
        self.assertAlmostEqual(table.loc["C", "bh_q_value"], 4 / 101)
        self.assertAlmostEqual(table.loc["A", "bh_q_value"], 6 / 101)
        self.assertTrue(table.loc["C", "bh_reject_0_05"])
        self.assertFalse(table.loc["A", "bh_reject_0_05"])
        self.assertEqual(table.loc["zero", "bh_q_value"], 1)
        self.assertAlmostEqual(table.loc["A", "excess_matches"], 1.94)
        self.assertAlmostEqual(table.loc["A", "excess_matches_per_kb"], 19.4)
        self.assertEqual(table.loc["A", "length_nt"], 100)

    def test_nonoverlapping_sites_use_distinct_greedy_10_nt_intervals(self):
        table = self.comparison().set_index("seq_id")
        self.assertEqual(table.loc["A", "real_nonoverlapping_human_sites"], 1)
        self.assertEqual(table.loc["B", "real_nonoverlapping_human_sites"], 1)
        self.assertEqual(table.loc["C", "real_nonoverlapping_human_sites"], 3)
        self.assertEqual(table.loc["zero", "real_nonoverlapping_human_sites"], 0)

    def test_bh_ties_preserve_original_order_and_reverse_monotonicity(self):
        probabilities = np.array([1.0, 0.04, 0.01, 0.01])
        actual = bh_adjust(probabilities)
        np.testing.assert_allclose(actual, [1.0, 0.04 * 4 / 3, 0.02, 0.02])
        np.testing.assert_array_equal(probabilities, [1.0, 0.04, 0.01, 0.01])

    def test_unknown_or_missing_transcript_ids_are_rejected(self):
        for bad_id in ("unknown", ""):
            with self.subTest(source="real", transcript_id=bad_id):
                self.real_rows = [[bad_id, 0, 10]]
                self.write_inputs()
                with self.assertRaises(ValueError):
                    self.count()
            self.real_rows = [["A", 0, 10]]
            with self.subTest(source="shuffle", transcript_id=bad_id):
                self.shuffle_rows = [[bad_id, 1]]
                self.write_inputs()
                with self.assertRaises(ValueError):
                    self.count()
            self.shuffle_rows = []

    def test_duplicate_transcript_universe_is_rejected(self):
        with self.assertRaises(ValueError):
            self.count(transcript_ids=["A", "A", "B", "C", "zero"])

    def test_invalid_shuffle_ids_are_rejected(self):
        for bad_shuffle_id in (0, 101, 1.5, "", "invalid"):
            with self.subTest(shuffle_id=bad_shuffle_id):
                self.shuffle_rows = [["A", bad_shuffle_id]]
                self.write_inputs()
                with self.assertRaises(ValueError):
                    self.count()

    def test_invalid_real_intervals_are_rejected(self):
        for start, end in ((0, 9), (0, 11), (-1, 9), (0.5, 10.5), ("", 10)):
            with self.subTest(start=start, end=end):
                self.real_rows = [["A", start, end]]
                self.write_inputs()
                with self.assertRaises(ValueError):
                    self.count()

    def test_shortlist_filters_and_ranks_deterministically_as_exploratory(self):
        table = pd.DataFrame(
            [
                ["A", 3, 0, 1 / 101, 10, 2],
                ["C", 4, 0, 1 / 101, 10, 2],
                ["B", 4, 0, 1 / 101, 10, 2],
                ["D", 9, 0, 2 / 101, 100, 3],
                ["too_few", 2, 0, 1 / 101, 100, 2],
                ["overlapping", 20, 0, 1 / 101, 100, 1],
                ["no_excess", 3, 3, 1 / 101, 100, 2],
                ["large_p", 8, 0, 0.1, 100, 2],
            ],
            columns=[
                "seq_id", "real_count", "shuffled_mean", "empirical_p",
                "excess_matches_per_kb", "real_nonoverlapping_human_sites",
            ],
        )
        shortlist = propose_shortlist(table, size=2)
        self.assertEqual(shortlist["seq_id"].tolist(), ["B", "C"])
        self.assertEqual(shortlist["selection_rank"].tolist(), [1, 2])
        self.assertEqual(shortlist["evidence_label"].tolist(), ["exploratory"] * 2)
        self.assertEqual(len(table), 8)
        self.assertNotIn("selection_rank", table.columns)

    def test_empty_shortlist_has_no_fabricated_candidate(self):
        table = self.comparison()
        shortlist = propose_shortlist(table, size=5)
        self.assertEqual(shortlist["seq_id"].tolist(), ["C"])
        empty = propose_shortlist(table.assign(real_nonoverlapping_human_sites=0), size=5)
        self.assertTrue(empty.empty)

    def verified_fixture(self):
        """Supply small checked inputs while exercising the normal output path."""
        human_path = self.directory / "humans.csv"
        summary_path = self.directory / "shuffle_summary.csv"
        run_path = self.directory / "step4_run.json"
        audit_path = self.directory / "step4_audit.json"
        self.transcripts.to_csv(human_path, index=False)
        totals = [1, 2, 3, 4] + [1] * 96
        matched_transcripts = [1, 2, 2, 2] + [1] * 96
        summary = pd.DataFrame(
            {
                "shuffle_id": range(1, 101),
                "total_exact_matches": totals,
                "human_sequences_with_matches": matched_transcripts,
            }
        )
        summary.to_csv(summary_path, index=False)
        run = {"real_output_replay": {"match_rows": len(self.real_rows)}}
        audit = {"total_corrected_match_rows": len(self.shuffle_rows)}
        run_path.write_text(json.dumps(run), encoding="utf-8")
        audit_path.write_text(json.dumps(audit), encoding="utf-8")
        paths = {
            "human_csv": human_path,
            "real_matches": self.real_path,
            "shuffled_matches": self.shuffle_path,
            "shuffle_summary": summary_path,
            "step4_run_status": run_path,
            "step4_audit": audit_path,
        }
        hashes = {
            name: hashlib.sha256(path.read_bytes()).hexdigest()
            for name, path in paths.items()
        }
        return paths, hashes, run, audit, self.transcripts, summary

    def run_entrypoint(self, output_dir):
        return runner.main(
            [
                "--project", str(self.directory),
                "--output-dir", str(output_dir),
                "--chunk-size", "2", "--shortlist-size", "2",
            ]
        )

    def test_entrypoint_saves_complete_table_shortlist_and_input_provenance(self):
        fixture = self.verified_fixture()
        output_dir = self.directory / "comparison"
        with mock.patch.object(runner, "verify_inputs", return_value=fixture):
            with redirect_stdout(io.StringIO()):
                returned = self.run_entrypoint(output_dir)
        status = json.loads((output_dir / "run_status.json").read_text())
        self.assertEqual(status, returned)
        self.assertEqual(status["status"], "complete")
        self.assertEqual(status["chunk_size"], 2)
        self.assertEqual(status["shortlist_size"], 2)
        self.assertFalse(list(output_dir.glob("*.partial")))
        saved = pd.read_csv(output_dir / "per_transcript_comparison.csv")
        pd.testing.assert_frame_equal(saved, self.comparison(), check_dtype=False)
        shortlist = pd.read_csv(output_dir / "exploratory_shortlist.csv")
        self.assertEqual(shortlist["seq_id"].tolist(), ["C"])
        self.assertEqual(shortlist["evidence_label"].tolist(), ["exploratory"])
        summary = json.loads((output_dir / "analysis_summary.json").read_text())
        self.assertEqual(summary["transcripts_tested"], 4)
        self.assertEqual(summary["shuffles_per_transcript"], 100)
        self.assertEqual(summary["real_match_rows"], 8)
        self.assertEqual(summary["shuffled_match_rows"], 106)
        self.assertEqual(summary["transcripts_with_zero_real_matches"], 1)
        self.assertEqual(summary["shuffled_sd_ddof"], 1)
        self.assertEqual(summary["bh_rejections_at_0_05"], 1)
        self.assertTrue(summary["all_100_count_and_matched_transcript_totals_verified"])
        text = (output_dir / "analysis_summary.md").read_text(encoding="utf-8")
        self.assertIn("1/101", text)
        self.assertIn("exploratory", text)
        self.assertIn("zero real matches", text)
        for name, path in fixture[0].items():
            expected_hash = fixture[1][name]
            self.assertEqual(status["inputs"][name]["sha256"], expected_hash)
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), expected_hash)

    def test_entrypoint_failure_keeps_partial_output_and_failed_status(self):
        fixture = self.verified_fixture()
        output_dir = self.directory / "interrupted_comparison"
        original_writer = pd.DataFrame.to_csv
        calls = 0

        def interrupt_second_write(frame, *args, **kwargs):
            nonlocal calls
            calls += 1
            if calls == 2:
                raise RuntimeError("forced output interruption")
            return original_writer(frame, *args, **kwargs)

        with mock.patch.object(runner, "verify_inputs", return_value=fixture):
            with mock.patch.object(pd.DataFrame, "to_csv", interrupt_second_write):
                with redirect_stdout(io.StringIO()):
                    with self.assertRaisesRegex(RuntimeError, "forced output interruption"):
                        self.run_entrypoint(output_dir)
        status = json.loads((output_dir / "run_status.json").read_text())
        self.assertEqual(status["status"], "failed")
        self.assertIn("forced output interruption", status["error"])
        self.assertTrue((output_dir / "per_transcript_comparison.csv.partial").exists())
        self.assertFalse((output_dir / "per_transcript_comparison.csv").exists())
        self.assertFalse((output_dir / "exploratory_shortlist.csv").exists())
        for name, path in fixture[0].items():
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), fixture[1][name])

    def test_entrypoint_rejects_existing_directory_without_changing_its_files(self):
        output_dir = self.directory / "existing_comparison"
        output_dir.mkdir()
        sentinels = {
            "per_transcript_comparison.csv": b"existing counts\r\n",
            "run_status.json": b'{"status":"complete"}\n',
            "sentinel.bin": b"\x00\xffpreserve these bytes",
        }
        for name, content in sentinels.items():
            (output_dir / name).write_bytes(content)
        with self.assertRaises(FileExistsError):
            self.run_entrypoint(output_dir)
        self.assertEqual(
            {path.name: path.read_bytes() for path in output_dir.iterdir()},
            sentinels,
        )


if __name__ == "__main__":
    unittest.main()
