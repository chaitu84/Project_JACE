"""Checks for full-sequence preparation, annotation mapping, and site redundancy."""

import csv
import gzip
import hashlib
import io
import json
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest import mock


PREPARATION_DIR = Path(__file__).resolve().parents[1] / "src/config_code/predict_interactions"
sys.path.insert(0, str(PREPARATION_DIR))

import prepare_exploratory as preparation  # noqa: E402


class ExploratoryPreparationTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name)

    def write_fasta(self, records, compressed=True):
        path = self.directory / ("transcripts.fa.gz" if compressed else "transcripts.fa")
        opener = gzip.open if compressed else open
        with opener(path, "wt", encoding="utf-8") as target:
            for header, sequence in records:
                target.write(f">{header}\n{sequence}\n")
        return path

    def test_extracts_exact_versioned_full_transcripts_in_original_orientation(self):
        source = self.write_fasta([
            ("TX1.1|GENE1.2|-|-|Example-201|Example|12|", "acgtttacgtaa"),
            ("TX1.2|GENE1.2|-|-|Example-202|Example|12|", "GGGGCCCCAAAA"),
            ("TX2.3|GENE2.1|-|-|Other-201|Other|10|", "TTTTACGTAC"),
        ])
        records = preparation.extract_transcripts(source, ["TX1.1", "TX2.3"])
        self.assertEqual(set(records), {"TX1.1", "TX2.3"})
        self.assertEqual(records["TX1.1"]["dna_sequence"], "ACGTTTACGTAA")
        self.assertEqual(records["TX1.1"]["rna_sequence"], "ACGUUUACGUAA")
        self.assertEqual(records["TX2.3"]["rna_sequence"], "UUUUACGUAC")
        self.assertEqual(len(records["TX1.1"]["rna_sequence"]), 12)
        self.assertEqual(
            records["TX1.1"]["header"],
            "TX1.1|GENE1.2|-|-|Example-201|Example|12|",
        )

    def test_missing_version_has_no_fallback_to_another_transcript(self):
        source = self.write_fasta([("TX1.1|GENE1.1|", "ACGTACGTACGT")])
        for requested in (["TX1.2"], ["TX1"], ["TX1.1", "absent.1"]):
            with self.subTest(requested=requested):
                with self.assertRaises(ValueError):
                    preparation.extract_transcripts(source, requested)

    def test_duplicate_source_record_for_requested_id_is_rejected(self):
        source = self.write_fasta([
            ("TX1.1|GENE1.1|", "ACGTACGTACGT"),
            ("TX1.1|GENE1.1|", "ACGTACGTACGT"),
        ])
        with self.assertRaises(ValueError):
            preparation.extract_transcripts(source, ["TX1.1"])

    def test_plain_fasta_and_rna_input_are_preserved(self):
        source = self.write_fasta([("TX1.1|GENE1.1|", "ACGUACGUACGU")], compressed=False)
        records = preparation.extract_transcripts(source, ["TX1.1"])
        self.assertEqual(records["TX1.1"]["rna_sequence"], "ACGUACGUACGU")

    def write_annotations(self):
        path = self.directory / "annotation.gtf.gz"
        rows = []
        for transcript_id, gene_id, gene_name, strand, exons in [
            ("TX1.1", "GENE1.2", "Example", "+", [(101, 110), (201, 215)]),
            ("TX2.3", "GENE1.2", "Example", "-", [(101, 110), (201, 215)]),
            ("TX1.2", "GENE_OTHER.1", "OtherVersion", "+", [(301, 325)]),
        ]:
            attributes = (
                f'gene_id "{gene_id}"; transcript_id "{transcript_id}"; '
                f'gene_name "{gene_name}"; transcript_name "{transcript_id}-name";'
            )
            rows.append("\t".join([
                "chrTest", "GENCODE", "transcript", str(exons[0][0]), str(exons[-1][1]),
                ".", strand, ".", attributes,
            ]))
            # Deliberately write exons in ascending genomic order on both strands.
            for start, end in exons:
                rows.append("\t".join([
                    "chrTest", "GENCODE", "exon", str(start), str(end),
                    ".", strand, ".", attributes,
                ]))
        with gzip.open(path, "wt", encoding="utf-8") as target:
            target.write("\n".join(rows) + "\n")
        return path

    def test_annotations_preserve_gene_versions_and_sort_exons_in_transcript_order(self):
        annotations = preparation.load_annotations(self.write_annotations(), ["TX1.1", "TX2.3"])
        self.assertEqual(set(annotations), {"TX1.1", "TX2.3"})
        for transcript_id in ("TX1.1", "TX2.3"):
            self.assertEqual(annotations[transcript_id]["gene_id"], "GENE1.2")
            self.assertEqual(annotations[transcript_id]["gene_name"], "Example")
            self.assertEqual(annotations[transcript_id]["transcript_name"], f"{transcript_id}-name")
            self.assertEqual(annotations[transcript_id]["chrom"], "chrTest")
        self.assertEqual(annotations["TX1.1"]["exons"], [(100, 110), (200, 215)])
        self.assertEqual(annotations["TX2.3"]["exons"], [(200, 215), (100, 110)])

    def annotation(self, strand):
        return {
            "chrom": "chrTest", "strand": strand,
            "exons": [(100, 110), (200, 215)] if strand == "+" else [(200, 215), (100, 110)],
        }

    def blocks(self, starts_ends, strand):
        return [
            {"chrom": "chrTest", "genomic_start0": start, "genomic_end0": end, "strand": strand}
            for start, end in starts_ends
        ]

    def test_plus_strand_interval_crossing_exon_junction(self):
        blocks = preparation.map_transcript_interval(self.annotation("+"), 8, 13)
        self.assertEqual(blocks, self.blocks([(108, 110), (200, 203)], "+"))
        self.assertEqual(sum(block["genomic_end0"] - block["genomic_start0"] for block in blocks), 5)

    def test_minus_strand_interval_crossing_exon_junction(self):
        blocks = preparation.map_transcript_interval(self.annotation("-"), 13, 18)
        self.assertEqual(blocks, self.blocks([(200, 202), (107, 110)], "-"))
        self.assertEqual(sum(block["genomic_end0"] - block["genomic_start0"] for block in blocks), 5)

    def test_exon_end_is_exclusive_and_full_minus_transcript_keeps_traversal_order(self):
        self.assertEqual(
            preparation.map_transcript_interval(self.annotation("+"), 0, 10),
            self.blocks([(100, 110)], "+"),
        )
        self.assertEqual(
            preparation.map_transcript_interval(self.annotation("-"), 0, 25),
            self.blocks([(200, 215), (100, 110)], "-"),
        )

    def test_invalid_mapping_bounds_and_strand_are_rejected(self):
        for start, end in [(-1, 10), (0, 26), (10, 10), (11, 10)]:
            with self.subTest(start=start, end=end):
                with self.assertRaises(ValueError):
                    preparation.map_transcript_interval(self.annotation("+"), start, end)
        with self.assertRaises(ValueError):
            preparation.map_transcript_interval(self.annotation("?"), 0, 10)

    def site_rows(self):
        return [
            {"human_start": str(start), "human_end": str(start + 10),
             "viral_start": str(number), "viral_end": str(number + 10),
             "viral_window": "A" * 10 if number < 4 else "C" * 10}
            for number, start in enumerate([0, 0, 1, 9, 19])
        ]

    def test_site_summary_distinguishes_repeated_rows_overlap_clusters_and_disjoint_sites(self):
        sites = preparation.describe_sites(self.site_rows())
        self.assertEqual(sites["match_rows"], 5)
        self.assertEqual(sites["unique_human_starts"], 4)
        self.assertEqual(sites["distinct_viral_window_sequences"], 2)
        self.assertEqual(sites["maximum_nonoverlapping_human_sites"], 2)
        self.assertEqual(sites["human_nt_covered"], 29)
        clusters = sites["human_overlap_clusters"]
        self.assertEqual(len(clusters), 2)
        for cluster, expected in zip(clusters, [(0, 19, 4, 3), (19, 29, 1, 1)]):
            self.assertEqual(
                (cluster["start0"], cluster["end0"], cluster["match_rows"], cluster["unique_human_starts"]),
                expected,
            )

    def test_empty_site_summary_has_no_fabricated_sites(self):
        sites = preparation.describe_sites([])
        for field in (
            "match_rows", "unique_human_starts", "distinct_viral_window_sequences",
            "maximum_nonoverlapping_human_sites", "human_nt_covered",
        ):
            self.assertEqual(sites[field], 0)
        self.assertEqual(sites["human_overlap_clusters"], [])

    def test_identical_full_rna_sequences_share_representative_in_requested_order(self):
        records = {
            "TX1.1": {"rna_sequence": "ACGUACGUACGU"},
            "TX2.3": {"rna_sequence": "ACGUACGUACGU"},
            "TX3.1": {"rna_sequence": "ACGUACGUACGA"},
        }
        mapping, representatives = preparation.group_identical_sequences(records, ["TX2.3", "TX3.1", "TX1.1"])
        self.assertEqual(mapping, {"TX2.3": "TX2.3", "TX3.1": "TX3.1", "TX1.1": "TX2.3"})
        self.assertEqual(representatives, ["TX2.3", "TX3.1"])

    def test_entrypoint_rejects_existing_output_directory_without_changing_files(self):
        output = self.directory / "prepared"
        output.mkdir()
        sentinel = output / "old_result.csv"
        sentinel.write_bytes(b"preserve existing outputs\r\n")
        with redirect_stdout(io.StringIO()):
            with self.assertRaises(FileExistsError):
                preparation.main(["--project", str(self.directory), "--output-dir", str(output)])
        self.assertEqual({path.name: path.read_bytes() for path in output.iterdir()},
                         {"old_result.csv": b"preserve existing outputs\r\n"})

    def test_complete_preparation_extracts_five_fastas_and_deduplicates_prediction_inputs(self):
        project = self.directory / "project"
        data = project / "data"
        data.mkdir(parents=True)
        comparison = project / "results/step4/transcript_comparison"
        comparison.mkdir(parents=True)
        targets = project / "src/config_code/hpv_step2/targets"
        targets.mkdir(parents=True)
        ids = [f"TX{number}.1" for number in range(1, 6)]
        sequences = ["A" * 22, "A" * 22, "A" * 25, "CC" + "A" * 22, "GG" + "A" * 22]
        genes = ["GENE1.2", "GENE1.2", "GENE1.2", "GENE2.1", "GENE3.1"]
        with gzip.open(data / "gencode.v50.lncRNA_transcripts.fa.gz", "wt", encoding="utf-8") as target:
            for identifier, sequence, gene in zip(ids, sequences, genes):
                target.write(f">{identifier}|{gene}|-|-|{identifier}-name|{gene}-name|{len(sequence)}|\n{sequence}\n")
        with gzip.open(data / "gencode.v50.long_noncoding_RNAs.gtf.gz", "wt", encoding="utf-8") as target:
            for number, (identifier, sequence, gene) in enumerate(zip(ids, sequences, genes)):
                start = 101 if number < 3 else 201 + 100 * number
                attrs = f'gene_id "{gene}"; transcript_id "{identifier}"; gene_name "{gene}-name"; transcript_name "{identifier}-name";'
                for feature in ("transcript", "exon"):
                    target.write(f"chrTest\tGENCODE\t{feature}\t{start}\t{start + len(sequence) - 1}\t.\t+\t.\t{attrs}\n")
        viral_id = "HPV16_E6starI_E7_spliced_region"
        viral_path = targets / f"{viral_id}.rna.fasta"
        viral_path.write_text(f">{viral_id} fixture\n{'U' * 12}\n", encoding="utf-8")

        real_path = project / "results/exact_complementarity_matches.csv"
        fields = [
            "hum_trspt_id", "viral_rna_id", "human_start", "human_end", "viral_start", "viral_end",
            "viral_window", "viral_reverse_complement", "human_matched_sequence", "matches", "match_length", "mismatches", "gaps",
        ]
        rows = []
        for identifier, sequence in zip(ids, sequences):
            for viral_start in range(3):
                for start in range(len(sequence) - 9):
                    if sequence[start:start + 10] == "A" * 10:
                        rows.append([identifier, viral_id, start, start + 10, viral_start, viral_start + 10,
                                     "U" * 10, "A" * 10, "A" * 10, 10, 10, 0, 0])
        with real_path.open("w", newline="", encoding="utf-8") as target:
            writer = csv.writer(target)
            writer.writerow(fields)
            writer.writerows(rows)
        (comparison / "run_status.json").write_text(json.dumps({
            "status": "complete",
            "inputs": {"real_matches": {"sha256": hashlib.sha256(real_path.read_bytes()).hexdigest()}},
        }), encoding="utf-8")
        with (comparison / "exploratory_shortlist.csv").open("w", newline="", encoding="utf-8") as target:
            writer = csv.writer(target)
            writer.writerow([
                "selection_rank", "seq_id", "length_nt", "real_count", "shuffled_mean", "shuffled_sd",
                "shuffles_ge_real", "empirical_p", "bh_q_value", "bh_reject_0_05", "excess_matches",
                "excess_matches_per_kb", "real_nonoverlapping_human_sites", "evidence_label",
            ])
            for rank, (identifier, sequence) in enumerate(zip(ids, sequences), start=1):
                count = sum(row[0] == identifier for row in rows)
                writer.writerow([rank, identifier, len(sequence), count, 0, 0, 0, 1 / 101, 1.0, False,
                                 count, count * 1000 / len(sequence), 2, "exploratory"])
        with (project / "results/lncrna_checks.csv").open("w", newline="", encoding="utf-8") as target:
            writer = csv.writer(target)
            writer.writerow(["seq_id", "seq", "length"])
            writer.writerows((identifier, sequence, len(sequence)) for identifier, sequence in zip(ids, sequences))

        output = project / "results/step5/exploratory_shortlist"
        source_bytes = {path: path.read_bytes() for path in project.rglob("*") if path.is_file()}
        with mock.patch("subprocess.run", side_effect=AssertionError("Preparation must not launch predictions")):
            with redirect_stdout(io.StringIO()):
                preparation.main(["--project", str(project), "--output-dir", str(output)])

        manifest = json.loads((output / "preparation_manifest.json").read_text(encoding="utf-8"))
        self.assertEqual(manifest["status"], "prepared")
        self.assertFalse(manifest["predictions_launched"])
        self.assertEqual(manifest["candidate_count"], 5)
        self.assertEqual(manifest["distinct_full_sequence_count"], 4)
        for identifier, expected in zip(ids, sequences):
            fasta_lines = (output / "human" / f"{identifier}.full.rna.fasta").read_text().splitlines()
            self.assertEqual(fasta_lines[0].split()[0], f">{identifier}")
            self.assertEqual("".join(fasta_lines[1:]), expected.replace("T", "U"))
        self.assertEqual(
            (output / f"{viral_id}.original.rna.fasta").read_bytes(),
            viral_path.read_bytes(),
        )
        with (output / "prediction_inputs.tsv").open(newline="", encoding="utf-8") as source:
            prediction_inputs = list(csv.DictReader(source, delimiter="\t"))
        self.assertEqual(len(prediction_inputs), 4)
        self.assertEqual([row["representative_id"] for row in prediction_inputs], [ids[0], *ids[2:]])
        self.assertTrue((output / "run_exploratory_intarna.sh").is_file())
        self.assertFalse((output / "predictions").exists())
        for path, expected in source_bytes.items():
            self.assertEqual(path.read_bytes(), expected)

        # Without a GTF, genuine header gene fields remain available and locus fields stay empty.
        (data / "gencode.v50.long_noncoding_RNAs.gtf.gz").unlink()
        fallback_output = project / "results/step5/header_only"
        with mock.patch("subprocess.run", side_effect=AssertionError("Preparation must not launch predictions")):
            with redirect_stdout(io.StringIO()):
                preparation.main(["--project", str(project), "--output-dir", str(fallback_output)])
        fallback = json.loads((fallback_output / "preparation_manifest.json").read_text(encoding="utf-8"))
        self.assertEqual(fallback["status"], "prepared")
        self.assertFalse(fallback["predictions_launched"])
        self.assertEqual(fallback["annotations"][ids[0]]["gene_id"], "GENE1.2")
        self.assertIsNone(fallback["annotations"][ids[0]]["chrom"])


if __name__ == "__main__":
    unittest.main()
