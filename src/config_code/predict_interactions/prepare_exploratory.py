"""Prepare the verified exploratory shortlist, annotations, and an unlaunched runner."""

import argparse
import csv
import gzip
import hashlib
import json
import shutil
import time
from datetime import datetime, timezone
from difflib import SequenceMatcher
from itertools import combinations
from pathlib import Path

from Bio import SeqIO

if __package__:
    from .prepare_step5 import require, validate_matches
else:
    from prepare_step5 import require, validate_matches


PROJECT_PATH = Path(__file__).resolve().parents[3]


def file_sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def extract_transcripts(fasta_path, transcript_ids):
    requested = set(transcript_ids)
    require(len(requested) == len(transcript_ids), "Shortlist transcript IDs must be unique")
    found = {}
    opener = gzip.open if Path(fasta_path).suffix == ".gz" else open
    with opener(fasta_path, "rt", encoding="utf-8") as source:
        for record in SeqIO.parse(source, "fasta"):
            identifier = record.id.split("|")[0]
            if identifier not in requested:
                continue
            require(identifier not in found, f"Duplicate FASTA transcript {identifier}")
            dna = str(record.seq).upper()
            rna = dna.replace("T", "U")
            require(bool(rna) and set(rna) <= set("ACGU"), f"Ambiguous or empty sequence for {identifier}")
            found[identifier] = {"header": record.description, "dna_sequence": dna, "rna_sequence": rna}
    require(set(found) == requested, f"Missing exact versioned FASTA IDs: {sorted(requested - set(found))}")
    return found


def load_annotations(gtf_path, transcript_ids):
    requested = set(transcript_ids)
    annotations = {}
    opener = gzip.open if Path(gtf_path).suffix == ".gz" else open
    with opener(gtf_path, "rt", encoding="utf-8") as source:
        for line in source:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in ("transcript", "exon"):
                continue
            attributes = {}
            for item in fields[8].split(";"):
                key, separator, value = item.strip().partition(" ")
                if separator:
                    attributes[key] = value.strip().strip('"')
            identifier = attributes.get("transcript_id")
            if identifier not in requested:
                continue
            annotation = annotations.setdefault(identifier, {
                "gene_id": attributes.get("gene_id"), "gene_name": attributes.get("gene_name"),
                "transcript_name": attributes.get("transcript_name"), "chrom": fields[0],
                "strand": fields[6], "exons": [],
            })
            require(annotation["chrom"] == fields[0] and annotation["strand"] == fields[6],
                    f"Inconsistent annotation location for {identifier}")
            require(annotation["gene_id"] == attributes.get("gene_id"), f"Inconsistent gene ID for {identifier}")
            if fields[2] == "exon":
                annotation["exons"].append((int(fields[3]) - 1, int(fields[4])))
    for annotation in annotations.values():
        annotation["exons"].sort(reverse=annotation["strand"] == "-")
    return annotations


def map_transcript_interval(annotation, start, end):
    require(annotation["strand"] in ("+", "-"), "Annotation strand must be + or -")
    length = sum(right - left for left, right in annotation["exons"])
    require(0 <= start < end <= length, "Transcript interval outside annotated exons")
    offset = 0
    segments = []
    for left, right in annotation["exons"]:
        exon_length = right - left
        overlap_start, overlap_end = max(start, offset), min(end, offset + exon_length)
        if overlap_start < overlap_end:
            local_start, local_end = overlap_start - offset, overlap_end - offset
            if annotation["strand"] == "+":
                genomic_start, genomic_end = left + local_start, left + local_end
            else:
                genomic_start, genomic_end = right - local_end, right - local_start
            segments.append({"chrom": annotation["chrom"], "genomic_start0": genomic_start,
                             "genomic_end0": genomic_end, "strand": annotation["strand"]})
        offset += exon_length
    return segments


def describe_sites(rows):
    intervals = sorted((int(row["human_start"]), int(row["human_end"])) for row in rows)
    clusters = []
    for start, end in intervals:
        require(start >= 0 and end - start == 10, "Expected zero-based, end-exclusive 10-nt sites")
        if clusters and start < clusters[-1]["end0"]:
            cluster = clusters[-1]
            cluster["end0"] = max(cluster["end0"], end)
        else:
            cluster = {"start0": start, "end0": end, "match_rows": 0, "starts": set()}
            clusters.append(cluster)
        cluster["match_rows"] += 1
        cluster["starts"].add(start)
    for cluster in clusters:
        cluster["unique_human_starts"] = len(cluster.pop("starts"))
    disjoint = 0
    previous_end = -1
    for start, end in sorted(set(intervals)):
        if start >= previous_end:
            disjoint += 1
            previous_end = end
    return {
        "match_rows": len(rows), "unique_human_starts": len({start for start, _ in intervals}),
        "distinct_viral_window_sequences": len({row["viral_window"] for row in rows}),
        "human_overlap_clusters": clusters, "maximum_nonoverlapping_human_sites": disjoint,
        "human_nt_covered": sum(cluster["end0"] - cluster["start0"] for cluster in clusters),
    }


def group_identical_sequences(records, ordered_ids):
    by_sequence = {}
    representatives = {}
    distinct = []
    for identifier in ordered_ids:
        sequence = records[identifier]["rna_sequence"]
        if sequence not in by_sequence:
            by_sequence[sequence] = identifier
            distinct.append(identifier)
        representatives[identifier] = by_sequence[sequence]
    return representatives, distinct


def shared_genomic_bases(first, second):
    return sum(
        max(0, min(a["genomic_end0"], b["genomic_end0"]) - max(a["genomic_start0"], b["genomic_start0"]))
        for a in first for b in second if a["chrom"] == b["chrom"] and a["strand"] == b["strand"]
    )


def relationship_reports(ids, records, annotations, match_rows, site_stats):
    pairs = []
    possible_blocks = set()
    genomic = {}
    for identifier in ids:
        annotation = annotations.get(identifier)
        genomic[identifier] = []
        for region in site_stats[identifier]["human_overlap_clusters"]:
            possible_blocks.add(records[identifier]["rna_sequence"][region["start0"]:region["end0"]])
            if annotation and annotation["exons"]:
                genomic[identifier].extend(map_transcript_interval(annotation, region["start0"], region["end0"]))
    for first, second in combinations(ids, 2):
        a, b = records[first]["rna_sequence"], records[second]["rna_sequence"]
        shared = SequenceMatcher(None, a, b, autojunk=False).find_longest_match()
        if shared.size >= 10:
            possible_blocks.add(a[shared.a:shared.a + shared.size])
        first_gene = annotations.get(first, {}).get("gene_id")
        second_gene = annotations.get(second, {}).get("gene_id")
        first_words = {row["human_matched_sequence"] for row in match_rows[first]}
        second_words = {row["human_matched_sequence"] for row in match_rows[second]}
        pairs.append({
            "transcript_a": first, "transcript_b": second, "identical_full_rna": a == b,
            "same_annotated_gene": first_gene == second_gene if first_gene and second_gene else None,
            "longest_shared_contiguous_sequence_nt": shared.size,
            "shared_block_a_start0": shared.a, "shared_block_a_end0": shared.a + shared.size,
            "shared_block_b_start0": shared.b, "shared_block_b_end0": shared.b + shared.size,
            "shared_human_matched_10mer_sequences": len(first_words & second_words),
            "shared_genomic_match_covered_nt": shared_genomic_bases(genomic[first], genomic[second])
                if genomic[first] and genomic[second] else None,
        })
    blocks = []
    for sequence in possible_blocks:
        occurrences = {}
        contained = {}
        for identifier in ids:
            human = records[identifier]["rna_sequence"]
            starts = []
            start = human.find(sequence)
            while start >= 0:
                starts.append(start)
                start = human.find(sequence, start + 1)
            if starts:
                occurrences[identifier] = [[start, start + len(sequence)] for start in starts]
                contained[identifier] = sum(
                    any(start <= int(row["human_start"]) and int(row["human_end"]) <= start + len(sequence)
                        for start in starts) for row in match_rows[identifier]
                )
        if sum(count > 0 for count in contained.values()) >= 2:
            blocks.append({"length_nt": len(sequence), "rna_sequence": sequence,
                           "transcript_intervals0": occurrences, "match_rows_fully_contained": contained})
    blocks.sort(key=lambda block: (-len(block["transcript_intervals0"]), -block["length_nt"], block["rna_sequence"]))
    return pairs, blocks


def write_csv(path, rows, columns=None):
    require(bool(rows) or columns is not None, f"No CSV schema for {path}")
    with path.open("x", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(output, fieldnames=columns or list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_manifest(output_dir, report):
    temporary = output_dir / "preparation_manifest.json.tmp"
    temporary.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    temporary.replace(output_dir / "preparation_manifest.json")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, default=PROJECT_PATH)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args(argv)
    project = args.project.resolve()
    output_dir = (args.output_dir or project / "results/step5/exploratory_shortlist").resolve()
    output_dir.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    report = {"status": "preparing", "started_utc": datetime.now(timezone.utc).isoformat(),
              "evidence_label": "exploratory", "predictions_launched": False,
              "script_sha256": file_sha256(Path(__file__))}
    write_manifest(output_dir, report)
    try:
        paths = {
            "shortlist": project / "results/step4/transcript_comparison/exploratory_shortlist.csv",
            "human_fasta": project / "data/gencode.v50.lncRNA_transcripts.fa.gz",
            "annotation_gtf": project / "data/gencode.v50.long_noncoding_RNAs.gtf.gz",
            "real_matches": project / "results/exact_complementarity_matches.csv",
            "viral_fasta": project / "src/config_code/hpv_step2/targets/HPV16_E6starI_E7_spliced_region.rna.fasta",
            "comparison_status": project / "results/step4/transcript_comparison/run_status.json",
        }
        hashes = {name: file_sha256(path) for name, path in paths.items() if path.is_file()}
        comparison = json.loads(paths["comparison_status"].read_text(encoding="utf-8"))
        require(comparison["status"] == "complete", "Transcript comparison is incomplete")
        require(hashes["real_matches"] == comparison["inputs"]["real_matches"]["sha256"], "Real-match checksum changed")
        with paths["shortlist"].open(newline="", encoding="utf-8") as source:
            shortlist = list(csv.DictReader(source))
        require(len(shortlist) == 5, "Expected the five exploratory shortlisted transcripts")
        ids = [row["seq_id"] for row in shortlist]
        require(all(row["evidence_label"] == "exploratory" and float(row["bh_q_value"]) > 0.05 for row in shortlist),
                "Expected exploratory candidates that did not pass BH FDR 0.05")
        print("Extracting full versioned GENCODE sequences", flush=True)
        records = extract_transcripts(paths["human_fasta"], ids)
        annotations = load_annotations(paths["annotation_gtf"], ids) if paths["annotation_gtf"].is_file() else {}
        viral_record = SeqIO.read(paths["viral_fasta"], "fasta")
        viral = str(viral_record.seq).upper().replace("T", "U")
        require(set(viral) <= set("ACGU"), "Viral sequence contains ambiguous bases")
        selected_rows = {identifier: [] for identifier in ids}
        with paths["real_matches"].open(newline="", encoding="utf-8") as source:
            reader = csv.DictReader(source)
            match_columns = reader.fieldnames
            for row in reader:
                if row["hum_trspt_id"] in selected_rows:
                    selected_rows[row["hum_trspt_id"]].append(row)
        stats = {}
        for row in shortlist:
            identifier = row["seq_id"]
            sequence = records[identifier]["rna_sequence"]
            require(len(sequence) == int(row["length_nt"]), f"Shortlist length differs for {identifier}")
            require(len(selected_rows[identifier]) == int(row["real_count"]), f"Match count differs for {identifier}")
            validate_matches(selected_rows[identifier], sequence, viral_record.id, viral)
            stats[identifier] = describe_sites(selected_rows[identifier])
            header_parts = records[identifier]["header"].split("|")
            annotation = annotations.get(identifier)
            if annotation:
                require(sum(end - start for start, end in annotation["exons"]) == len(sequence), f"Exon length differs for {identifier}")
                require(annotation["gene_id"] == header_parts[1], f"FASTA/GTF gene ID differs for {identifier}")
            else:
                annotations[identifier] = {
                    "gene_id": header_parts[1] if len(header_parts) > 1 else None,
                    "gene_name": header_parts[5] if len(header_parts) > 5 else None,
                    "transcript_name": header_parts[4] if len(header_parts) > 4 else None,
                    "chrom": None, "strand": None, "exons": [],
                }
            print(f"Validated {identifier}: {len(sequence)} nt, {len(selected_rows[identifier])} exact-match rows", flush=True)
        representatives, distinct_ids = group_identical_sequences(records, ids)
        pairs, blocks = relationship_reports(ids, records, annotations, selected_rows, stats)
        annotation_headers = []
        if paths["annotation_gtf"].is_file():
            with gzip.open(paths["annotation_gtf"], "rt", encoding="utf-8") as source:
                for line in source:
                    if not line.startswith("#"):
                        break
                    annotation_headers.append(line.rstrip("\n"))

        (output_dir / "human").mkdir()
        candidates = []
        regions = []
        motifs = []
        for row in shortlist:
            identifier = row["seq_id"]
            record = records[identifier]
            sequence = record["rna_sequence"]
            annotation = annotations[identifier]
            dominant_region = max(stats[identifier]["human_overlap_clusters"], key=lambda region: region["match_rows"])
            fasta_name = f"human/{identifier}.full.rna.fasta"
            with (output_dir / fasta_name).open("x", encoding="utf-8", newline="\n") as target:
                target.write(f">{identifier} GENCODE_v50 full_transcript length={len(sequence)} orientation=5to3 exploratory=true\n")
                for start in range(0, len(sequence), 70):
                    target.write(sequence[start:start + 70] + "\n")
            candidates.append({
                "selection_rank": int(row["selection_rank"]), "transcript_id": identifier,
                "length_nt": len(sequence), "gene_id": annotation["gene_id"], "gene_name": annotation["gene_name"],
                "transcript_name": annotation["transcript_name"], "chrom": annotation["chrom"], "strand": annotation["strand"],
                "annotation_source": "GTF and FASTA header" if annotation["exons"] else "FASTA header only",
                "representative_id": representatives[identifier], "human_fasta": fasta_name,
                "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
                "T_to_U_count": record["dna_sequence"].count("T"),
                "real_match_rows": stats[identifier]["match_rows"],
                "unique_human_starts": stats[identifier]["unique_human_starts"],
                "distinct_viral_window_sequences": stats[identifier]["distinct_viral_window_sequences"],
                "overlap_cluster_count": len(stats[identifier]["human_overlap_clusters"]),
                "maximum_nonoverlapping_human_sites": stats[identifier]["maximum_nonoverlapping_human_sites"],
                "match_covered_human_nt": stats[identifier]["human_nt_covered"],
                "dominant_overlap_cluster_nt": dominant_region["end0"] - dominant_region["start0"],
                "dominant_overlap_cluster_rows": dominant_region["match_rows"],
                "bh_q_value": float(row["bh_q_value"]), "evidence_label": "exploratory",
            })
            for number, region in enumerate(stats[identifier]["human_overlap_clusters"], start=1):
                genomic = map_transcript_interval(annotation, region["start0"], region["end0"]) if annotation["exons"] else []
                regions.append({"transcript_id": identifier, "region_number": number, **region,
                                "rna_sequence": sequence[region["start0"]:region["end0"]],
                                "genomic_segments0": json.dumps(genomic)})
            for motif in sorted({match["viral_window"] for match in selected_rows[identifier]}):
                subset = [match for match in selected_rows[identifier] if match["viral_window"] == motif]
                motifs.append({"transcript_id": identifier, "viral_window": motif,
                               "human_matched_sequence": subset[0]["human_matched_sequence"], "match_rows": len(subset),
                               "unique_human_starts": len({match["human_start"] for match in subset}),
                               "distinct_viral_starts": len({match["viral_start"] for match in subset})})
        write_csv(output_dir / "candidate_manifest.csv", candidates)
        write_csv(output_dir / "candidate_real_matches.csv", [row for identifier in ids for row in selected_rows[identifier]], match_columns)
        write_csv(output_dir / "match_regions.csv", regions)
        write_csv(output_dir / "match_motifs.csv", motifs)
        write_csv(output_dir / "pairwise_relationships.csv", pairs)
        (output_dir / "shared_sequence_blocks.json").write_text(json.dumps(blocks, indent=2) + "\n", encoding="utf-8")
        (output_dir / "annotations.json").write_text(json.dumps(annotations, indent=2) + "\n", encoding="utf-8")
        viral_name = "HPV16_E6starI_E7_spliced_region.original.rna.fasta"
        shutil.copyfile(paths["viral_fasta"], output_dir / viral_name)
        shutil.copyfile(paths["shortlist"], output_dir / "exploratory_shortlist.source.csv")
        runner = Path(__file__).with_name("run_exploratory_intarna.sh")
        shutil.copyfile(runner, output_dir / runner.name)
        ubuntu_runner = (output_dir / runner.name).as_posix()
        if output_dir.drive:
            ubuntu_runner = f"/mnt/{output_dir.drive[0].lower()}{ubuntu_runner[2:]}"
        with (output_dir / "prediction_inputs.tsv").open("x", newline="", encoding="utf-8") as output:
            writer = csv.writer(output, delimiter="\t", lineterminator="\n")
            writer.writerow(["representative_id", "human_fasta", "viral_fasta", "sequence_sha256", "transcript_ids"])
            for representative in distinct_ids:
                members = [identifier for identifier in ids if representatives[identifier] == representative]
                candidate = next(candidate for candidate in candidates if candidate["transcript_id"] == representative)
                writer.writerow([representative, candidate["human_fasta"], viral_name, candidate["sequence_sha256"], ",".join(members)])

        lines = ["All five candidates are exploratory. None passed BH FDR 0.05.", "",
                 "Full GENCODE v50 transcripts were extracted by exact versioned ID. RNA FASTAs retain original 5-prime to 3-prime base order with T converted to U. "
                 "The original viral RNA FASTA was copied byte for byte. Gene names below reproduce the available GENCODE annotations.", "",
                 "| Transcript | Length nt | Gene ID | Gene name | Match rows | Overlap clusters | Maximum disjoint 10-nt sites |",
                 "| --- | --- | --- | --- | --- | --- | --- |"]
        for candidate in candidates:
            lines.append(f"| {candidate['transcript_id']} | {candidate['length_nt']} | {candidate['gene_id']} | {candidate['gene_name']} | "
                         f"{candidate['real_match_rows']} | {candidate['overlap_cluster_count']} | {candidate['maximum_nonoverlapping_human_sites']} |")
        lines.extend(["", f"There are {len(distinct_ids)} distinct full RNA sequences, so the runner prepares {len(distinct_ids)} predictions. "
                      "Sequence grouping uses exact full RNA equality. Gene membership and shared genomic match coverage are recorded in pairwise_relationships.csv.", "",
                      "Overlap clusters merge 10-nt windows that share at least one human base. Adjacent windows form separate clusters. "
                      "The maximum disjoint count uses sorted unique human starts and accepts a site when its start is at least the previous accepted end. "
                      "These counts describe geometry. Repeated sequences can recur at multiple positions. "
                      "match_motifs.csv records repeated human and viral 10-mers, and match_regions.csv records cluster sequences and genomic mappings.", ""])
        for candidate in candidates:
            lines.append(f"{candidate['transcript_id']} has {candidate['dominant_overlap_cluster_rows']} of its "
                         f"{candidate['real_match_rows']} match rows within one {candidate['dominant_overlap_cluster_nt']}-nt overlap cluster.")
            lines.append("")
        lines.extend([
                      "Viral GA-rich motifs can produce many overlapping hits to complementary CCU-rich human repeats. "
                      "The raw row count includes each human start. Motif multiplicities are tabulated below.", "",
                      "| Transcript | Viral 10-mer | Human complementary 10-mer | Match rows | Distinct viral starts |",
                      "| --- | --- | --- | --- | --- |"])
        for identifier in ids:
            prominent = sorted([motif for motif in motifs if motif["transcript_id"] == identifier],
                               key=lambda motif: (-motif["match_rows"], motif["viral_window"]))[:2]
            for motif in prominent:
                lines.append(f"| {identifier} | {motif['viral_window']} | {motif['human_matched_sequence']} | "
                             f"{motif['match_rows']} | {motif['distinct_viral_starts']} |")
        lines.extend(["",
                      "| Shared RNA block nt | Transcripts containing the block | Fully contained match rows by transcript |",
                      "| --- | --- | --- |"])
        for block in blocks[:6]:
            lines.append(f"| {block['length_nt']} | {', '.join(block['transcript_intervals0'])} | "
                         f"{json.dumps(block['match_rows_fully_contained'])} |")
        lines.extend(["", "Shared sequences and overlapping repeat windows make these candidates related observations. "
                      "The complete pairwise comparisons and shared block positions are available in the accompanying CSV and JSON files. "
                      "All transcript and genomic interval fields ending in 0 are zero-based and end-exclusive. "
                      "Genomic positions come from the provided GENCODE GTF. IntaRNA result coordinates will be one-based and inclusive.", "",
                      "The prepared runner uses full human RNA as target/sequence 1 and original viral RNA as query/sequence 2. "
                      "Local settings are model=X, mode=H, acc=C, tAccW=qAccW=150, tAccL=qAccL=100, seedBP=7, outMode=C, outNumber=1. "
                      "Output columns are id1,start1,end1,id2,start2,end2,E,ED1,ED2. A version check requires IntaRNA 3.4.1. "
                      "Each invocation creates a fresh predictions/run_* directory with individual results, stdout/stderr logs and run status.", "",
                      "Run from Ubuntu after reviewing these exploratory inputs. The command below uses the standard WSL drive mount for a Windows project.",
                      "", "```bash", "conda activate binf-intarna", f"bash {ubuntu_runner}", "```"])
        (output_dir / "preparation_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
        for name, expected in hashes.items():
            require(file_sha256(paths[name]) == expected, f"Source changed during preparation: {name}")
        report.update({
            "status": "prepared", "finished_utc": datetime.now(timezone.utc).isoformat(),
            "elapsed_seconds": round(time.monotonic() - started, 2),
            "inputs": {name: {"path": str(paths[name]), "sha256": digest} for name, digest in hashes.items()},
            "source_headers": {identifier: records[identifier]["header"] for identifier in ids},
            "candidate_count": len(ids), "distinct_full_sequence_count": len(distinct_ids),
            "representatives": representatives, "viral_id": viral_record.id, "viral_length_nt": len(viral),
            "coordinate_convention": "Transcript and genomic audit intervals zero-based, end-exclusive",
            "all_candidate_match_rows_reproduced_by_fresh_exact_search": True,
            "annotations": annotations, "site_statistics": stats,
            "annotation_header_lines": annotation_headers,
            "runner_sha256": file_sha256(output_dir / runner.name),
            "ubuntu_command": f"conda activate binf-intarna\nbash {ubuntu_runner}",
            "output_sha256": {path.relative_to(output_dir).as_posix(): file_sha256(path)
                              for path in sorted(output_dir.rglob("*")) if path.is_file() and path.name != "preparation_manifest.json"},
        })
        write_manifest(output_dir, report)
        print(f"Prepared {len(distinct_ids)} distinct full-sequence prediction inputs in {output_dir}")
        return report
    except BaseException as error:
        report.update({"status": "failed", "error": f"{type(error).__name__}: {error}"})
        write_manifest(output_dir, report)
        raise


if __name__ == "__main__":
    main()
