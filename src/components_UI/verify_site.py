"""Verify the website's compact data and local assets against saved results."""

import bisect
import csv
import hashlib
from html.parser import HTMLParser
import json
import math
from pathlib import Path
import statistics
from urllib.parse import unquote, urlsplit


SITE = Path(__file__).resolve().parent
PROJECT = SITE.parents[1]
CONTROL_RUN = PROJECT / (
    "results/step5/exploratory_controls/results/invocations/"
    "run_20261116T051139789910Z_0558ca3f"
)
ORIGINAL_RUN = PROJECT / (
    "results/step5/exploratory_shortlist/predictions/"
    "run_20261116T044527Z_4D0EGI"
)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_csv(path, delimiter=","):
    with path.open(newline="", encoding="utf-8-sig") as source:
        return list(csv.DictReader(source, delimiter=delimiter))


def close(actual, expected, label):
    require(math.isclose(float(actual), float(expected), rel_tol=1e-9, abs_tol=1e-10),
            f"{label} differs: {actual!r} versus {expected!r}")


def digest(path):
    result = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


class PageAssets(HTMLParser):
    def __init__(self):
        super().__init__()
        self.ids = set()
        self.links = []
        self.resources = []

    def handle_starttag(self, tag, attributes):
        attrs = dict(attributes)
        if "id" in attrs:
            require(attrs["id"] not in self.ids, f"Duplicate HTML ID {attrs['id']}")
            self.ids.add(attrs["id"])
        if tag == "a" and attrs.get("href"):
            self.links.append(attrs["href"])
        if tag in ("script", "img") and attrs.get("src"):
            self.resources.append(attrs["src"])
        if tag == "link" and attrs.get("href"):
            self.resources.append(attrs["href"])


def verify():
    website = json.loads((SITE / "data/website-data.json").read_text(encoding="utf-8"))
    metadata = website["metadata"]
    candidates = website["candidates"]
    require(len(candidates) == 5 and len({row["id"] for row in candidates}) == 5,
            "Website must contain the five unique shortlisted transcript IDs")
    require(metadata["verified"] is True, "Website data is not marked verified")
    require(metadata["candidate_count"] == 5 and metadata["control_count"] == 500,
            "Website overview totals differ")
    require(metadata["controls_per_candidate"] == 100, "Wrong controls per candidate")
    for name in ("failed_controls", "missing_controls", "no_interaction_controls"):
        require(metadata[name] == 0, f"Unexpected {name}")
    require(metadata["software"]["intarna"] == "3.4.1", "Wrong IntaRNA version")
    require(metadata["software"]["viennarna"] == "2.7.0", "Wrong ViennaRNA version")
    require(metadata["target"]["length_nt"] == 573, "Wrong viral RNA length")

    original_job_rows = read_csv(CONTROL_RUN / "job_status.csv")
    require(len(original_job_rows) == 500, "Original control run does not have 500 jobs")
    actual_pairs = {}
    for job in original_job_rows:
        key = (job["transcript_id"], int(job["shuffle_id"]))
        require(key not in actual_pairs, f"Duplicate original transcript/shuffle pair {key}")
        require(job["status"] == "completed", f"Non-interaction or failed original job {key}")
        actual_pairs[key] = float(job["energy"])
        require(math.isfinite(actual_pairs[key]), f"Nonfinite saved control energy {key}")

    original_summaries = {row["transcript_id"]: row for row in read_csv(CONTROL_RUN / "summary.csv")}
    annotations = {
        row["transcript_id"]: row
        for row in read_csv(PROJECT / "results/step5/exploratory_shortlist/candidate_manifest.csv")
    }
    shortlist = {
        row["seq_id"]: row
        for row in read_csv(PROJECT / "results/step4/transcript_comparison/exploratory_shortlist.csv")
    }
    require({row["id"] for row in candidates} == set(shortlist), "Website shortlist IDs differ")
    by_id = {}
    for candidate in candidates:
        identifier = candidate["id"]
        by_id[identifier] = candidate
        energies = [actual_pairs[(identifier, number)] for number in range(1, 101)]
        require(len(candidate["energies"]) == 100, f"Missing website energies for {identifier}")
        for number, (actual, expected) in enumerate(zip(candidate["energies"], energies), start=1):
            close(actual, expected, f"{identifier} shuffle {number}")
        original_rows = read_csv(ORIGINAL_RUN / f"{identifier}.csv", delimiter=";")
        require(len(original_rows) == 1, f"Original prediction count differs for {identifier}")
        original = original_rows[0]
        close(candidate["original_energy"], original["E"], f"{identifier} original energy")
        require(candidate["predicted_human"] == [int(original["start1"]), int(original["end1"])],
                f"Human prediction coordinates differ for {identifier}")
        require(candidate["predicted_viral"] == [int(original["start2"]), int(original["end2"])],
                f"Viral prediction coordinates differ for {identifier}")
        for field in ("ED1", "ED2"):
            close(candidate[field], original[field], f"{identifier} {field}")
        mean, median = statistics.mean(energies), statistics.median(energies)
        count = sum(energy <= float(original["E"]) for energy in energies)
        close(candidate["shuffled_mean"], mean, f"{identifier} shuffled mean")
        close(candidate["shuffled_median"], median, f"{identifier} shuffled median")
        require(candidate["controls_le_original"] == count, f"Wrong tail count for {identifier}")
        close(candidate["tail_fraction"], (count + 1) / 101, f"{identifier} empirical tail fraction")
        for field in ("shuffled_mean", "shuffled_median"):
            close(candidate[field], original_summaries[identifier][field], f"Source summary {identifier} {field}")
        require(candidate["length_nt"] == int(annotations[identifier]["length_nt"]),
                f"Length differs for {identifier}")
        require(candidate["gene_id"] == annotations[identifier]["gene_id"], f"Gene ID differs for {identifier}")
        require(candidate["gene_name"] == annotations[identifier]["gene_name"], f"Gene name differs for {identifier}")
        require(candidate["overlap_clusters"] == int(annotations[identifier]["overlap_cluster_count"]),
                f"Overlap clusters differ for {identifier}")
        require(candidate["max_nonoverlapping_sites"] == int(annotations[identifier]["maximum_nonoverlapping_human_sites"]),
                f"Disjoint sites differ for {identifier}")
        require(candidate["dominant_repeat"]["length_nt"] == int(annotations[identifier]["dominant_overlap_cluster_nt"]),
                f"Repeat length differs for {identifier}")
        require(candidate["dominant_repeat"]["match_rows"] == int(annotations[identifier]["dominant_overlap_cluster_rows"]),
                f"Repeat match count differs for {identifier}")
        for field in ("real_count", "shuffled_mean", "shuffled_sd", "shuffles_ge_real", "empirical_p", "bh_q_value"):
            close(candidate["step4"][field], shortlist[identifier][field], f"Step 4 {identifier} {field}")
        require(candidate["tail_fraction"] > 0.05, f"Unexpected nominal energy-tail result for {identifier}")
        require(candidate["step4"]["bh_q_value"] > 0.05, f"Unexpected significant Step 4 result for {identifier}")
        require(candidate["original_energy"] < candidate["shuffled_mean"], f"Unexpected original-versus-mean ordering for {identifier}")
    require(sum(row["related_group"] == "shared_626nt" for row in candidates) == 4,
            "The four related candidates must retain their group")
    require(min(candidates, key=lambda row: row["tail_fraction"])["id"] == "ENST00000778197.1",
            "LINC00869 must have the smallest control-relative tail fraction")

    edges = metadata["common_bin_edges"]
    require(len(edges) >= 3 and all(a < b for a, b in zip(edges, edges[1:])),
            "Histogram bin edges must be strictly increasing")
    highest_bin = 0
    for candidate in candidates:
        require(edges[0] <= candidate["original_energy"] <= edges[-1], "Histogram omits an original energy")
        bins = [0] * (len(edges) - 1)
        for energy in candidate["energies"]:
            require(edges[0] <= energy <= edges[-1], "Histogram omits a control energy")
            index = min(bisect.bisect_right(edges, energy) - 1, len(bins) - 1)
            bins[index] += 1
        require(sum(bins) == 100, "Histogram binning drops a control")
        highest_bin = max(highest_bin, max(bins))
    require(metadata["common_histogram_max_count"] >= highest_bin,
            "Common histogram count axis clips bars")
    javascript = "\n".join(path.read_text(encoding="utf-8") for path in SITE.glob("*.js"))
    require("common_bin_edges" in javascript and "common_histogram_max_count" in javascript,
            "Frontend must consume the common histogram axes")

    downloaded_controls = read_csv(SITE / "data/control-energies.csv")
    require(len(downloaded_controls) == 500, "Control download must have 500 actual energies")
    download_pairs = {}
    for row in downloaded_controls:
        key = (row["transcript_id"], int(row["shuffle_id"]))
        require(key not in download_pairs, f"Duplicate downloaded pair {key}")
        require(row["status"] == "completed", f"Wrong downloaded outcome for {key}")
        download_pairs[key] = float(row["energy_kcal_mol"])
    require(download_pairs == actual_pairs, "Downloaded control energies differ from actual results")
    downloaded_candidates = read_csv(SITE / "data/candidate-comparison.csv")
    require(len(downloaded_candidates) == 5, "Candidate download must contain five rows")
    require({row["transcript_id"] for row in downloaded_candidates} == set(by_id),
            "Candidate download omits or repeats a transcript")
    for row in downloaded_candidates:
        candidate = by_id[row["transcript_id"]]
        require(row["gene_id"] == candidate["gene_id"] and row["gene_name"] == candidate["gene_name"],
                f"Downloaded gene annotation differs for {candidate['id']}")
        require(int(row["length_nt"]) == candidate["length_nt"] and int(row["total_controls"]) == 100,
                f"Downloaded length or control count differs for {candidate['id']}")
        require(row["related_group"] == candidate["related_group"],
                f"Downloaded relationship label differs for {candidate['id']}")
        for field, source_field in (
            ("original_energy_kcal_mol", "original_energy"),
            ("shuffled_mean_kcal_mol", "shuffled_mean"),
            ("shuffled_median_kcal_mol", "shuffled_median"),
            ("controls_energy_le_original", "controls_le_original"),
            ("exploratory_energy_tail_fraction", "tail_fraction"),
        ):
            close(row[field], candidate[source_field], f"Download {candidate['id']} {field}")
        close(row["step4_empirical_p"], candidate["step4"]["empirical_p"], "Downloaded Step 4 p")
        close(row["step4_bh_q_value"], candidate["step4"]["bh_q_value"], "Downloaded Step 4 q")

    hashes_checked = 0
    for source in metadata["sources"]:
        path = PROJECT / source["path"]
        require(path.is_file(), f"Missing provenance source {source['path']}")
        require(digest(path) == source["sha256"], f"Original scientific source changed: {source['path']}")
        hashes_checked += 1

    page = PageAssets()
    page.feed((SITE / "index.html").read_text(encoding="utf-8"))
    for url in page.resources + page.links:
        parsed = urlsplit(url)
        if parsed.scheme or parsed.netloc:
            continue
        if parsed.path:
            require((SITE / unquote(parsed.path)).is_file(), f"Broken local link or asset {url}")
        elif parsed.fragment:
            require(parsed.fragment in page.ids, f"Missing navigation destination {url}")
    report = {
        "status": "verified",
        "candidate_count": 5,
        "unique_transcript_shuffle_pairs": 500,
        "all_control_energy_values_match_original_job_status": True,
        "all_means_medians_inclusive_tail_counts_and_fractions_verified": True,
        "all_original_prediction_coordinates_and_accessibility_costs_verified": True,
        "all_candidate_genes_lengths_repeats_and_step4_statistics_verified": True,
        "histograms_share_verified_common_bin_edges_and_count_axis": True,
        "all_500_controls_counted_in_histogram_bins": True,
        "both_csv_downloads_match_verified_data": True,
        "original_source_hashes_checked": hashes_checked,
        "all_provenance_source_hashes_preserved": True,
        "local_links_assets_and_navigation_destinations_exist": True,
        "energy_tail_cutoff_0_05_passes": 0,
        "step4_bh_cutoff_0_05_passes": 0,
        "browser_checks_included": False,
    }
    (SITE / "qa_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return report


if __name__ == "__main__":
    verify()
