
import csv
import io
import tempfile
import json
import sys
import platform
from datetime import datetime, timezone
from importlib.metadata import version
from collections import Counter, defaultdict

import streamlit as st
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

# --------------------------------------------------
# PAGE CONFIGURATION
# --------------------------------------------------

st.set_page_config(
    page_title="HPV RNA Analysis",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="collapsed"
)

# --------------------------------------------------
# WEBSITE DESIGN
# --------------------------------------------------

st.markdown("""
<style>
    .stApp {
        background:
            radial-gradient(
                ellipse at 85% 25%,
                rgba(34, 75, 145, 0.25),
                transparent 50%
            ),
            linear-gradient(
                135deg,
                #030712 0%,
                #081329 50%,
                #030712 100%
            );
        color: #e8f2ff;
    }

    .block-container {
        max-width: 1250px;
        padding-top: 2rem;
        padding-bottom: 4rem;
    }

    h1, h2, h3 {
        color: #e8f2ff !important;
    }

    h1 {
        font-size: 2.7rem !important;
        font-weight: 800 !important;
    }

    p, label {
        color: #d7e6fa;
    }

    div[data-testid="stMetric"] {
        background: rgba(15, 30, 57, 0.85);
        border: 1px solid rgba(75, 175, 255, 0.25);
        padding: 18px;
        border-radius: 15px;
    }

    div[data-testid="stMetricValue"] {
        color: #5de4ff;
    }

    div[data-testid="stDataFrame"] {
        border: 1px solid rgba(75, 175, 255, 0.25);
        border-radius: 12px;
        overflow: hidden;
    }

    .stButton > button[kind="primary"] {
        background: linear-gradient(
            90deg,
            #176bd7,
            #13b8ca
        );
        color: white;
        border: none;
        border-radius: 10px;
        font-weight: 700;
        padding: 0.7rem 1.5rem;
    }

    .stButton > button[kind="primary"]:hover {
        background: linear-gradient(
            90deg,
            #2582f2,
            #20d4e5
        );
        color: white;
    }

    .hero-card {
        position: relative;
        overflow: hidden;
        min-height: 280px;
        border: 1px solid rgba(87, 175, 255, 0.3);
        border-radius: 22px;
        background:
            linear-gradient(
                110deg,
                rgba(9, 22, 48, 0.98) 0%,
                rgba(9, 22, 48, 0.92) 48%,
                rgba(15, 40, 75, 0.7) 100%
            );
        padding: 40px;
        margin-bottom: 28px;
    }

    .hero-title {
        position: relative;
        z-index: 2;
        font-size: 2.5rem;
        font-weight: 800;
        line-height: 1.15;
        color: #f0f8ff;
        max-width: 610px;
    }

    .hero-subtitle {
        position: relative;
        z-index: 2;
        color: #b5cbe5;
        font-size: 1.05rem;
        line-height: 1.7;
        max-width: 530px;
        margin-top: 16px;
    }

    .hero-tag {
        position: relative;
        z-index: 2;
        color: #5de4ff;
        letter-spacing: 2px;
        font-size: 0.8rem;
        font-weight: 700;
        margin-bottom: 18px;
    }

    .dna-scene {
        position: absolute;
        right: 5%;
        top: 8%;
        width: 220px;
        height: 250px;
        perspective: 700px;
        opacity: 0.9;
        pointer-events: none;
    }

    .dna-helix {
        position: absolute;
        left: 50%;
        top: 50%;
        transform-style: preserve-3d;
        animation: helix-spin 9s linear infinite;
    }

    .dna-pair {
        position: absolute;
        width: 150px;
        height: 8px;
        left: -75px;
        top: 0;
        transform-style: preserve-3d;
    }

    .dna-bar {
        position: absolute;
        left: 10px;
        right: 10px;
        top: 3px;
        height: 2px;
        background: linear-gradient(
            90deg,
            #31e7ff,
            #9b77ff
        );
        opacity: 0.8;
    }

    .dna-dot {
        position: absolute;
        top: -5px;
        width: 17px;
        height: 17px;
        border-radius: 50%;
        box-shadow: 0 0 14px currentColor;
    }

    .dna-dot.left {
        left: 0;
        color: #30dfff;
        background: #30dfff;
    }

    .dna-dot.right {
        right: 0;
        color: #a078ff;
        background: #a078ff;
    }

    @keyframes helix-spin {
        from {
            transform: rotateY(0deg);
        }
        to {
            transform: rotateY(360deg);
        }
    }

    @media (max-width: 800px) {
        .dna-scene {
            opacity: 0.2;
            right: -20px;
        }

        .hero-card {
            padding: 24px;
        }

        .hero-title {
            font-size: 1.9rem;
        }
    }

    @media (prefers-reduced-motion: reduce) {
        .dna-helix {
            animation: none;
        }
    }
    div[data-testid="stTabs"] button[role="tab"] {color:#bcd4ed;font-weight:650;padding:14px 18px;}
    div[data-testid="stTabs"] button[aria-selected="true"] {color:#5de4ff;border-bottom-color:#5de4ff;}
    div[data-testid="stExpander"] {border:1px solid rgba(93,228,255,.2);border-radius:14px;}
    div[data-testid="stMetric"] {box-shadow:0 10px 26px rgba(0,0,0,.15);}
    .research-panel {background:linear-gradient(110deg,#112a48,#0a192d);border:1px solid #285879;border-radius:16px;padding:20px 24px;margin:12px 0 24px;}
</style>
""", unsafe_allow_html=True)

# --------------------------------------------------
# FILE PATHS
# --------------------------------------------------

project = Path(__file__).resolve().parent

viral_path = (
    project
    / "src/config_code/hpv_step2/targets"
    / "HPV16_E6starI_E7_spliced_region.rna.fasta"
)

human_path = project / "results/lncrna_checks.csv"

# --------------------------------------------------
# DNA ANIMATION
# --------------------------------------------------

dna_pairs = ""

for i in range(22):
    angle = i * 32
    height = (i - 11) * 10

    dna_pairs += (
        f'<div class="dna-pair" '
        f'style="transform: translateY({height}px) rotateY({angle}deg);">'
        '<div class="dna-bar"></div>'
        '<div class="dna-dot left"></div>'
        '<div class="dna-dot right"></div>'
        '</div>'
    )


# --------------------------------------------------
# HERO SECTION
# --------------------------------------------------

hero_html = f"""
<div class="hero-card">
    <div class="hero-tag">BIOINFORMATICS RESEARCH</div>
    <div class="hero-title">
        HPV RNA<br>
        Complementarity Analysis
    </div>
    <div class="hero-subtitle">
        Explore exact reverse-complementary sequence
        matches between HPV16 RNA and human long
        non-coding RNA transcripts.
    </div>
    <div class="dna-scene">
        <div class="dna-helix">
            {dna_pairs}
        </div>
    </div>
</div>
"""

import textwrap

st.markdown(
    textwrap.dedent(hero_html).strip(),
    unsafe_allow_html=True
)


# --------------------------------------------------
# DATA SOURCES AND ANALYSIS UTILITIES
# --------------------------------------------------

RESULT_COLUMNS = [
    "Human Transcript", "Viral Target", "Human Start", "Human End",
    "Viral Start", "Viral End", "Viral Sequence", "Reverse Complement",
    "Human Matched Sequence", "Match Length", "Mismatches", "Gaps"
]
SAVED_RESULTS = project / "results" / "exact_complementarity_matches.csv"
MAX_PREVIEW = 10000
MAX_DISPLAY_IDS = 1000


def normalize_rna(sequence):
    return str(sequence).strip().upper().replace("T", "U")


@st.cache_data(show_spinner=False)
def read_viral_target(path, modified_time):
    record = SeqIO.read(path, "fasta")
    return record.id, normalize_rna(record.seq)


def reverse_complement_rna(sequence):
    return str(Seq(sequence).reverse_complement_rna())


def build_target_index(viral_sequence, k):
    """Index complementary human k-mers by all matching viral positions."""
    lookup = defaultdict(list)
    for start in range(len(viral_sequence) - k + 1):
        viral_kmer = viral_sequence[start:start + k]
        if set(viral_kmer) <= set("ACGU"):
            lookup[reverse_complement_rna(viral_kmer)].append((start, viral_kmer))
    return lookup


def scan_transcript(seq_id, human_sequence, target_id, k, target_index):
    """Yield all exact matches, including overlapping matches and repeated viral sites."""
    for human_start in range(len(human_sequence) - k + 1):
        human_kmer = human_sequence[human_start:human_start + k]
        for viral_start, viral_kmer in target_index.get(human_kmer, ()):
            yield {
                "Human Transcript": seq_id,
                "Viral Target": target_id,
                "Human Start": human_start,
                "Human End": human_start + k,
                "Viral Start": viral_start,
                "Viral End": viral_start + k,
                "Viral Sequence": viral_kmer,
                "Reverse Complement": human_kmer,
                "Human Matched Sequence": human_kmer,
                "Match Length": k,
                "Mismatches": 0,
                "Gaps": 0,
            }


def run_search(target_id, viral_sequence, human_csv, k, transcript_limit, progress, status):
    target_index = build_target_index(viral_sequence, k)
    if not target_index:
        raise ValueError("The viral target contains no valid A/C/G/U windows of this length.")

    output = tempfile.NamedTemporaryFile(
        mode="w", suffix=".csv", prefix="hpv_matches_", newline="", delete=False,
        encoding="utf-8"
    )
    output_path = Path(output.name)
    preview = []
    transcript_counts = Counter()
    position_counts = Counter()
    total_matches = 0
    analyzed = 0
    valid_transcripts = 0
    started_at = datetime.now(timezone.utc).isoformat()
    try:
        writer = csv.DictWriter(output, fieldnames=RESULT_COLUMNS)
        writer.writeheader()
        # Chunking keeps the entire input transcript dataset out of memory.
        reader = pd.read_csv(
            human_csv, usecols=["seq_id", "seq"], dtype=str,
            chunksize=500, keep_default_na=False
        )
        for chunk in reader:
            for row in chunk.itertuples(index=False):
                if analyzed >= transcript_limit:
                    break
                analyzed += 1
                seq_id = str(row.seq_id).strip()
                human_sequence = normalize_rna(row.seq)
                if seq_id and human_sequence:
                    valid_transcripts += 1
                    for match in scan_transcript(
                        seq_id, human_sequence, target_id, k, target_index
                    ):
                        writer.writerow(match)
                        total_matches += 1
                        transcript_counts[seq_id] += 1
                        position_counts[match["Viral Start"]] += 1
                        if len(preview) < MAX_PREVIEW:
                            preview.append(match)
                if analyzed % 500 == 0 or analyzed == transcript_limit:
                    progress.progress(min(analyzed / transcript_limit, 1.0))
                    status.caption(f"Processed {analyzed:,} transcript records")
            if analyzed >= transcript_limit:
                break
        output.flush()
    except Exception:
        output.close()
        output_path.unlink(missing_ok=True)
        raise
    finally:
        output.close()

    return {
        "preview": pd.DataFrame(preview, columns=RESULT_COLUMNS),
        "path": str(output_path),
        "total_matches": total_matches,
        "transcript_counts": dict(transcript_counts),
        "position_counts": dict(position_counts),
        "unique_transcripts": len(transcript_counts),
        "analyzed": analyzed,
        "valid_transcripts": valid_transcripts,
        "k": k,
        "target_id": target_id,
        "viral_length": len(viral_sequence),
        "source": "new",
        "started_at_utc": started_at,
        "input_viral_file": str(viral_path.relative_to(project)),
        "input_human_file": str(human_csv.relative_to(project)),
        "python_version": platform.python_version(),
        "packages": {name: version(name) for name in ("streamlit", "pandas", "biopython")},
    }


@st.cache_data(show_spinner=False)
def saved_results_preview(path, modified_time):
    frame = pd.read_csv(path, nrows=MAX_PREVIEW)
    # Count rows without loading the entire CSV into memory.
    with open(path, "r", encoding="utf-8", newline="") as handle:
        row_count = max(sum(1 for _ in handle) - 1, 0)
    return frame, row_count


def alignment_html(viral_seq, human_seq):
    """Viral 5'→3' aligned to antiparallel human 3'→5'."""
    colors = {"A": "#147d9b", "U": "#7058b9", "G": "#177c67", "C": "#b37a28"}
    human_antiparallel = human_seq[::-1]
    def render(seq):
        return "".join(
            f'<span style="display:inline-block;min-width:24px;padding:6px 3px;'
            f'margin:2px;text-align:center;border-radius:5px;background:{colors.get(base, "#555")};'
            f'color:white;font-weight:700">{base}</span>'
            for base in seq
        )
    return (
        '<div style="background:#0b1b31;padding:18px;border-radius:12px;overflow-x:auto;'
        'font-family:monospace;color:#e8f2ff;white-space:nowrap">'
        f'<div>Viral 5′ &nbsp; {render(viral_seq)} &nbsp; 3′</div>'
        f'<div style="padding:7px 0;color:#7de0e8">Human 3′ {render(human_antiparallel)} &nbsp; 5′</div>'
        '</div>'
    )


def results_dashboard(frame, total, counts=None, positions=None, viral_length=None, source_label="Results"):
    st.subheader(source_label)
    if frame is None:
        st.info("Run an analysis to explore results.")
        return
    if "Human Transcript" not in frame.columns:
        st.error("The results CSV is missing the 'Human Transcript' column.")
        return

    unique_count = len(counts) if counts is not None else frame["Human Transcript"].nunique()
    metric1, metric2, metric3, metric4 = st.columns(4)
    metric1.metric("Exact Matches", f"{total:,}")
    metric2.metric("Unique Transcripts" + (" (preview)" if counts is None and total > len(frame) else ""), f"{unique_count:,}")
    k_values = frame["Match Length"].dropna().unique() if "Match Length" in frame else []
    metric3.metric("Match Length", f"{int(k_values[0])} nt" if len(k_values) == 1 else "Mixed / unknown")
    metric4.metric("Rows in Preview", f"{len(frame):,}")
    if total > len(frame):
        st.caption(f"Showing the first {len(frame):,} of {total:,} rows. Charts and transcript counts "
                   + ("use the complete new run." if counts is not None else "use the preview only for saved results."))
    if frame.empty:
        st.info("No exact matches were found for these settings.")
        return

    st.markdown("#### Explore matching transcripts")
    if total > len(frame):
        st.info("Search and filters below apply to the preview only, not the entire CSV. Download the full results to inspect all matches.")
    q = st.text_input("Search transcript ID", placeholder="e.g., ENST000...", key=f"q_{source_label}")
    available_ids = sorted(frame["Human Transcript"].dropna().astype(str).unique().tolist())[:MAX_DISPLAY_IDS]
    if frame["Human Transcript"].nunique() > MAX_DISPLAY_IDS:
        st.caption("Transcript selector shows the first 1,000 IDs; use search to find other IDs in the preview.")
    chosen = st.selectbox("Select a transcript", ["All transcripts"] + available_ids,
                          key=f"transcript_{source_label}")
    filtered = frame
    if q:
        filtered = filtered[filtered["Human Transcript"].astype(str).str.contains(q, case=False, regex=False, na=False)]
    if chosen != "All transcripts":
        filtered = filtered[filtered["Human Transcript"].astype(str) == chosen]
    if "Viral Start" in filtered.columns and not filtered.empty:
        min_pos = int(pd.to_numeric(filtered["Viral Start"], errors="coerce").min())
        max_pos = int(pd.to_numeric(filtered["Viral Start"], errors="coerce").max())
        if min_pos < max_pos:
            interval = st.slider("Viral start coordinate (0-based)", min_pos, max_pos,
                                 (min_pos, max_pos), key=f"pos_{source_label}")
            filtered = filtered[pd.to_numeric(filtered["Viral Start"], errors="coerce").between(*interval)]
    st.caption(f"{len(filtered):,} matching rows in the current preview")
    st.dataframe(filtered, use_container_width=True, hide_index=True, height=360)

    st.markdown("#### Sequence alignment viewer")
    if filtered.empty:
        st.info("Adjust the filters to select a match.")
    elif all(c in filtered for c in ("Viral Sequence", "Human Matched Sequence")):
        # Selection by preview row, rather than huge selectbox options.
        # Do not persist a stale row index when filtering shrinks the table.
        row_number = st.number_input("Match row (within filtered results)", min_value=1,
                                     max_value=len(filtered), value=1, step=1)
        selected = filtered.iloc[int(row_number) - 1]
        viral = normalize_rna(selected["Viral Sequence"])
        human = normalize_rna(selected["Human Matched Sequence"])
        if viral and human and len(viral) == len(human):
            st.markdown(alignment_html(viral, human), unsafe_allow_html=True)
            st.caption(
                f"Transcript: {selected['Human Transcript']} | "
                f"Human coordinates: {selected.get('Human Start', '?')}–{selected.get('Human End', '?')} | "
                f"Viral coordinates: {selected.get('Viral Start', '?')}–{selected.get('Viral End', '?')}"
            )
            st.code(f"Viral 5′→3′: {viral}\nHuman 3′→5′: {human[::-1]}", language="text")
            if reverse_complement_rna(viral) != human:
                st.warning("These selected sequences are not exact reverse complements. Check the source CSV.")
        else:
            st.warning("Sequence data is missing or invalid for this row.")
    else:
        st.info("The selected CSV does not contain sequence columns for alignment.")

    st.markdown("#### Match distributions")
    left, right = st.columns(2)
    with left:
        st.write("**Top 15 matching transcripts**")
        if counts is not None:
            top = pd.Series(counts).sort_values(ascending=False).head(15)
        else:
            top = frame.groupby("Human Transcript").size().sort_values(ascending=False).head(15)
        st.bar_chart(top, horizontal=True)
    with right:
        st.write("**Matches by viral start position**")
        if positions is not None:
            pos_series = pd.Series(positions, dtype="int64")
        elif "Viral Start" in frame:
            pos_series = pd.to_numeric(frame["Viral Start"], errors="coerce").dropna().astype(int).value_counts()
        else:
            pos_series = pd.Series(dtype="int64")
        if not pos_series.empty:
            pos_series.index = pos_series.index.astype(int)
            max_start = max(int(pos_series.index.max()), (viral_length or 1) - 1)
            bin_size = max(1, (max_start + 1 + 19) // 20)
            bins = defaultdict(int)
            for pos, count in pos_series.items():
                start = (int(pos) // bin_size) * bin_size
                bins[f"{start}–{start + bin_size - 1}"] += int(count)
            chart_df = pd.DataFrame(
                [(int(label.split("–")[0]), label, count) for label, count in bins.items()],
                columns=["sort", "Position", "Matches"]
            ).sort_values("sort")
            st.bar_chart(chart_df.set_index("Position")["Matches"])
        else:
            st.info("No viral positions available.")

    st.markdown("#### Match length and transcript coverage")
    chart_left, chart_right = st.columns(2)
    with chart_left:
        st.write("**Matches by length (available rows)**")
        if "Match Length" in frame:
            lengths = pd.to_numeric(frame["Match Length"], errors="coerce").dropna().astype(int)
            st.bar_chart(lengths.value_counts().sort_index())
            st.caption("For a new run, all matches use the selected k-mer length. "
                       "Compare multiple runs to study different lengths.")
    with chart_right:
        st.write("**Covered nucleotides in selected transcript (preview)**")
        if all(c in filtered for c in ("Human Start", "Human End")) and not filtered.empty:
            first_id = str(filtered.iloc[0]["Human Transcript"])
            rows = filtered[filtered["Human Transcript"].astype(str) == first_id]
            intervals = sorted((int(a), int(b)) for a, b in zip(rows["Human Start"], rows["Human End"]))
            covered = 0
            if intervals:
                left_edge, right_edge = intervals[0]
                for start, end in intervals[1:]:
                    if start <= right_edge:
                        right_edge = max(right_edge, end)
                    else:
                        covered += right_edge - left_edge
                        left_edge, right_edge = start, end
                covered += right_edge - left_edge
            st.metric("Distinct matched bases", f"{covered:,} nt")
            st.caption(f"Transcript {first_id}. Counts overlapping positions once; "
                       "this is not percent coverage without transcript length.")
        else:
            st.info("Select a transcript with position columns to see covered bases.")

    st.markdown("#### Export results")
    st.download_button("Download filtered preview CSV", filtered.to_csv(index=False).encode("utf-8"),
                       "hpv_filtered_preview.csv", "text/csv", key=f"filtered_dl_{source_label}")


# --------------------------------------------------
# NAVIGATION
# --------------------------------------------------

if "run" not in st.session_state:
    st.session_state.run = None

overview_tab, run_tab, results_tab, methodology_tab = st.tabs(
    ["🏠 Overview", "▶ Run Analysis", "📊 Results", "📖 Methodology"]
)

with overview_tab:
    st.header("About the Project")
    st.write(
        "This research tool investigates exact reverse-complementary RNA sequence "
        "matches between the HPV16 E6*I/E7 spliced target and human long non-coding RNA (lncRNA) transcripts."
    )
    st.write(
        "Use the analysis tab to search sequences, then explore transcript IDs, "
        "viral positions, and antiparallel alignments in the results tab."
    )
    st.info(
        "An exact sequence match does not establish physical RNA binding, "
        "gene suppression, biological activity, or therapeutic effectiveness."
    )
    st.markdown('<div class="research-panel"><strong style="color:#5de4ff">ANALYSIS PIPELINE</strong><h3 style="color:white">Viral RNA → k-mer index → human lncRNA matches</h3><p>Inspect the position distributions and antiparallel alignments for every selected match.</p></div>', unsafe_allow_html=True)
    st.subheader("Workflow")
    a, b, c = st.columns(3)
    a.markdown("**1. Select**  \nChoose a viral target, k-mer length, and transcript count.")
    b.markdown("**2. Analyze**  \nScan human lncRNAs for reverse-complementary windows.")
    c.markdown("**3. Explore**  \nFilter matches, inspect alignments, visualize distributions, and export.")
    st.subheader("Data availability")
    st.write(f"HPV16 target FASTA: {'✅ Found' if viral_path.is_file() else '❌ Missing'}")
    st.write(f"Human transcript CSV: {'✅ Found' if human_path.is_file() else '❌ Missing'}")
    st.write(f"Previously saved matches: {'✅ Found' if SAVED_RESULTS.is_file() else 'Not found'}")

with run_tab:
    st.header("Run Exact Complementarity Analysis")
    st.write("**HPV target:** HPV16 E6*I/E7 spliced RNA")
    col1, col2 = st.columns(2)
    with col1:
        match_length = st.slider("Match length (nucleotides)", 8, 20, 10)
    with col2:
        transcript_limit = st.selectbox(
            "Human transcript records to analyze",
            options=[100, 500, 1000, 10000, 50000, 197195],
            format_func=lambda n: f"Up to {n:,} records"
        )
    st.caption("Longer exact matches are typically less common. The transcript limit is an upper bound. "
               "Each analysis is an exact-match screen, not an RNA-binding prediction.")
    if st.button("▶ Run Analysis", type="primary", use_container_width=True):
        if not viral_path.is_file():
            st.error(f"HPV target FASTA not found: {viral_path}")
        elif not human_path.is_file():
            st.error(f"Human transcript CSV not found: {human_path}. Run the sequence export first.")
        else:
            progress = st.progress(0)
            status = st.empty()
            try:
                target_id, target_sequence = read_viral_target(
                    str(viral_path), viral_path.stat().st_mtime
                )
                previous = st.session_state.run
                with st.spinner("Searching exact complementary sequences..."):
                    new_run = run_search(
                        target_id, target_sequence, human_path, match_length,
                        transcript_limit, progress, status
                    )
                st.session_state.run = new_run
                if previous and previous.get("source") == "new":
                    Path(previous["path"]).unlink(missing_ok=True)
                st.success(f"Finished: {new_run['total_matches']:,} matches across "
                           f"{new_run['analyzed']:,} transcript records.")
                st.info("Open the Results tab to explore the new analysis.")
            except Exception as exc:
                st.error(f"Analysis failed: {exc}")
            finally:
                progress.empty()
                status.empty()

with results_tab:
    st.header("Analysis Results")
    if st.session_state.run is not None:
        run = st.session_state.run
        st.caption(f"Current run: {run['target_id']} | {run['k']} nt | "
                   f"{run['analyzed']:,} records analyzed")
        results_dashboard(run["preview"], run["total_matches"],
                          run["transcript_counts"], run["position_counts"],
                          run["viral_length"], "Current analysis")
        if Path(run["path"]).is_file():
            size_mb = Path(run["path"]).stat().st_size / (1024 * 1024)
            if size_mb <= 100:
                with open(run["path"], "rb") as handle:
                    st.download_button("Download complete current analysis CSV", handle.read(),
                                       "hpv_full_analysis.csv", "text/csv", key="full_current")
            else:
                st.warning(f"Full results are {size_mb:,.1f} MB; to avoid exhausting browser memory, "
                           f"the CSV is stored on the server at {run['path']}.")
        st.download_button("Download run metadata (JSON)",
                           json.dumps({k: v for k, v in run.items() if k not in
                                       ("preview", "transcript_counts", "position_counts")}, indent=2),
                           "hpv_analysis_metadata.json", "application/json")
    else:
        st.info("No analysis has been run in this session yet. Open Run Analysis to begin.")

    if SAVED_RESULTS.is_file():
        with st.expander("Previously saved analysis results", expanded=st.session_state.run is None):
            try:
                saved_frame, saved_total = saved_results_preview(
                    str(SAVED_RESULTS), SAVED_RESULTS.stat().st_mtime
                )
                results_dashboard(saved_frame, saved_total, source_label="Saved analysis")
                if SAVED_RESULTS.stat().st_size <= 100 * 1024 * 1024:
                    with open(SAVED_RESULTS, "rb") as handle:
                        st.download_button("Download complete saved analysis CSV", handle.read(),
                                           "hpv_saved_analysis.csv", "text/csv", key="full_saved")
                else:
                    st.caption("Saved CSV exceeds 100 MB; access it directly from the results folder.")
            except Exception as exc:
                st.error(f"Could not load saved results: {exc}")

with methodology_tab:
    st.header("Methodology & Scientific Limitations")
    st.markdown('<div class="research-panel"><strong style="color:#5de4ff">RESEARCH PROTOCOL</strong><h3 style="color:white">From sequences to inspectable matches</h3><p>Data provenance, matching assumptions, and reproducibility.</p></div>', unsafe_allow_html=True)
    st.subheader("Input datasets")
    st.write("**Viral target:** HPV16 E6*I/E7 spliced RNA sequence in FASTA format.")
    st.code(str(viral_path.relative_to(project)), language="text")
    (st.success if viral_path.is_file() else st.error)("Viral FASTA found" if viral_path.is_file() else "Viral FASTA missing")
    st.write("**Human input:** CSV containing `seq_id` and `seq` columns for human lncRNA transcripts.")
    st.code(str(human_path.relative_to(project)), language="text")
    (st.success if human_path.is_file() else st.error)("Human CSV found" if human_path.is_file() else "Human CSV missing")
    st.caption("Confirm reference assembly, GENCODE release, sequence accession, and target provenance before publishing results.")
    st.subheader("Search algorithm")
    st.markdown(
        "1. Convert input sequences to uppercase RNA alphabet (`T` → `U`).\n"
        "2. Generate every valid A/C/G/U k-mer window of the HPV16 target.\n"
        "3. Build an index from each reverse-complementary k-mer to its viral positions.\n"
        "4. Read human transcripts in chunks and scan each possible k-mer position.\n"
        "5. Report all exact matches, including overlapping positions and repeated target sites."
    )
    st.subheader("Coordinate conventions")
    st.write("All start coordinates are **zero-based** and all end coordinates are **exclusive**. "
             "For example, 0–10 represents the first 10 nucleotides. The human strand is "
             "shown antiparallel (3′→5′) in the alignment viewer.")
    st.subheader("Interpretation")
    st.warning(
        "Exact complementarity is only a computational screening signal. "
        "This pipeline does not model RNA accessibility, folding, expression, "
        "subcellular localization, binding energy, off-target interactions, "
        "or experimental effects. Match counts are not evidence of therapeutic efficacy."
    )
    st.subheader("Reproducibility")
    st.write("Record the input dataset releases, target accession, analysis date, "
             "match length, transcript limit, and software versions alongside exported results.")
    st.write("**Install packages:**")
    st.code("python3 -m pip install streamlit pandas biopython", language="bash")
    st.write("**Launch website:**")
    st.code("python3 -m streamlit run app_research_dashboard.py", language="bash")
    st.caption("For a published analysis, manually supply the actual dataset releases, "
               "reference accession and provenance; they cannot be inferred from file names.")
    with st.expander("Download reproducibility metadata", expanded=True):
        accession = st.text_input("Verified viral accession (optional)")
        release = st.text_input("Verified GENCODE release (optional)")
        assembly = st.text_input("Verified human assembly (optional)")
        run = st.session_state.run
        metadata = {
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "viral_file": str(viral_path.relative_to(project)),
            "human_file": str(human_path.relative_to(project)),
            "accession_user_supplied": accession or None,
            "gencode_user_supplied": release or None,
            "assembly_user_supplied": assembly or None,
            "analysis": {key: run[key] for key in ("target_id", "k", "analyzed", "total_matches", "started_at_utc")}
                        if run else None,
            "python": platform.python_version(),
            "packages": {name: version(name) for name in ("streamlit", "pandas", "biopython")},
            "coordinates": "0-based, end-exclusive",
            "interpretation": "Exact matches only; not evidence of binding or therapeutic efficacy"
        }
        st.download_button("Download metadata JSON", json.dumps(metadata, indent=2).encode(),
                           "hpv_analysis_metadata.json", "application/json")

st.divider()
st.markdown(
    '<div style="text-align:center;color:#8da9c9;padding:20px;font-size:13px">'
    'HPV RNA Complementarity Analysis · Bioinformatics Research Project</div>',
    unsafe_allow_html=True
)
