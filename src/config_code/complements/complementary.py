import heapq
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import pandas as pd


WINDOW_SIZE = 10
TEST_LIMIT = None  # Change to None to process all human transcripts
PROJECT_PATH = Path(__file__).resolve().parents[3]


def rev_comp(strand: str) -> str:
    return str(Seq(strand).reverse_complement_rna())


def read_fasta(path: Path) -> tuple[str, str]:
    record = SeqIO.read(path, "fasta")
    sequence = str(record.seq).upper().replace("T", "U")

    return record.id, sequence


def normalize_sequence(sequence: str) -> str:
    return sequence.strip().upper().replace("T", "U")


def load_human_sequences(path: Path, test_limit=None) -> pd.DataFrame:
    human_data = pd.read_csv(path, usecols=["seq_id", "seq"], dtype=str)
    input_rows = len(human_data)
    human_data = human_data.dropna(subset=["seq_id", "seq"]).copy()
    human_data["seq"] = human_data["seq"].map(normalize_sequence)
    human_data = human_data[human_data["seq"].str.len() > 0]
    valid_rows = len(human_data)

    if test_limit is not None:
        if test_limit <= 0:
            raise ValueError("test_limit must be a positive integer or None.")
        human_data = human_data.head(test_limit)

    human_data = human_data.reset_index(drop=True)
    human_data.attrs.update(
        input_rows=input_rows,
        valid_rows_before_limit=valid_rows,
        excluded_rows=input_rows - valid_rows,
        test_limit=test_limit,
    )
    return human_data


def prepare_windows(sequence: str, window_size=10, window_step=1) -> list[tuple]:
    if window_size <= 0 or window_step <= 0:
        raise ValueError("Window size and step must be positive integers.")

    sequence = normalize_sequence(sequence)
    windows = []
    for start in range(0, len(sequence) - window_size + 1, window_step):
        end = start + window_size
        chunk = sequence[start:end]
        if all(base in "ACGU" for base in chunk):
            windows.append((start, end, chunk, rev_comp(chunk)))
    return windows


def build_match_index(human_sequences, queries, window_size=10) -> dict:
    """Index requested RNA words, including every overlapping occurrence.

    Each value has two uint32 columns containing the human row and start.
    A separator between transcripts prevents matches crossing their boundaries.
    """
    if not 1 <= window_size <= 15:
        raise ValueError("The compact index supports window sizes from 1 to 15.")

    base_codes = np.full(256, 255, dtype=np.uint8)
    for code, base in enumerate(b"ACGU"):
        base_codes[base] = code

    query_by_code = {}
    for query in queries:
        if len(query) != window_size or any(base not in "ACGU" for base in query):
            raise ValueError("Every indexed query must be an unambiguous RNA window.")
        code = 0
        for base in query:
            code = (code << 2) | int(base_codes[ord(base)])
        query_by_code[code] = query

    index = {
        query: np.empty((0, 2), dtype=np.uint32)
        for query in query_by_code.values()
    }
    if not query_by_code:
        return index

    query_codes = np.array(list(query_by_code), dtype=np.uint32)
    wanted = None
    if window_size <= 10:
        wanted = np.zeros(1 << (2 * window_size), dtype=bool)
        wanted[query_codes] = True

    hit_code_batches = []
    hit_position_batches = []
    batch_sequences = []
    batch_starts = []
    batch_rows = []
    batch_length = 0

    def scan_batch():
        joined = "N".join(batch_sequences).encode("ascii", errors="replace")
        bases = base_codes[np.frombuffer(joined, dtype=np.uint8)]
        window_count = len(bases) - window_size + 1
        if window_count <= 0:
            return

        codes = np.zeros(window_count, dtype=np.uint32)
        valid = np.ones(window_count, dtype=bool)
        for offset in range(window_size):
            current = bases[offset:offset + window_count]
            codes <<= 2
            codes |= current & 3
            valid &= current < 4

        if wanted is None:
            selected = valid & np.isin(codes, query_codes, kind="sort")
        else:
            selected = valid & wanted[codes]

        positions = np.flatnonzero(selected)
        if not len(positions):
            return

        starts = np.asarray(batch_starts, dtype=np.int64)
        rows = np.asarray(batch_rows, dtype=np.uint32)
        local_rows = np.searchsorted(starts, positions, side="right") - 1
        hit_pairs = np.empty((len(positions), 2), dtype=np.uint32)
        hit_pairs[:, 0] = rows[local_rows]
        hit_pairs[:, 1] = positions - starts[local_rows]
        hit_code_batches.append(codes[positions])
        hit_position_batches.append(hit_pairs)

    for row_index, sequence in enumerate(human_sequences):
        if row_index >= 2 ** 32 or len(sequence) >= 2 ** 32:
            raise ValueError("Input coordinates exceed the uint32 index capacity.")
        batch_sequences.append(sequence)
        batch_starts.append(batch_length)
        batch_rows.append(row_index)
        batch_length += len(sequence) + 1

        if batch_length >= 2_000_000:
            scan_batch()
            batch_sequences.clear()
            batch_starts.clear()
            batch_rows.clear()
            batch_length = 0

        if (row_index + 1) % 10_000 == 0:
            print(f"Indexed human sequences: {row_index + 1}", flush=True)

    if batch_sequences:
        scan_batch()
    if not hit_code_batches:
        return index

    all_codes = np.concatenate(hit_code_batches)
    all_positions = np.concatenate(hit_position_batches)
    order = np.argsort(all_codes, kind="stable")
    sorted_codes = all_codes[order]
    sorted_positions = all_positions[order]
    boundaries = np.concatenate((
        [0], np.flatnonzero(np.diff(sorted_codes)) + 1, [len(sorted_codes)]
    ))
    for start, end in zip(boundaries[:-1], boundaries[1:]):
        query = query_by_code[int(sorted_codes[start])]
        index[query] = sorted_positions[start:end]
    return index


def iter_indexed_matches(windows, index):
    """Yield matches in human-row, viral-start, human-start order."""
    def window_matches(window):
        viral_start, viral_end, chunk, comp = window
        for human_row, human_start in index[comp]:
            yield (
                int(human_row), viral_start, viral_end, chunk, comp,
                int(human_start)
            )

    streams = [window_matches(window) for window in windows]
    yield from heapq.merge(*streams, key=lambda match: (match[0], match[1], match[5]))


def find_matches(query: str, human_sequence: str) -> list[int]:
    positions = []
    start = 0

    while True:
        position = human_sequence.find(query, start)

        if position == -1:
            break

        positions.append(position)

        # Move one base forward to include overlapping occurrences.
        start = position + 1

    return positions


def main():
    project = PROJECT_PATH

    viral_path = (
        project
        / "src/config_code/hpv_step2/targets"
        / "HPV16_E6starI_E7_spliced_region.rna.fasta"
    )
    human_path = project / "results/lncrna_checks.csv"
    output_path = project / "results/exact_complementarity_matches.csv"

    viral_id, viral_sequence = read_fasta(viral_path)

    human_rna_csv = load_human_sequences(human_path, TEST_LIMIT)
    windows = prepare_windows(viral_sequence, WINDOW_SIZE)
    index = build_match_index(
        human_rna_csv["seq"], (window[3] for window in windows), WINDOW_SIZE
    )

    columns = [
        "hum_trspt_id",
        "viral_rna_id",
        "human_start",
        "human_end",
        "viral_start",
        "viral_end",
        "viral_window",
        "viral_reverse_complement",
        "human_matched_sequence",
        "matches",
        "match_length",
        "mismatches",
        "gaps"
    ]

    rows = []
    human_ids = human_rna_csv["seq_id"].tolist()
    for human_row, viral_start, viral_end, chunk, comp, human_start in (
        iter_indexed_matches(windows, index)
    ):
        human_end = human_start + WINDOW_SIZE
        rows.append({
            "hum_trspt_id": human_ids[human_row],
            "viral_rna_id": viral_id,
            "human_start": human_start,
            "human_end": human_end,
            "viral_start": viral_start,
            "viral_end": viral_end,
            "viral_window": chunk,
            "viral_reverse_complement": comp,
            "human_matched_sequence": comp,
            "matches": WINDOW_SIZE,
            "match_length": WINDOW_SIZE,
            "mismatches": 0,
            "gaps": 0
        })

    sequences = pd.DataFrame(rows, columns=columns)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    sequences.to_csv(output_path, index=False)

    print(sequences.head().to_string(index=False))
    print(f"\nTotal exact matches: {len(sequences)}")
    print(f"Saved to: {output_path}")

    return sequences


if __name__ == "__main__":
    results = main()
