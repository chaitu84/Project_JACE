#!/usr/bin/env bash
# Run this copy from a prepared exploratory bundle using IntaRNA 3.4.1.
set -Eeuo pipefail

if (( $# != 0 )); then
    printf 'Usage: bash %s\n' "$0" >&2
    exit 2
fi

bundle_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
intarna_bin="${INTARNA_BIN:-IntaRNA}"
python_bin="${PYTHON_BIN:-python3}"
if ! command -v -- "$python_bin" >/dev/null 2>&1; then
    printf 'Python executable unavailable: %s\n' "$python_bin" >&2
    exit 127
fi

mkdir -p -- "$bundle_dir/predictions"
run_dir="$(mktemp -d "$bundle_dir/predictions/run_$(date -u +%Y%m%dT%H%M%SZ)_XXXXXX")"
total=0
finalized=0
final_state=failed
printf 'representative_id\tstatus\tintarna_exit_code\tprediction_rows\ttranscript_ids\n' > "$run_dir/job_status.tsv"

write_status() {
    "$python_bin" - "$run_dir" "$1" "$total" "$SECONDS" <<'PY'
import csv
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

run_dir, state, total, elapsed = sys.argv[1:]
run_dir = Path(run_dir)
with (run_dir / "job_status.tsv").open(newline="", encoding="utf-8") as source:
    jobs = list(csv.DictReader(source, delimiter="\t"))
report = {
    "status": state,
    "updated_utc": datetime.now(timezone.utc).isoformat(),
    "expected_intarna_version": "3.4.1",
    "requested_distinct_sequences": int(total),
    "completed_with_interaction": sum(job["status"] == "completed" for job in jobs),
    "completed_without_interaction": sum(job["status"] == "completed_no_interaction" for job in jobs),
    "failed_jobs": sum(job["status"].startswith("failed") for job in jobs),
    "elapsed_seconds": int(elapsed),
    "evidence_label": "exploratory",
    "statistical_context": "None of the shortlisted transcripts passed BH FDR 0.05",
    "jobs": jobs,
}
temporary = run_dir / "run_status.json.tmp"
temporary.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
temporary.replace(run_dir / "run_status.json")
PY
}

on_exit() {
    local exit_code=$?
    if (( finalized == 0 )); then
        write_status "$final_state" || true
        printf 'Run ended with status %s. Inspect %s\n' "$final_state" "$run_dir" >&2
    fi
    return "$exit_code"
}
trap on_exit EXIT
trap 'final_state=interrupted; exit 130' INT
trap 'final_state=interrupted; exit 143' TERM
write_status running
printf 'Exploratory IntaRNA run directory: %s\n' "$run_dir"

"$python_bin" - "$bundle_dir" "$run_dir/prediction_inputs.tsv" 2> "$run_dir/input_validation.log" <<'PY'
import csv
import hashlib
import re
import sys
from pathlib import Path

bundle = Path(sys.argv[1]).resolve()
manifest = bundle / "prediction_inputs.tsv"
fields = ["representative_id", "human_fasta", "viral_fasta", "sequence_sha256", "transcript_ids"]

def read_fasta(path):
    lines = path.read_text(encoding="utf-8").splitlines()
    headers = [line for line in lines if line.startswith(">")]
    if len(headers) != 1 or not lines or not lines[0].startswith(">"):
        raise ValueError(f"Expected one FASTA record in {path}")
    sequence = "".join(line.strip() for line in lines[1:])
    if not sequence or set(sequence) - set("ACGU"):
        raise ValueError(f"Expected unambiguous uppercase RNA in {path}")
    return headers[0][1:].split()[0], sequence

def input_path(relative):
    path = Path(relative)
    if path.is_absolute():
        raise ValueError("Manifest FASTA paths must be relative to the bundle")
    resolved = (bundle / path).resolve()
    try:
        resolved.relative_to(bundle)
    except ValueError:
        raise ValueError("Manifest FASTA path leaves the bundle") from None
    if not resolved.is_file():
        raise ValueError(f"Missing FASTA {resolved}")
    return resolved

with manifest.open(newline="", encoding="utf-8") as source:
    reader = csv.DictReader(source, delimiter="\t")
    if reader.fieldnames != fields:
        raise ValueError("Unexpected prediction_inputs.tsv header")
    rows = list(reader)
if not rows:
    raise ValueError("The prediction manifest contains no distinct sequences")
seen_ids = set()
seen_sequences = set()
for row in rows:
    if set(row) != set(fields) or any(not value or any(char in value for char in "\t\r\n") for value in row.values()):
        raise ValueError("Every manifest field must contain one nonempty TSV value")
    representative = row["representative_id"]
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", representative) or representative in seen_ids:
        raise ValueError("Representative IDs must be unique and safe as filenames")
    seen_ids.add(representative)
    human_id, sequence = read_fasta(input_path(row["human_fasta"]))
    read_fasta(input_path(row["viral_fasta"]))
    digest = hashlib.sha256(sequence.encode("ascii")).hexdigest()
    if human_id != representative or digest != row["sequence_sha256"]:
        raise ValueError(f"Human FASTA identity or checksum differs for {representative}")
    if sequence in seen_sequences:
        raise ValueError("The prediction manifest contains duplicate full human sequences")
    seen_sequences.add(sequence)
with Path(sys.argv[2]).open("x", newline="", encoding="utf-8") as destination:
    writer = csv.DictWriter(destination, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
PY
total=$(( $(wc -l < "$run_dir/prediction_inputs.tsv") - 1 ))
write_status running

if ! command -v -- "$intarna_bin" >/dev/null 2>&1; then
    printf 'IntaRNA executable unavailable: %s\n' "$intarna_bin" >&2
    exit 127
fi
if ! "$intarna_bin" --version > "$run_dir/intarna_version.log" 2>&1; then
    printf 'IntaRNA version check failed. Inspect %s\n' "$run_dir/intarna_version.log" >&2
    exit 1
fi
if ! "$python_bin" - "$run_dir/intarna_version.log" 2> "$run_dir/version_validation.log" <<'PY'
import re
import sys
from pathlib import Path

text = Path(sys.argv[1]).read_text(encoding="utf-8", errors="replace")
match = re.search(r"IntaRNA\b[^\r\n0-9]*([0-9]+\.[0-9]+\.[0-9]+)", text, re.IGNORECASE)
if match is None or match.group(1) != "3.4.1":
    raise ValueError("This runner requires IntaRNA 3.4.1")
PY
then
    printf 'Expected IntaRNA 3.4.1. Inspect %s\n' "$run_dir/intarna_version.log" >&2
    exit 1
fi

number=0
failed=0
while IFS=$'\t' read -r representative human_relative viral_relative sequence_hash transcript_ids; do
    number=$((number + 1))
    human_path="$bundle_dir/$human_relative"
    viral_path="$bundle_dir/$viral_relative"
    prediction="$run_dir/${representative}.csv"
    partial="$prediction.partial"
    command=(
        "$intarna_bin" "--target=$human_path" "--query=$viral_path"
        --model=X --mode=H --acc=C --tAccW=150 --tAccL=100
        --qAccW=150 --qAccL=100 --seedBP=7 --outNumber=1 --outMode=C
        --outCsvCols=id1,start1,end1,id2,start2,end2,E,ED1,ED2 "--out=$partial"
    )
    printf '%q ' "${command[@]}" > "$run_dir/${representative}.command.sh"
    printf '\n' >> "$run_dir/${representative}.command.sh"
    printf '[%d/%d] Starting %s (%s)\n' "$number" "$total" "$representative" "$transcript_ids"
    if "${command[@]}" > "$run_dir/${representative}.stdout.log" 2> "$run_dir/${representative}.stderr.log"; then
        if prediction_rows=$("$python_bin" - "$partial" "$human_path" "$viral_path" 2> "$run_dir/${representative}.validation.log" <<'PY'
import csv
import math
import sys
from pathlib import Path

def fasta_identity(path):
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    return lines[0][1:].split()[0], len("".join(line.strip() for line in lines[1:]))

human_id, human_length = fasta_identity(sys.argv[2])
viral_id, viral_length = fasta_identity(sys.argv[3])
with Path(sys.argv[1]).open(newline="", encoding="utf-8-sig") as source:
    reader = csv.DictReader(source, delimiter=";")
    expected = "id1,start1,end1,id2,start2,end2,E,ED1,ED2".split(",")
    if reader.fieldnames != expected:
        raise ValueError("Prediction CSV has an unexpected or missing header")
    rows = list(reader)
if len(rows) > 1:
    raise ValueError("Expected at most one interaction with outNumber=1")
for row in rows:
    if set(row) != set(expected) or any(value is None for value in row.values()):
        raise ValueError("Prediction CSV row has missing or extra fields")
    if row["id1"] != human_id or row["id2"] != viral_id:
        raise ValueError("Prediction CSV IDs differ from the full input FASTAs")
    start1, end1, start2, end2 = [int(row[key]) for key in ("start1", "end1", "start2", "end2")]
    if not (1 <= start1 <= end1 <= human_length and 1 <= start2 <= end2 <= viral_length):
        raise ValueError("Prediction coordinates leave the full input sequences")
    if not all(math.isfinite(float(row[key])) for key in ("E", "ED1", "ED2")):
        raise ValueError("Prediction energy or accessibility value is nonfinite")
print(len(rows))
PY
        ); then
            mv -- "$partial" "$prediction"
            job_state=completed
            if (( prediction_rows == 0 )); then
                job_state=completed_no_interaction
            fi
            printf '%s\t%s\t0\t%s\t%s\n' "$representative" "$job_state" "$prediction_rows" "$transcript_ids" >> "$run_dir/job_status.tsv"
            printf '[%d/%d] %s: %s\n' "$number" "$total" "$representative" "$job_state"
        else
            failed=$((failed + 1))
            printf '%s\tfailed_output_validation\t0\tunknown\t%s\n' "$representative" "$transcript_ids" >> "$run_dir/job_status.tsv"
            printf '[%d/%d] %s: output validation failed\n' "$number" "$total" "$representative" >&2
        fi
    else
        intarna_exit=$?
        failed=$((failed + 1))
        printf '%s\tfailed_intarna\t%d\tunknown\t%s\n' "$representative" "$intarna_exit" "$transcript_ids" >> "$run_dir/job_status.tsv"
        printf '[%d/%d] %s: IntaRNA exited %d\n' "$number" "$total" "$representative" "$intarna_exit" >&2
    fi
    write_status running
done < <(tail -n +2 "$run_dir/prediction_inputs.tsv")

final_state=complete
if (( failed > 0 )); then
    final_state=completed_with_failures
fi
write_status "$final_state"
finalized=1
printf 'Finished %d distinct sequences with %d failures. Results: %s\n' "$total" "$failed" "$run_dir"
if (( failed > 0 )); then
    exit 1
fi
