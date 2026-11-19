"""Serve the original research website plus an exact RNA complementarity API."""
import csv
import gzip
import json
import os
import tempfile
import threading
import uuid
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

from flask import Flask, jsonify, request, send_file, send_from_directory

ROOT = Path(__file__).resolve().parents[2]
WEB = Path(__file__).resolve().parent
TARGET = ROOT / 'src/config_code/hpv_step2/targets/HPV16_E6starI_E7_spliced_region.rna.fasta'
HUMAN_CSV = ROOT / 'results/lncrna_checks.csv'
HUMAN_FASTA = ROOT / 'data/gencode.v50.lncRNA_transcripts.fa.gz'
MAX_TRANSCRIPTS = 197195
MAX_PREVIEW = 200
RESULTS = {}
LOCK = threading.Lock()
app = Flask(__name__, static_folder=None)


def normalize(value):
    return str(value).strip().upper().replace('T', 'U')


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGU', 'UGCA'))[::-1]


def viral_target():
    if not TARGET.is_file():
        raise FileNotFoundError(f'HPV16 target FASTA missing: {TARGET}')
    identifier = None
    chunks = []
    for line in TARGET.read_text(encoding='utf-8').splitlines():
        if line.startswith('>'):
            if identifier is not None:
                raise ValueError('Expected a single viral target FASTA record')
            identifier = line[1:].split()[0]
        elif line.strip():
            chunks.append(normalize(line))
    sequence = ''.join(chunks)
    if not identifier or not sequence:
        raise ValueError('Viral target FASTA is empty')
    return identifier, sequence


def human_records():
    if HUMAN_CSV.is_file():
        with HUMAN_CSV.open(newline='', encoding='utf-8-sig') as stream:
            reader = csv.DictReader(stream)
            if not {'seq_id', 'seq'}.issubset(reader.fieldnames or []):
                raise ValueError('Human CSV requires seq_id and seq columns')
            for row in reader:
                yield row['seq_id'], row['seq']
    elif HUMAN_FASTA.is_file():
        with gzip.open(HUMAN_FASTA, 'rt', encoding='utf-8') as stream:
            identifier, chunks = None, []
            for line in stream:
                if line.startswith('>'):
                    if identifier is not None:
                        yield identifier, ''.join(chunks)
                    identifier, chunks = line[1:].split()[0], []
                elif line.strip():
                    chunks.append(line.strip())
            if identifier is not None:
                yield identifier, ''.join(chunks)
    else:
        raise FileNotFoundError('Neither results/lncrna_checks.csv nor the GENCODE lncRNA FASTA was found')


@app.get('/')
def home():
    return send_from_directory(WEB, 'index.html')


@app.get('/api/status')
def status():
    return jsonify({'target_available': TARGET.is_file(), 'human_available': HUMAN_CSV.is_file() or HUMAN_FASTA.is_file(),
                    'human_source': 'CSV' if HUMAN_CSV.is_file() else ('compressed FASTA' if HUMAN_FASTA.is_file() else None),
                    'target': 'HPV16 E6*I/E7 spliced RNA'})


@app.get('/api/target')
def target():
    try:
        identifier, sequence = viral_target()
        return jsonify({'id': identifier, 'length_nt': len(sequence), 'sequence': sequence})
    except (ValueError, OSError) as exc:
        return jsonify({'error': str(exc)}), 404


@app.post('/api/analyze')
def analyze():
    payload = request.get_json(silent=True) or {}
    try:
        k = int(payload.get('match_length', 10))
        limit = int(payload.get('transcript_limit', 100))
    except (ValueError, TypeError):
        return jsonify({'error': 'Match length and transcript limit must be integers'}), 400
    if not 8 <= k <= 20 or limit not in (100, 500, 1000, 10000, 50000, 197195):
        return jsonify({'error': 'Invalid analysis settings'}), 400
    try:
        target_id, viral = viral_target()
        index = defaultdict(list)
        for start in range(len(viral) - k + 1):
            subsequence = viral[start:start+k]
            if set(subsequence) <= set('ACGU'):
                index[reverse_complement(subsequence)].append((start, subsequence))
        if not index:
            return jsonify({'error': 'Target contains no valid windows of this length'}), 400
        columns = ['Human Transcript','Viral Target','Human Start','Human End','Viral Start','Viral End',
                   'Viral Sequence','Reverse Complement','Human Matched Sequence','Match Length','Mismatches','Gaps']
        fd, output_path = tempfile.mkstemp(prefix='hpv_analysis_', suffix='.csv')
        total, analyzed, valid = 0, 0, 0
        transcript_counts, position_counts = Counter(), Counter()
        preview = []
        try:
            with os.fdopen(fd, 'w', newline='', encoding='utf-8') as out:
                writer = csv.DictWriter(out, fieldnames=columns)
                writer.writeheader()
                for transcript_id, raw_sequence in human_records():
                    if analyzed >= limit:
                        break
                    analyzed += 1
                    sequence = normalize(raw_sequence)
                    transcript_id = transcript_id.strip()
                    if not transcript_id or not sequence:
                        continue
                    valid += 1
                    for human_start in range(len(sequence) - k + 1):
                        human_kmer = sequence[human_start:human_start+k]
                        for viral_start, viral_kmer in index.get(human_kmer, ()):
                            row = {'Human Transcript': transcript_id, 'Viral Target': target_id,
                                   'Human Start': human_start, 'Human End': human_start+k,
                                   'Viral Start': viral_start, 'Viral End': viral_start+k,
                                   'Viral Sequence': viral_kmer, 'Reverse Complement': human_kmer,
                                   'Human Matched Sequence': human_kmer, 'Match Length': k, 'Mismatches': 0, 'Gaps': 0}
                            writer.writerow(row)
                            total += 1
                            transcript_counts[transcript_id] += 1
                            position_counts[viral_start] += 1
                            if len(preview) < MAX_PREVIEW:
                                preview.append(row)
        except Exception:
            Path(output_path).unlink(missing_ok=True)
            raise
        run_id = uuid.uuid4().hex
        with LOCK:
            RESULTS[run_id] = output_path
        return jsonify({'run_id': run_id, 'target_id': target_id, 'viral_length': len(viral),
                        'match_length': k, 'analyzed': analyzed, 'valid_transcripts': valid,
                        'total_matches': total, 'unique_transcripts': len(transcript_counts),
                        'preview': preview, 'transcript_counts': transcript_counts.most_common(15),
                        'position_counts': sorted(position_counts.items()),
                        'started_at_utc': datetime.now(timezone.utc).isoformat(),
                        'coordinates': 'zero-based, end-exclusive',
                        'note': 'Exact reverse-complement matches only; not binding predictions.'})
    except (ValueError, OSError, csv.Error) as exc:
        return jsonify({'error': str(exc)}), 400


@app.get('/api/download/<run_id>')
def download(run_id):
    with LOCK:
        path = RESULTS.get(run_id)
    if not path or not Path(path).is_file():
        return jsonify({'error': 'Run not found; rerun the analysis'}), 404
    return send_file(path, mimetype='text/csv', as_attachment=True, download_name='hpv_full_analysis.csv')


@app.get('/<path:filename>')
def assets(filename):
    return send_from_directory(WEB, filename)


if __name__ == '__main__':
    app.run(host='127.0.0.1', port=5050, debug=False)
