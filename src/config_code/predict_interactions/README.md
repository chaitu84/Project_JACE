# Step 5: one worked example

**ENST00000832824.1**, GENCODE v50 **DDX11L16-260** (gene ENSG00000290825.2),
is the requested demonstration candidate. It was found in the real-match CSV;
it has not been ranked as the strongest candidate.

## Existing files and their connections

`src/main.py` calls `src/config_code/read.py` for sequence quality checks.
Its `check_fasta()` reads plain/gzipped FASTA and takes the versioned transcript
ID before the first `|` in a GENCODE header. The lncRNA path is:

```text
data/gencode.v50.lncRNA_transcripts.fa.gz
    -> read.py -> results/lncrna_checks.csv
    -> complementary.py -> results/exact_complementarity_matches.csv
    -> shuffled.py      -> shuffled control CSVs
```

Both searches use the original viral RNA at
`src/config_code/hpv_step2/targets/HPV16_E6starI_E7_spliced_region.rna.fasta`.
Step 2's `scripts/build_targets.py` builds targets from its `reference/` files;
`scripts/verify_targets.py` provides an independent audit. This spliced target
is a **573-nt reference-derived partial region**, not an established complete
viral transcript.

Also present: the all-transcript GENCODE FASTA and `all_transcript_checks.csv`,
E6/E7 quality-check CSVs, DNA and RNA versions of spliced/intron-retained targets,
and `results/shuffled_seq/{shuffled_matches.csv,shuffle_summary.csv}`.
The all-transcript CSV is not the complementarity search input.
`Predict.py` is an unfinished stub; the new code does not overwrite it.

`lncrna_checks.csv` has **197,195** unique transcript rows. The real-match CSV
has **164,014** match rows across **77,645** transcript IDs. Both search scripts
currently set `TEST_LIMIT=None`, and the real CSV reaches the last input
transcript. An independent full replay checked every FASTA record against the
CSV and reproduced all 164,014 real-match rows with **zero missing or extra
rows**, confirming full-dataset output rather than a test subset. Evidence is
in `results/step5/worked_example/coverage_audit.json`.

## Validated candidate

The complete GENCODE header is:
`ENST00000832824.1|ENSG00000290825.2|-|-|DDX11L16-260|DDX11L16|1379|`.
The sequence contains **1,379 nt**. All three original rows were collected:

| Row | Human [start,end) | Viral [start,end) | Human 5′→3′ | Viral 5′→3′ |
|---|---|---|---|---|
| 1 | [642,652) | [6,16) | GUGGGUCCUG | CAGGACCCAC |
| 2 | [52,62) | [95,105) | UUGCUUGCAG | CUGCAAGCAA |
| 3 | [1188,1198) | [232,242) | AUCUGCAACA | UGUUGCAGAU |

These are **zero-based, end-exclusive transcript/target coordinates**, not
genomic coordinates. Human [642,652) means one-based inclusive bases 643–652.
Each human slice equals the viral slice's RNA reverse complement. A fresh
exact search of this full candidate reproduces exactly these three rows.

Prepared under `results/step5/worked_example/`:

- `ENST00000832824.1.full.rna.fasta`: full human sequence, original 5′→3′.
- `HPV16_E6starI_E7_spliced_region.original.rna.fasta`: byte-identical copy
  of the original viral RNA FASTA, preserving its header and sequence.
- `candidate_real_matches.csv`: all original candidate rows/columns.
- `candidate_validation.json`: source paths, header, lengths, sequence hashes,
  validation and actual prediction status.
- `intarna_command.txt`: portable command to run in this output folder.
- `coverage_audit.json`: independent coverage/control audit.

The human FASTA changes its **321 T characters to U** to express RNA. Base order,
length and coordinates do not change. Neither input is reverse-complemented.
No essential input is missing. Large source files and original results were
not rewritten. `results/` remains Git-ignored; `.gitignore` now covers future
GENCODE `.fa` and `.fa.gz` inputs too. The lncRNA gz file was **already tracked**;
ignore rules do not untrack it or remove history. No Git index edit, commit or
push was performed.

## IntaRNA and the next action

Official sources checked 2026-11-13:
[manual](https://backofenlab.github.io/IntaRNA/),
[browser interface](https://backofenlab.github.io/IntaRNA/gui/),
[server help](https://rna.informatik.uni-freiburg.de/Help.jsp?toolName=IntaRNA).

Use the human RNA as **target / sequence 1**, virus as **query / sequence 2**.
These are software input roles; the virus remains the project's biological
target. Both inputs stay 5′→3′. The browser permits 2–5,000 nt each and a length
product ≤1,000,000. Here 1,379 × 573 = **790,167**, so full sequences fit.
The Freiburg server also accepts both (4–2,000 nt each). No cropping is needed.
This is suitable for exploratory prediction; validation of IntaRNA in bacterial
systems does not establish accuracy for this human/viral pair.

No local `IntaRNA` executable was found. WSL could not enumerate distributions
in this session, and Docker's engine was unavailable. Preparation succeeded;
prediction, energy, structure and predicted-match overlap remain **pending**.

**Next action without installation:** open the browser interface linked above.
Select **Prediction**. Upload the viral FASTA into **Query RNA** and human FASTA
into **Target RNA**. Keep start indices **1**. Request **1** interaction, use
heuristic prediction, accessibility enabled for both RNAs, folding window
**150**, maximum intramolecular base-pair span **100**, and seed length **7**.
Leave seed positions unrestricted; do not force an original 10-mer.
Choose output columns:
`id1,start1,end1,id2,start2,end2,subseq1,subseq2,hybridDP,bpList,E,E_hybrid,ED1,ED2`.
Click **Compute interactions**, then download the session ZIP to retain
parameters, engine versions and actual outputs. Save its native semicolon CSV
as `results/step5/worked_example/intarna_prediction.csv`.

For standalone execution, IntaRNA must first be installed. The official conda
installation is Linux-only (for example in working WSL with conda):

```bash
conda install -c conda-forge -c bioconda intarna
```

Run this from the worked-example folder containing the two prepared FASTAs:

```bash
IntaRNA --target=ENST00000832824.1.full.rna.fasta --query=HPV16_E6starI_E7_spliced_region.original.rna.fasta --model=X --mode=H --acc=C --tAccW=150 --tAccL=100 --qAccW=150 --qAccL=100 --seedBP=7 --outNumber=1 --outMode=C --outCsvCols=id1,start1,end1,id2,start2,end2,subseq1,subseq2,hybridDP,bpList,E,E_hybrid,ED1,ED2 --out=intarna_prediction.csv
```

If IntaRNA becomes available in Windows, this command from the repository root
prepares the example, runs it, and saves actual version/stdout/stderr/CSV:

```powershell
& "$env:LOCALAPPDATA/Programs/Python/Python313/python.exe" src/config_code/predict_interactions/prepare_step5.py --run-intarna
```

Add `--intarna 'C:/path/to/IntaRNA.exe'` if needed. It refuses to overwrite an
existing prediction CSV. The installed Python was found at this location but
is not on PATH. Preparation reuses `rev_comp`, `read_fasta`, `find_matches`, and
`WINDOW_SIZE` from `complementary.py`.

After saving a real prediction, validate it and check match overlap:

```powershell
& "$env:LOCALAPPDATA/Programs/Python/Python313/python.exe" src/config_code/predict_interactions/prepare_step5.py --prediction results/step5/worked_example/intarna_prediction.csv
```

This writes `prediction_overlap.json`, checking IDs, bounds, subsequences,
overlap on **both** RNAs, and how many original antiparallel base pairs are
present in `bpList`. Interval overlap alone does not establish identical
pairing. No overlap is a legitimate possible prediction outcome.

## Reading output and folding

Native CSV uses `;`. `start1/end1` refer to human, `start2/end2` to viral RNA.
Positions are **one-based inclusive**: convert `[start,end]` to Python
`[start-1:end]`. With full inputs there is no crop offset.

`hybridDP` uses `.` for unpaired bases, `(` on the human side paired to `)` on
the viral side, and `&` between strands. Both subsequences retain original
order; pairing is antiparallel. It shows intermolecular pairing, not the whole
internal fold of each RNA. `bpList` lists explicit coordinate pairs.

`E = E_hybrid + ED1 + ED2`, in kcal/mol. Hybridization can lower energy;
ED1/ED2 penalize opening folded human/viral sites. More negative E is more
favorable within the same model/settings, not a probability, p-value, measured
affinity or evidence of HPV inhibition.

The command searches full sequences but estimates accessibility using **local**
150-nt folding windows and a 100-nt intramolecular base-pair span, following the
documented server settings. It omits long-range folds. Heuristic mode can miss
the global best interaction. An exact 10-mer does not by itself show that a
site will be accessible in a cell.

## Shuffled-control review

| Setting | Real search | Shuffled search |
|---|---|---|
| Window size | 10 | 10 |
| Viral step | 1 (564 windows) | 10 (57 windows) |
| Transcript CSV | lncrna_checks.csv | Same |
| Current test limit | None | None |
| Human matches | All overlapping occurrences | Same imported function |
| Repeated viral windows | Count at each position | Same |
| Match count | Human/viral position pairs | Same per shuffle |

Raw counts are unfair because the real run searches almost ten times as many
viral starting positions. A future fair control needs `WINDOW_STEP=1` in
`shuffled.py` and fresh output files. This example does not launch that rerun
or change existing controls.

Other conditional differences: real skips non-ACGU viral windows, shuffled
does not; real applies the test limit before missing-row removal, shuffled
filters missing/empty rows first. With the present unambiguous target, valid
inputs and no limit, these do not change coverage. Align them before future
subset/ambiguous-input comparisons.

Existing files are in `results/shuffled_seq/`; current code writes
`results/step4/`, which did not exist at inspection. The summary records **18
completed** shuffles, each with 197,195 transcripts, 57 windows, step 10 and
seed 42. Shuffle **19** also has match rows but no completed summary. Current
`NUM_SHUFFLES=100` is not evidence of 100 completed runs. Outputs predate the
current scripts, so current settings alone are not historical run provenance.

Shuffling viral bases preserves length and single-base composition, not
dinucleotides or folding. It supplies a **chance baseline**, not the real
candidate shortlist. One empty shuffled run would not prove significance.
The existing runs are nonempty, and a fair comparison still needs matched
search rules and completed replicates.
