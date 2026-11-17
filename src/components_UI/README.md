# BINF human-lncRNA–HPV16 project website

This self-contained static website presents the verified Step 5 results for five exploratory transcripts. HTML, CSS, JavaScript, prepared data, and verification utilities are contained in this directory. A browser loads the compact local data file, with no framework installation or external asset requests.

## View locally

Install Python 3 and use a modern browser. The prepared data is included with the website, and the server uses Python's standard library.

For a new checkout, clone the repository and enter its root directory.

```text
git clone https://github.com/chaitu84/Project_JACE.git
cd Project_JACE
```

For an existing checkout, open a terminal in its root directory and run `git pull` to get the latest version.

On Windows, run this command from the repository root in PowerShell.

```powershell
py -3 -m http.server 8765 --bind 127.0.0.1 --directory src/components_UI
```

Open **http://127.0.0.1:8765/** in your browser. Keep the server terminal open while viewing the site. Press Ctrl+C to stop it.

On macOS, Linux, or Ubuntu/WSL, run this command from the repository root.

```bash
python3 -m http.server 8765 --bind 127.0.0.1 --directory src/components_UI
```

If your Python 3 installation uses the command `python`, substitute `python` for `py -3` or `python3` above. If port 8765 is already in use, choose another port such as 8766 and open the matching local URL.

If your terminal is already inside `src/components_UI/`, use `py -3 -m http.server 8765 --bind 127.0.0.1` on Windows or `python3 -m http.server 8765 --bind 127.0.0.1` on macOS/Linux.

## Explore the results

- Hover over histogram bars, or focus a chart and use the left/right arrow keys, to inspect actual control counts and example energies. Home/End select the first/last bin, and Escape closes the tooltip.
- Click a comparison-table heading to sort the results. Click a transcript to open its details.
- Choose a transcript to inspect its annotations, repeat-match context, opening costs, and original predicted coordinates.
- Download the comparison CSV, all 500 control energies, the prepared JSON, or SVG charts. The coordinate tracks depict positions along the full RNAs.

Histograms use the same 13 two-kcal/mol bins from −36 to −10 and a shared 0–25 control-count axis. Bins include the lower edge and exclude the upper edge, with the final upper edge included. Each plot contains all 100 saved energies for that candidate. Original energies have labeled vertical markers. Hover examples identify actual shuffle numbers.

## Scientific interpretation

All 500 controls completed with reported interactions. Each control contributes one best reported energy. The exploratory energy tail fraction is `(count(E_shuffle <= E_original) + 1) / 101`. All five fractions exceed 0.05. These fractions describe the saved control comparison and differ from the Step 4 exact-match p-values and BH-adjusted values.

Step 4 tested all 197,195 transcripts and produced zero BH rejections at FDR 0.05. The minimum adjusted value was approximately 0.679814. Candidate selection used the same exploratory screening data. Four candidates share a 626-nt block and the same original predicted span sequences, while their full sequences and control-energy distributions differ. Repeated and overlapping CCU-rich windows account for much of the exact-match count.

The shuffle procedure preserves nucleotide counts. Dinucleotide frequencies and natural structure can change. Interaction energies include hybridization and accessibility/opening costs. Biological binding and therapeutic activity require experimental evidence.

## Files and verification

`data/website-data.json` contains 500 actual control energies, five candidate summaries, annotations, coordinates, definitions, software versions, common chart axes, and relative scientific source paths with hashes. `data/candidate-comparison.csv` and `data/control-energies.csv` are the downloadable tables. No full GENCODE sequences or large Step 4 match tables are copied into the website assets.

`prepare_data.py` reads scientific sources and writes only this directory's prepared data. It verifies all 500 unique transcript/shuffle pairs against native CSVs, completed receipts, input identities, commands, versions, artifact hashes, and saved summaries. It independently recalculates means, medians, inclusive tail counts, and fractions. All 4,135 inspected scientific source-file hashes were preserved.

An older `prediction_summary.csv` was written before control completion and retains a historical `pending` field. The completed invocation `run_20261116T051139789910Z_0558ca3f` supplies the current control results. This reconciliation appears in the data provenance. There were no numeric discrepancies.

The optional full-source checker below is for the research workspace containing the original `results/` files. Those large generated results are Git-ignored. Website viewing uses the bundled JSON and CSV files, and its startup steps are complete above.

With the original research results available, run the checker from the project root.

```bash
python3 src/components_UI/verify_site.py
```

On Windows, use `py -3 src/components_UI/verify_site.py`.

`qa_report.json` records independent data, download, local-link, and chart-domain checks. JavaScript histogram checks verified all 500 values, bin-edge handling, original markers, and common axes. Local HTTP requests returned the expected bytes for HTML, CSS, JavaScript, JSON, and both CSV downloads. Text and chart-color contrast were checked numerically.

Browser automation was blocked when Computer Use could not verify the current Chrome URL. Desktop/mobile visual inspection, interactive browser execution, and browser-console inspection remain unverified. Responsive layouts, keyboard handlers, focus indicators, reduced-motion rules, and local download links are implemented.

Scientific scripts and original results were preserved. The website was prepared for local viewing.
