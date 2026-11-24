# BINF Human lncRNA–HPV16 Research Website

## Overview

This research website presents computational analyses of potential RNA interactions between human long noncoding RNAs (lncRNAs) and the HPV16 E6*I/E7 spliced RNA target.

The website combines the verified Step 5 results for five exploratory transcripts with an interactive RNA Explorer that allows users to perform exact reverse-complementarity screening.

The interface is built using HTML, CSS, and JavaScript, with a Python Flask backend supporting live computational analysis.

## Website Features

- Interactive RNA interaction energy histograms
- Comparison of original sequences against randomized RNA controls
- Sortable candidate comparison tables
- Detailed transcript annotations
- Interactive RNA Explorer
- Adjustable k-mer length and transcript limits
- Exact reverse-complementarity screening
- HPV16 target sequence information
- Conceptual RNA splicing visualization
- Color-coded nucleotide alignments
- Interactive viral position map
- Downloadable computational results and charts

## Installation and Local Setup

### Requirements

- Python 3
- Git
- A modern web browser
- Flask
- Pandas
- Biopython

### Clone the Repository

For a new checkout, clone the repository and enter its root directory.

```bash
git clone https://github.com/chaitu84/Project_JACE.git
cd Project_JACE
```

For an existing checkout, retrieve the latest changes.

```bash
git pull
```

### Windows

Open PowerShell in the repository root directory.

Install the required Python dependencies.

```powershell
py -3 -m pip install flask pandas biopython
```

Start the website server.

```powershell
py -3 src/components_UI/server.py
```

### macOS, Linux, or Ubuntu/WSL

Open a terminal in the repository root directory.

Install the required Python dependencies.

```bash
python3 -m pip install flask pandas biopython
```

Start the website server.

```bash
python3 src/components_UI/server.py
```

### Access the Website

Once the server is running, open the following address in your browser:

http://127.0.0.1:5050/

Keep the terminal open while using the website. Press Ctrl+C to stop the server.

If your Python installation uses the command `python`, substitute it for `py -3` or `python3`.

**Important:** The website must be launched through `server.py` to enable the interactive RNA Explorer. Python's standard HTTP server does not support the backend analysis functionality.

## RNA Explorer

The RNA Explorer provides an interactive computational environment for identifying regions of exact reverse complementarity between human lncRNAs and the HPV16 E6*I/E7 spliced RNA target.

### Analysis Parameters

**Match Length (k-mer)**

Determines the number of consecutive complementary nucleotides required for a reported match.

A match length of 10 nucleotides serves as an exploratory screening parameter. Shorter sequences are more likely to produce random matches, while longer sequences impose more restrictive matching requirements.

This parameter is adjustable and is not intended to represent a biologically validated threshold for gene silencing or therapeutic activity.

**Transcript Limit**

Determines how many human lncRNA records are included in the computational analysis. Larger transcript limits increase the scope of the search but may require additional processing time.

### Computational Methodology

The Python backend performs the following operations:

1. Loads the HPV16 target sequence and available human lncRNA sequences.
2. Normalizes nucleotide sequences into the RNA alphabet.
3. Constructs an index of viral target k-mers.
4. Searches human lncRNA sequences for exact reverse-complementary matches.
5. Records transcript identifiers, complementary sequences, and matching coordinates.
6. Returns computational results to the interactive website.

### Understanding the Results

The RNA Explorer reports:

- **Exact matches:** The total number of complementary sequence windows identified.
- **Unique matching transcripts:** The number of distinct transcripts containing at least one match.
- **Records scanned:** The number of transcript records processed.
- **Match length:** The selected nucleotide window length.

Individual transcripts may contain multiple matching windows. Overlapping windows and repeated matching sites are counted separately.

Coordinates are reported using zero-based, end-exclusive indexing.

### Sequence Alignment

Users can select an individual match to inspect the corresponding human and viral RNA sequences.

The alignment demonstrates canonical RNA base pairing:

- Adenine (A) pairs with uracil (U).
- Guanine (G) pairs with cytosine (C).

Complementary RNA sequences are displayed in opposite orientations to illustrate reverse complementarity.

### Viral Position Map

The Viral Position Map visualizes the distribution of complementary matches across the 573-nucleotide HPV16 target.

- Near-white indicates intervals without identified matches.
- Light teal indicates lower match density.
- Dark teal indicates higher match density.
- Slate blue highlights the currently selected complementary match.

Match density represents the number of complementary windows overlapping a target interval. It does not represent binding strength or thermodynamic stability.

### Data Export

Users can download the complete computational match results as a CSV file for additional analysis and reproducibility.

### Required Research Data

The RNA Explorer requires access to the human lncRNA sequence dataset and HPV16 target FASTA file.

Expected paths include:

```text
results/lncrna_checks.csv

src/config_code/hpv_step2/targets/
    HPV16_E6starI_E7_spliced_region.rna.fasta
```

Large generated research datasets may be excluded from Git tracking. These files must be obtained or generated separately if they are not included in the local checkout.

The prepared Step 5 visualizations use the bundled website JSON and CSV files.

## Explore the Step 5 Results

- Hover over histogram bars, or focus a chart and use the left/right arrow keys, to inspect actual control counts and example energies.
- Home/End select the first/last bin, and Escape closes the tooltip.
- Click a comparison-table heading to sort the results.
- Click a transcript to inspect its annotations, repeat-match context, opening costs, and original predicted coordinates.
- Download the comparison CSV, all 500 control energies, prepared JSON, or SVG charts.

Histograms use the same 13 two-kcal/mol bins from −36 to −10 and a shared 0–25 control-count axis.

Bins include the lower edge and exclude the upper edge, with the final upper edge included. Each plot contains all 100 saved energies for that candidate. Original energies have labeled vertical markers. Hover examples identify actual shuffle numbers.

## Randomized RNA Controls

RNA sequences are shuffled to establish a randomized baseline for evaluating predicted interaction energies.

The shuffle procedure preserves nucleotide counts while rearranging their order. This allows comparisons between original sequences and randomized sequences with the same nucleotide composition.

Each of the five exploratory transcripts was compared against 100 shuffled controls, producing 500 control-energy measurements.

Dinucleotide frequencies and natural RNA structure may change during shuffling.

## Scientific Interpretation

All 500 controls completed with reported interactions. Each control contributes one best reported energy.

The exploratory energy tail fraction is calculated as:

```text
(count(E_shuffle <= E_original) + 1) / 101
```

All five fractions exceed 0.05. These fractions describe the saved control comparison and differ from the Step 4 exact-match p-values and BH-adjusted values.

Step 4 tested all 197,195 transcripts and produced zero BH rejections at FDR 0.05. The minimum adjusted value was approximately 0.679814.

Candidate selection used the same exploratory screening data. Four candidates share a 626-nt block and the same original predicted span sequences, while their full sequences and control-energy distributions differ.

Repeated and overlapping CCU-rich windows account for much of the exact-match count.

Interaction energies include hybridization and accessibility/opening costs.

## Scientific Limitations

Exact reverse complementarity does not independently establish functional RNA binding, gene silencing, or therapeutic activity.

Additional considerations include:

- RNA secondary structure
- Target accessibility
- Thermodynamic stability
- Cellular conditions
- Longer complementary regions
- Statistical significance
- Experimental validation

The RNA Explorer is intended as an initial computational screening tool rather than a validated therapeutic prediction system.

## Website Files

```text
src/components_UI/
├── index.html
├── styles.css
├── app.js
├── rna-explorer.js
├── server.py
├── prepare_data.py
├── verify_site.py
├── qa_report.json
└── data/
    ├── website-data.json
    ├── candidate-comparison.csv
    └── control-energies.csv
```

### Frontend

- `index.html` defines the website structure.
- `styles.css` provides the visual design and responsive layout.
- `app.js` handles the prepared research visualizations.
- `rna-explorer.js` manages interactive RNA analysis and result visualization.

### Backend

- `server.py` launches the Flask server and supports live computational analysis.

### Prepared Research Data

- `website-data.json` contains 500 actual control energies, five candidate summaries, annotations, coordinates, definitions, software versions, common chart axes, and relative scientific source paths with hashes.
- `candidate-comparison.csv` contains the candidate comparison results.
- `control-energies.csv` contains the randomized control-energy results.

No full GENCODE sequences or large Step 4 match tables are copied into the prepared website assets.

## Data Preparation and Verification

`prepare_data.py` reads scientific sources and writes only this directory's prepared data.

It verifies all 500 unique transcript/shuffle pairs against native CSVs, completed receipts, input identities, commands, versions, artifact hashes, and saved summaries.

It independently recalculates means, medians, inclusive tail counts, and fractions. All 4,135 inspected scientific source-file hashes were preserved during the original data preparation.

An older `prediction_summary.csv` was written before control completion and retains a historical `pending` field.

The completed invocation `run_20261116T051139789910Z_0558ca3f` supplies the current control results. This reconciliation appears in the data provenance. There were no numeric discrepancies.

### Full-Source Verification

The optional full-source checker requires the original research results.

On macOS or Linux:

```bash
python3 src/components_UI/verify_site.py
```

On Windows:

```powershell
py -3 src/components_UI/verify_site.py
```

### Quality Assurance

The original `qa_report.json` records independent data, download, local-link, and chart-domain checks.

JavaScript histogram checks verified all 500 values, bin-edge handling, original markers, and common axes. Local HTTP requests returned the expected bytes for the original HTML, CSS, JavaScript, JSON, and CSV downloads.

Text and chart-color contrast were checked numerically.

The original browser automation was blocked, leaving automated desktop/mobile visual inspection, interactive browser execution, and browser-console inspection unverified.

These original verification results do not constitute independent verification of the subsequently added RNA Explorer or Flask backend.

## Project Purpose

This project combines bioinformatics, computational RNA sequence analysis, statistical controls, and interactive data visualization to investigate potential interactions between human lncRNAs and HPV16.

The website is designed to make computational research more accessible, interpretable, and reproducible while maintaining appropriate scientific limitations.
