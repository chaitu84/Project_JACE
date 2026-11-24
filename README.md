# Project_JACE

BY J.A.C.E (Joint Attachment for Carcinoma Eradication)

## Run the project website

The website presents the completed exploratory human-lncRNA–HPV16 analysis, including five candidates, all 500 saved control energies, and an interactive RNA Explorer for live sequence complementarity analysis.

The website uses HTML, CSS, and JavaScript for its interface and a Python Flask backend for computational analysis.

### Installation

Install Python 3, then clone the repository and enter its root directory.

```bash
git clone https://github.com/chaitu84/Project_JACE.git
cd Project_JACE
```

For an existing checkout, open a terminal in the repository root and run:

```bash
git pull
```

### Windows PowerShell

Install the required dependencies.

```powershell
py -3 -m pip install flask pandas biopython
```

Start the server from the repository root.

```powershell
py -3 src/components_UI/server.py
```

### macOS or Linux

Install the required dependencies.

```bash
python3 -m pip install flask pandas biopython
```

Start the server from the repository root.

```bash
python3 src/components_UI/server.py
```

### Access the Website

Open [**http://127.0.0.1:5050/**](http://127.0.0.1:5050/) in your browser.

Keep the server terminal open while viewing the website. Press Ctrl+C to stop it.

The website must be launched through `server.py` to support the interactive RNA Explorer. Python's standard HTTP server does not support live computational analysis.

### Website Features

- Interactive RNA interaction-energy histograms
- Original versus shuffled RNA comparisons
- Sortable candidate comparison tables
- Detailed transcript annotations
- Interactive RNA Explorer
- Adjustable k-mer lengths and transcript limits
- Exact reverse-complementarity screening
- HPV16 target sequence information
- Conceptual RNA splicing visualization
- Color-coded nucleotide alignments
- Interactive viral position mapping
- Downloadable computational results and charts

### Required Research Data

The prepared research visualizations are bundled in this repository. Live RNA Explorer analysis requires the HPV16 target sequence:

```text
src/config_code/hpv_step2/targets/HPV16_E6starI_E7_spliced_region.rna.fasta
```

It also requires one of the following human lncRNA datasets:

```text
results/lncrna_checks.csv
```

or

```text
data/gencode.v50.lncRNA_transcripts.fa.gz
```

Large research datasets may need to be obtained or generated separately if they are not included in the local checkout.

Full instructions are available in the [website README](https://github.com/chaitu84/Project_JACE/blob/main/src/components_UI/README.md).

## Overview

This is a computational biology project looking at whether human non-coding RNA sequences could point toward future RNA-based ways of targeting high-risk HPV.

Right now we're focused on **HPV16** and **HPV18**, specifically the **E6** and **E7** oncogenes, which interfere with human tumor-suppressor pathways.

The project is still in progress. We haven't locked in the final approach, tools, or implementation details yet.

## Problem

High-risk HPV types express E6 and E7, and those genes play a role in how cervical and other HPV-driven cancers develop. The question we want to look at is whether comparing public human non-coding RNA data against HPV reference sequences turns up any sequence relationships worth following up on.

## Goal

Build a **Python tool** that compares public HPV reference sequences to human non-coding RNA datasets and summarizes the matches and predicted RNA interactions that come out of it.

## What the Final Product Could Look Like

A tool that can let someone explore and compare possible RNA-based strategies for targeting HPV oncogene transcripts.

Features we're considering:

- Comparing HPV16 and HPV18 E6/E7 reference sequences
- Pulling in human non-coding RNA datasets from public databases
- Finding and summarizing sequence matches that look relevant
- Running RNA structure or interaction predictions, depending on which tools we go with
- Putting the results into a readable report or a simple interactive visualization

## Data Sources

- **NCBI Nucleotide** or **PaVE** — HPV reference sequences
- **GENCODE** or **Ensembl** — human non-coding RNA annotations

## Tools Required

Python is the main language we're planning on. Beyond that:

- **Biopython** — handling sequences
- **pandas** — organizing and analyzing the data
- **BLAST+** (or something similar) — sequence comparison
- **ViennaRNA** — RNA structure and interaction predictions
- Some kind of visualization tool or lightweight dashboard if it helps show the results

## Scope and Limitations

A sequence match or a predicted RNA interaction doesn't mean a molecule will actually suppress HPV gene expression, or that it'd be safe, or that it'd work as a treatment.

The project only looks at **computational evidence**. Experimental validation, clinical effectiveness, and anything involving delivery through an engineered viral vector are all outside what this software can show.

## Current Status

**HPV16 exploratory analysis completed.** The website presents the verified exact-match screen and interaction-energy controls. All candidates remain exploratory, and none passed Step 4 BH FDR 0.05.

## Team

| **Name** | **Role** |
|---|---|
| Jacky | Coder |
| Arnav | Researcher |
| Cleo. | Coder |
| Ekijah | Researcher |
