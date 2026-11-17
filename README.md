# Project_JACE
BY J.A.C.E (Joint Attachment for Carcinoma Eradication)

## Run the project website

The website presents the completed exploratory human-lncRNA–HPV16 analysis, including five candidates and all 500 saved control energies. Its prepared data is bundled in this repository.

Install Python 3, then clone the repository and start the local server.

```text
git clone https://github.com/chaitu84/Project_JACE.git
cd Project_JACE
```

Windows PowerShell, from the repository root.

```powershell
py -3 -m http.server 8765 --bind 127.0.0.1 --directory src/components_UI
```

On macOS or Linux, use `python3` in place of `py -3`. Open **http://127.0.0.1:8765/**. Full instructions are in [the website README](src/components_UI/README.md).

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

| Name | Major | Role |
|---|---|---|
| Jacky | Computer Science | Coder |
| Arnav | Data Science | Researcher |
| Cleo. | Computer | Coder |
| Ekijah| Mathematics & Economics | Researcher |
