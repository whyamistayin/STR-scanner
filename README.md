# Polymorphic STR detection pipeline

This repository contains a pipeline for preprocessing genome assemblies, extracting chromosome-specific sequences,
running pangenome analysis using Cactus, Tantan predictions and final analysis using output of previous steps.

## Overview

The workflow consists of three main steps:

1. **Download assemblies** (via `download.sh`)
2. **Preprocess assemblies** (via `preprocessing.sh`)
3. **Run analysis** (via `run.sh`)

---

## Requirements

Make sure the following tools are installed and available in your environment:

- `bash`
- `python3`
- `seqtk`
- `awk`, `grep`, `zcat`
- `docker` (required for running Cactus)
- `Tantan`
- paternal and maternal hg002 reference (named `hg002.pat.fasta` and `hg002.mat.fasta` with chromosomes named "chr{i}")
- Python and Bash scripts:
  - `extract_names.py`
  - `processing.py`
  - `create_namefile.sh`
  - `preprocessing.sh` (!)
  - `run.sh` (!)
  - `download.sh` (!)

You will also need access to the Docker image:
quay.io/comparative-genomics-toolkit/cactus:latest

---

## Step 1: Download Assemblies

A script called `download.sh` is expected to download genome assemblies into the following structure:

./raw_data/*/*/*/*/*_genomic.fna.gz

Run:

bash download.sh

---

## Step 2: Preprocessing

The `preprocessing.sh` script:

- Extracts contig names corresponding to chromosomes
- Filters assemblies to keep only chromosome sequences
- Renames chromosome headers (e.g., `chr1`, `chr2`, ..., `chr22`)
- Splits each filtered assembly into per-chromosome files
- Removes empty outputs

### Run preprocessing:

bash preprocessing.sh

### Output:

- Filtered assemblies:
  *_genomic.fna.filtered

- Per-chromosome files:
  *.chr1, *.chr2, ..., *.chr22

---

## Step 3: Run Analysis

The `run.sh` script performs chromosome-wise pangenome analysis using Cactus.

### Usage:

bash run.sh 

Pangenome analysis refers to maternal hg002 by default.

### What it does:

For each chromosome (1–22):

1. Generates input path files using `create_namefile.sh`
2. Runs `cactus-pangenome` via Docker
3. Runs Tantan with maternal hg002 reference.
4. Processes the resulting graph (`.gfa.gz`) and prediction (`.bed`) using `processing.py`

### Outputs:

- Cactus output directories:
  out_chr1/, out_chr2/, ..., out_chr22/

- Processed logs:
  log.chr1.txt, log.chr2.txt, ..., log.chr22.txt

---

## Notes

- The pipeline assumes diploid reference assemblies (`hg002.mat.fasta` and `hg002.pat.fasta`) are present and already split by chromosome.
- Chromosomes processed: **1–22 only** (autosomes).
- Ensure sufficient disk space and memory, as Cactus can be resource-intensive.
- Docker must be ready before executing `run.sh`.

---

## Pipeline Summary

download.sh  →  preprocessing.sh  →  run.sh
     ↓               ↓                 ↓
 raw_data      filtered + split     analysis results

---

## Troubleshooting

- **Missing chromosome files**: Ensure preprocessing completed successfully.
- **Docker errors**: Verify Docker is installed and running.
- **Empty outputs**: Assemblies may not contain all chromosomes; missing chromosomes are automatically skipped.

---
