# SRmdup
# SRmdup
## Overview
SRmdup is a read-level duplicate filtering method for Single-end low-coverage ancient human DNA.

This repository contains scripts for:
1. read-level duplicate filtering of BAM files,
2. extraction of read-level tables before and after filtering,
3. drawing a figure visualizing the filtering statistics.

SRmdup performs read-level duplicate filtering using corrected outer boundaries.  
The current implementation removes duplicates in four categories: same start and same end, same start, same end, and start and end within 1 bp.  

## Environment
SRmdup requires the following environment:
- Python 3.9+
- R 4.2+
- samtools
- bash
Python package:
- `pysam`
R packages:
- `data.table`
- `ggplot2`
- `scales`
- `patchwork`

## Installation

Clone the repository and create a conda environment:
```bash
mkdir SRmup
git clone <your-repo-url>
cd SRmdup
conda create -n srmdup python=3.10 -y
conda activate srmdup
conda install -c conda-forge -c bioconda pysam samtools r-base r-data.table r-ggplot2 r-scales r-patchwork

## Scripts
- `SRmdup.py`: main script for read-level duplicate filtering.
- `SRmdup_pipeline.sh`: wrapper script for running the full workflow.
- `Reads_lengths.py`: extracts read-level tables before and after filtering.
- `Properties_plot.R`: draws a summary figure from intermediate files.

## Usage
Run the full pipeline:
```bash
bash scripts/SRmdup_pipeline.sh input.before.bam SAMPLE output_dir 4
Run duplicate filtering only:
```bash
python3 scripts/SRmdup.py --bam input.before.bam --out-bam output.after.bam --stats-tsv stats.tsv
More example commands are available in examples/example_commands.txt.

## Input and output
### Input
- a BAM file before duplicate filtering. 
### Output
- filtered BAM file 
- a table with the count of deleted reads by strand and category.
- two .tsv files with all reads and their length and cigar string before and after filtering 
- two .tsv files with the depth of all covered positions before and after filtering 
- a plot showing strands, categories of the deleted reads, read length distribution and coverage depth before and after filtering.
### Optional outputs from SRmdup.py
- deleted-read tables 
- kept-read table 
## Example
Additional example commands for SRmdup
### 1. Duplicate filtering with optional deleted-read and kept-read tables
```bash
python3 scripts/SRmdup.py \
  --bam examples/example.bam \
  --out-bam examples/example.after.bam \
  --stats-tsv examples/stats.tsv \
  --del-prefix examples/deleted_reads \
  --keep-tsv examples/kept_reads.tsv

### 2. Duplicate filtering with custom temporary directory and threads
```bash
python3 scripts/SRmdup.py \
  --bam examples/example.bam \
  --out-bam examples/example.after.bam \
  --stats-tsv examples/stats.tsv \
  --tmp-dir examples/tmp \
  --read-threads 4 \
  --write-threads 4 \
  --progress-every 50000

### 3. Allow filtering on an unsorted BAM file
```bash
python3 scripts/SRmdup.py \
  --bam examples/example.bam \
  --out-bam examples/example.after.unsorted.bam \
  --stats-tsv examples/stats.unsorted.tsv \
  --accept-unsorted
