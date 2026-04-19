#!/usr/bin/env bash
set -euo pipefail

# Usage:
# bash scripts/SRmdup_pipeline.sh input.before.bam SAMPLE_NAME [OUTDIR] [THREADS]
#
# Example:
# bash scripts/SRmdup_pipeline.sh examples/example.before.bam EXAMPLE property_out 4

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 INPUT_BAM SAMPLE_NAME [OUTDIR] [THREADS]"
  exit 1
fi

INPUT_BAM="$1"
SAMPLE_NAME="$2"
OUTDIR="${3:-property_out}"
THREADS="${4:-4}"

mkdir -p "$OUTDIR"

# Resolve script directory so this pipeline works from the repository root or from any other working directory.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRMDUP_PY="$SCRIPT_DIR/SRmdup.py"
READS_LEN_PY="$SCRIPT_DIR/Read_lengths.py"
PLOT_R="$SCRIPT_DIR/Properties_plot.R"

# -------- output files --------
AFTER_BAM="$OUTDIR/${SAMPLE_NAME}.after.bam"
STATS_TSV="$OUTDIR/stats.tsv"

READS_BEFORE_TSV="$OUTDIR/read_lengths.before.tsv"
READS_AFTER_TSV="$OUTDIR/read_lengths.after.tsv"

COV_BEFORE_TSV="$OUTDIR/reads_coverage.before.tsv"
COV_AFTER_TSV="$OUTDIR/reads_coverage.after.tsv"

echo "[1/4] Run duplicate filtering"
python3 "$SRMDUP_PY" \
  --bam "$INPUT_BAM" \
  --out-bam "$AFTER_BAM" \
  --stats-tsv "$STATS_TSV" \
  --read-threads "$THREADS" \
  --write-threads "$THREADS"

echo "[2/4] Extract read lengths (before/after)"
python3 "$READS_LEN_PY" "$INPUT_BAM" "$READS_BEFORE_TSV"
python3 "$READS_LEN_PY" "$AFTER_BAM" "$READS_AFTER_TSV"

echo "[3/4] Calculate depth (before/after)"
samtools depth -@ "$THREADS" "$INPUT_BAM" > "$COV_BEFORE_TSV"
samtools depth -@ "$THREADS" "$AFTER_BAM" > "$COV_AFTER_TSV"

echo "[4/4] Draw a plot"
Rscript "$PLOT_R" \
  --sample_label "$SAMPLE_NAME" \
  --qc_dir "$OUTDIR" \
  --out_dir "$OUTDIR" \
  --out_prefix "$SAMPLE_NAME"

echo "All done."
echo "Before BAM : $INPUT_BAM"
echo "After BAM  : $AFTER_BAM"
echo "Stats TSV  : $STATS_TSV"
echo "Output dir : $OUTDIR"