#!/usr/bin/env bash
set -euo pipefail

# Example Oxford Nanopore end-to-end run (tumor + matched normal).

TUMOR_FASTQ="tumor.fastq.gz"
NORMAL_FASTQ="normal.fastq.gz"
REF="GRCh38.fa"
OUTDIR="results_ont"
THREADS=16

gbmreadsep nanopore endtoend \
  --tumor-fastq "$TUMOR_FASTQ" \
  --normal-fastq "$NORMAL_FASTQ" \
  --ref "$REF" \
  --outdir "$OUTDIR" \
  --threads "$THREADS" \
  -v

echo "Open: $OUTDIR/report.html"
