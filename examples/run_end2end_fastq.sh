#!/usr/bin/env bash
set -euo pipefail

TUMOR_FASTQ=${1:-tumor.fq.gz}
NORMAL_FASTQ=${2:-normal.fq.gz}
REF=${3:-ref.fa}
OUTDIR=${4:-ont_run}

gbmreadsep nanopore endtoend \
  --tumor-fastq "${TUMOR_FASTQ}" \
  --normal-fastq "${NORMAL_FASTQ}" \
  --ref "${REF}" \
  --outdir "${OUTDIR}"

echo "Report: ${OUTDIR}/report.html"
