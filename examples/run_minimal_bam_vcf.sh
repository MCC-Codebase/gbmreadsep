#!/usr/bin/env bash
set -euo pipefail

OUTDIR="${1:-toy}"
RESULTS="${2:-toy_out}"

# Generate tiny inputs
gbmreadsep make-toy-data --outdir "${OUTDIR}"

# Run assignment
gbmreadsep assign \
  --bam "${OUTDIR}/tumor.bam" \
  --vcf "${OUTDIR}/anchors.vcf.gz" \
  --outdir "${RESULTS}"

echo "Report: ${RESULTS}/report.html"
