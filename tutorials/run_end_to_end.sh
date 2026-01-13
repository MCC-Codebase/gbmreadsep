#!/usr/bin/env bash
set -euo pipefail

# Example end-to-end run. Adjust paths for your environment.

BAM="tumor.sorted.bam"
VCF="somatic.vcf.gz"
SAMPLE="TUMOR"
OUTDIR="results"

gbmreadsep assign \
  --bam "${BAM}" \
  --vcf "${VCF}" \
  --sample "${SAMPLE}" \
  --outdir "${OUTDIR}" \
  --purity 0.6 \
  --write-tagged-bam "${OUTDIR}/tagged.bam" \
  --split-bams \
  -v
