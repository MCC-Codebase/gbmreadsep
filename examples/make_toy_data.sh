#!/usr/bin/env bash
set -euo pipefail

OUTDIR="${1:-toy}"

gbmreadsep make-toy-data --outdir "${OUTDIR}"

echo "Toy data written to: ${OUTDIR}"
