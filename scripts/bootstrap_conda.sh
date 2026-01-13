#!/usr/bin/env bash
set -euo pipefail

command -v mamba >/dev/null 2>&1 || {
  echo "mamba not found. Install micromamba/mamba first."; exit 1;
}

mamba env create -f environment.yml || mamba env update -f environment.yml

# Activate (works for conda and mamba)
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gbmreadsep

gbmreadsep doctor
pytest -q
