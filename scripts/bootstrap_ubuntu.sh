#!/usr/bin/env bash
set -euo pipefail

# System deps (Ubuntu/Debian)
sudo apt-get update
sudo apt-get install -y \
  python3-venv python3-pip \
  samtools minimap2 bcftools tabix pigz

# Python env
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -e ".[dev]"

# Verify
gbmreadsep doctor
pytest -q
