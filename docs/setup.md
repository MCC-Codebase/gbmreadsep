# Environment setup (Codex-friendly)

This repository includes reproducible environment definitions and bootstrap scripts.

## Conda/Mamba (recommended)

```bash
mamba env create -f environment.yml || mamba env update -f environment.yml
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gbmreadsep

gbmreadsep doctor
pytest -q
```

## Ubuntu/Debian venv + apt

```bash
bash scripts/bootstrap_ubuntu.sh
```

## CI

A GitHub Actions workflow is included in `.github/workflows/ci.yml` that provisions the environment
via `environment.yml` and runs:

- `gbmreadsep doctor`
- `pytest -q`

## Definition of done

A correct setup satisfies:

1. `gbmreadsep --help` works
2. `gbmreadsep doctor` shows OK for python, minimap2, samtools, bcftools, tabix (and docker if installed)
3. `pytest -q` passes
