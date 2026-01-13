SHELL := /bin/bash

.PHONY: env install doctor test

env:
	@command -v mamba >/dev/null 2>&1 || (echo "mamba not found. Install micromamba/mamba first." && exit 1)
	mamba env create -f environment.yml || mamba env update -f environment.yml

install:
	pip install -e .[dev]

doctor:
	gbmreadsep doctor

test:
	pytest -q
