# GBMReadSep

**GBMReadSep** assigns bulk DNA sequencing reads to **tumor** vs **non-tumor (TME/normal)** using a
Bayesian evidence model driven by **somatic SNV anchor variants**.

It produces:

- **Per-read posterior** probability `P(tumor | evidence)` and a discrete class (`T`, `N`, `U`).
- Optionally: a **tagged BAM** and/or **split BAMs** (`tumor.bam`, `tme.bam`, `uncertain.bam`).
- QC plots + a self-contained **HTML report** (counts, uncertainty, evidence coverage).

---

## What this tool can and cannot do

### Can do

- Enrich for tumor reads **that overlap tumor-specific somatic anchors** (SNVs) and quantify uncertainty.
- Provide **per-read** posterior probabilities and assignment tags.

### Cannot do (fundamental limitations)

- Perfectly separate all reads by cell type in bulk DNA-seq.
  Reads that do **not** overlap informative somatic loci are intrinsically ambiguous and will default to the
  purity-based prior.
- Separate immune subtypes (T cell vs myeloid, etc.) from DNA alone.

---

## Installation

From the repository root:

```bash
pip install .
# or editable during development
pip install -e ".[dev]"
```

---


## Installation (recommended)

Two supported paths are provided for easy, reproducible setup.

### Option A: Conda/Mamba (recommended for non-computational users)

```bash
# from the repo root
mamba env create -f environment.yml || mamba env update -f environment.yml
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gbmreadsep

# verify
gbmreadsep doctor
pytest -q
```

### Option B: Ubuntu/Debian venv + apt (minimal)

```bash
bash scripts/bootstrap_ubuntu.sh
```

You can also use:

```bash
bash scripts/bootstrap_conda.sh
```

### Makefile shortcuts

```bash
make env
make doctor
make test
```


## Quickstart

You need:

1. A coordinate-sorted, indexed BAM (`.bam` + `.bai`)
2. A somatic VCF (`.vcf`/`.vcf.gz`) containing tumor-specific variants (preferably tumor-normal calling)

Run:

```bash
gbmreadsep assign \
  --bam tumor.sorted.bam \
  --vcf somatic.vcf.gz \
  --sample TUMOR \
  --outdir results \
  --purity 0.6 \
  --write-tagged-bam results/tagged.bam \
  --split-bams
```

Open:

- `results/report.html`

---

## Oxford Nanopore WGS helpers

GBMReadSep core functionality only needs **BAM + somatic VCF**.
For **Oxford Nanopore WGS**, the package also provides optional helper commands to:

- align FASTQ → BAM (minimap2 + samtools)
- call somatic variants from ONT BAM(s) using **ClairS** (tumor/normal) or **ClairS-TO** (tumor-only)
- run an end-to-end pipeline FASTQ/BAM → VCF → read-separation report

### Prerequisites (Nanopore end-to-end)

- `minimap2` and `samtools` in your `PATH` (only needed if you start from FASTQ)
- Docker (recommended) for ClairS/ClairS-TO, or a native install of those tools

Check your setup:

```bash
gbmreadsep doctor
```

### End-to-end: FASTQ → VCF → read separation

Tumor-only:

```bash
gbmreadsep nanopore endtoend \
  --tumor-fastq tumor.fastq.gz \
  --ref GRCh38.fa \
  --outdir results_ont \
  --threads 16
```

Tumor + matched normal:

```bash
gbmreadsep nanopore endtoend \
  --tumor-fastq tumor.fastq.gz \
  --normal-fastq normal.fastq.gz \
  --ref GRCh38.fa \
  --outdir results_ont \
  --threads 16
```

If you already have aligned BAM(s):

```bash
gbmreadsep nanopore endtoend \
  --tumor-bam tumor.sorted.bam \
  --normal-bam normal.sorted.bam \
  --ref GRCh38.fa \
  --outdir results_ont
```

Top-level report:

- `results_ont/report.html` (links to the read-separation report)

### Just get a somatic VCF

Tumor-only (ClairS-TO):

```bash
gbmreadsep nanopore call-vcf \
  --tumor-bam tumor.sorted.bam \
  --ref GRCh38.fa \
  --outdir somatic_vcf \
  --threads 16
```

Tumor/normal (ClairS):

```bash
gbmreadsep nanopore call-vcf \
  --tumor-bam tumor.sorted.bam \
  --normal-bam normal.sorted.bam \
  --ref GRCh38.fa \
  --outdir somatic_vcf \
  --threads 16
```

The command prints the VCF path and also writes `call_vcf_summary.json`.

---

## Outputs

Core (`assign`) output directory contains:

- `report.html` – main human-readable report
- `assignments.tsv.gz` – per-read assignments (`TP`, class, evidence count, etc.)
- `summary.json` – machine-readable summary
- `anchor_stats.json` – anchor selection stats
- `plots/*.png` – QC plots
- `tagged.bam` – if requested
- `tumor.bam`, `tme.bam`, `uncertain.bam` – if `--split-bams`

Nanopore (`nanopore endtoend`) output directory contains:

- `alignment/` – BAMs if you started from FASTQ
- `somatic_vcf/` – somatic VCF(s)
- `readsep/` – read separation outputs (including `readsep/report.html`)
- `report.html` – top-level pipeline report linking everything

---

## Interpretation tips

- Expect a large fraction of reads with `TE=0` (no anchors overlapped). Those reads default to the purity prior.
- Treat `T` as **tumor-enriched**, not perfectly tumor-only.
- Treat `N` as **non-tumor-enriched** (TME/normal), not immune-specific.

---

## Documentation

- `docs/tutorial.md` – step-by-step starting from BAM + somatic VCF
- `docs/nanopore_tutorial.md` – step-by-step ONT WGS workflow (BAM/FASTQ → VCF → assignment)
- `docs/limitations.md` – interpretation caveats

---

## License

MIT.
