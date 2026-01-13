# Nanopore WGS tutorial (FASTQ/BAM → somatic VCF → read separation)

This tutorial is written for wet-lab users who want a predictable, copy/paste workflow.

GBMReadSep **needs** a tumor BAM and a somatic VCF (SNVs). For Oxford Nanopore WGS, GBMReadSep
also provides helper commands to (optionally) create these inputs.

---

## 0) One-time setup

### Install GBMReadSep

In a clean Python environment:

```bash
pip install gbmreadsep  # once you publish it
# or from the cloned repo
pip install .
```

### (Recommended) Install alignment tools

If you want GBMReadSep to run alignment from FASTQ, you need `minimap2` and `samtools`:

```bash
mamba create -n gbmreadsep -c conda-forge -c bioconda python=3.11 minimap2 samtools
mamba activate gbmreadsep
pip install .
```

### Install Docker (recommended for VCF calling)

Somatic variant calling is executed through Docker by default (ClairS / ClairS-TO). Make sure:

```bash
docker --version
```

Then run:

```bash
gbmreadsep doctor
```

---

## 1) Prepare inputs

### Reference FASTA

You need a reference FASTA (e.g. GRCh38) **matching the reference used for alignment and variant calling**.

Your reference folder should contain:

- `GRCh38.fa`
- `GRCh38.fa.fai` (FASTA index)

GBMReadSep will create `.fai` automatically if missing.

### Tumor and normal sequencing data

Supported starting points:

- **FASTQ/FASTQ.GZ** from Nanopore basecalling
- **sorted + indexed BAM** from Dorado/minimap2 alignment

---

## 2) Align FASTQ → BAM (optional)

If you already have sorted/indexed BAM files, skip this step.

Tumor alignment:

```bash
gbmreadsep nanopore align \
  --fastq tumor.fastq.gz \
  --ref GRCh38.fa \
  --out-bam tumor.sorted.bam \
  --threads 16 \
  --preset lr:hq
```

Normal alignment (optional):

```bash
gbmreadsep nanopore align \
  --fastq normal.fastq.gz \
  --ref GRCh38.fa \
  --out-bam normal.sorted.bam \
  --threads 16 \
  --preset lr:hq
```

Notes:

- `lr:hq` is recommended for modern Q20+ ONT reads.
- For noisier reads, try `--preset map-ont`.

---

## 3) Call a somatic VCF (BAM → VCF)

### Option A: Tumor + matched normal (recommended)

```bash
gbmreadsep nanopore call-vcf \
  --tumor-bam tumor.sorted.bam \
  --normal-bam normal.sorted.bam \
  --ref GRCh38.fa \
  --outdir somatic_vcf \
  --threads 16 \
  --platform ont_r10_dorado_sup_5khz
```

This uses **ClairS** (paired calling). The main output is:

- `somatic_vcf/output.vcf.gz`

### Option B: Tumor-only

```bash
gbmreadsep nanopore call-vcf \
  --tumor-bam tumor.sorted.bam \
  --ref GRCh38.fa \
  --outdir somatic_vcf \
  --threads 16 \
  --platform ont_r10_dorado_sup_5khz
```

This uses **ClairS-TO**. Outputs:

- `somatic_vcf/snv.vcf.gz`
- `somatic_vcf/indel.vcf.gz`

Use the SNV VCF for read separation.

---

## 4) Assign reads: tumor vs TME (BAM + VCF → report)

```bash
gbmreadsep assign \
  --bam tumor.sorted.bam \
  --vcf somatic_vcf/output.vcf.gz \
  --outdir readsep \
  --purity 0.6 \
  --split-bams \
  --write-tagged-bam readsep/tagged.bam
```

Open:

- `readsep/report.html`

---

## 5) One-command end-to-end pipeline

Tumor-only:

```bash
gbmreadsep nanopore endtoend \
  --tumor-fastq tumor.fastq.gz \
  --ref GRCh38.fa \
  --outdir results_ont \
  --threads 16
```

Tumor + normal:

```bash
gbmreadsep nanopore endtoend \
  --tumor-fastq tumor.fastq.gz \
  --normal-fastq normal.fastq.gz \
  --ref GRCh38.fa \
  --outdir results_ont \
  --threads 16
```

Then open:

- `results_ont/report.html`

---

## Interpretation guardrails

- **Expect uncertainty.** Reads not overlapping anchors are not classifiable from DNA alone.
- **Class T = tumor-enriched**, not guaranteed tumor-only.
- **Class N = non-tumor-enriched**, not immune-specific.

If you have a good purity estimate (histology / flow / copy-number based), set `--purity`.

