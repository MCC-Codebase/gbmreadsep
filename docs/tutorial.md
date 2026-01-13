# Tutorial: read assignment from BAM + somatic VCF

This tutorial assumes you have a **coordinate-sorted, indexed BAM** for a bulk glioblastoma DNA-seq sample
(WGS/WES/panel) and a **somatic SNV VCF** containing tumor-specific variants.

For Oxford Nanopore WGS inputs, see: `docs/nanopore_tutorial.md`.

> Important: bulk DNA-seq cannot perfectly separate reads by cell type. GBMReadSep provides a Bayesian
> *best-effort* tumor-vs-non-tumor assignment using only the subset of reads that overlap somatic anchors.

---

## 0) Install GBMReadSep

From the repository root:

```bash
pip install -e .
```

---

## 1) Prepare inputs

### 1.1 Coordinate-sort and index your BAM

If not already done:

```bash
samtools sort -o tumor.sorted.bam tumor.bam
samtools index tumor.sorted.bam
```

### 1.2 Obtain a somatic SNV VCF

Recommended: tumor-normal calling (e.g., GATK Mutect2 or Strelka2).

If you only have tumor, you can still produce a tumor-only callset, but anchor quality will be worse.
Use strict filtering (PASS only) and prefer known GBM driver mutations if available.

Your VCF should contain, for the tumor sample:

- `AF` or `AD` (allelic depths) and ideally `DP`.

---

## 2) Run read assignment

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

If you do **not** know purity:

```bash
gbmreadsep assign \
  --bam tumor.sorted.bam \
  --vcf somatic.vcf.gz \
  --sample TUMOR \
  --outdir results
```

GBMReadSep will estimate purity heuristically from VAFs.

---

## 3) Inspect outputs

- `results/report.html` – main human-readable report
- `results/assignments.tsv.gz` – per-read assignments (posterior, class, evidence count)
- `results/tagged.bam` – BAM annotated with:
  - `TP` = posterior P(tumor)
  - `TC` = class {T,N,U}
  - `TE` = number of anchor observations used
  - `TL` = summed log-likelihood ratio
- `results/tumor.bam`, `results/tme.bam`, `results/uncertain.bam` – if `--split-bams`

---

## 4) Interpretation guidelines

- A large fraction of reads may have `TE=0` (no anchors overlapped) and thus `TP≈purity`.
- Treat `T` as a **tumor-enriched** subset, not a perfect tumor-only set.
- Treat `N` as **non-tumor-enriched** (TME/normal), not immune-specific.

---

## 5) Common parameter tuning

- Increase stringency (fewer but cleaner reads):
  - `--tumor-threshold 0.95`
  - `--normal-threshold 0.05`
  - `--min-baseq 30`
- Increase yield (more reads, more uncertainty):
  - `--tumor-threshold 0.8`
  - `--normal-threshold 0.2`

---

## 6) Reproducibility

All parameters and counts are recorded in `summary.json` and `anchor_stats.json`.

For Oxford Nanopore WGS (FASTQ/BAM → somatic VCF → assignment), see `docs/nanopore_tutorial.md`.
