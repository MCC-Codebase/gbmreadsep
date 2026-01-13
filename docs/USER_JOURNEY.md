# GBMReadSep user journey

This document describes the **primary workflows**, user expectations, and minimum prerequisites.

---

## Minimum prerequisites

**You need:**
- Python 3.10+ with GBMReadSep installed
- A **tumor BAM** (sorted + indexed)
- A **somatic VCF** (bgzip + tabix indexed recommended)
- Reference FASTA if you plan to align FASTQ or call variants

**Check prerequisites**:

```bash
gbmreadsep doctor
```

---

## Workflow A: BAM + VCF → read assignment (most common)

**When to use**: you already have a tumor BAM and a somatic VCF.

**Command**
```bash
gbmreadsep assign \
  --bam tumor.bam \
  --vcf somatic.vcf.gz \
  --outdir results/
```

**What you get**
- `results/report.html` — summary with plots
- `results/assignments.tsv.gz` — per‑read assignments
- `results/summary.json` — machine‑readable summary

**Key expectation**: Many reads will be **UNASSIGNED** because they do not overlap anchors.

---

## Workflow B: FASTQ end‑to‑end (ONT)

**When to use**: you have FASTQ and want a one‑command pipeline.

```bash
gbmreadsep nanopore endtoend \
  --tumor-fastq tumor.fq.gz \
  --normal-fastq normal.fq.gz \
  --ref ref.fa \
  --outdir ont_run/
```

**What happens internally**
1) Align FASTQ → BAM (minimap2 + samtools)
2) Call somatic VCF (ClairS / ClairS‑TO)
3) Assign reads + render report

---

## Workflow C: “I already have a VCF”

**When to use**: You have a somatic VCF from another caller (Mutect2, Strelka, etc.).

Use the same command as workflow A. GBMReadSep only needs **BAM + VCF**.

---

## Workflow D: “I only have tumor”

**When to use**: No matched normal is available.

```bash
gbmreadsep nanopore endtoend \
  --tumor-fastq tumor.fq.gz \
  --ref ref.fa \
  --outdir ont_tumor_only/
```

**Important caution**: Tumor‑only calling can retain germline variants. This is expected and should be handled with downstream filtering or caution in interpretation.

---

## What users should expect

- **UNASSIGNED reads are normal.** Most reads do not overlap anchor SNVs.
- **Tumor‑enriched BAM** (if enabled) means **enriched**, not pure tumor.
- **Purity matters.** If unknown, GBMReadSep will estimate it, but providing a realistic value improves stability.

---

## Troubleshooting guidance

- Missing BAM index → `samtools index tumor.bam`
- VCF not bgzip/tabix indexed → `bgzip -c file.vcf > file.vcf.gz; tabix -p vcf file.vcf.gz`
- Contig mismatch (chr1 vs 1) → re‑run with `--contig-style ucsc` or `--contig-style ensembl`
- Check `logs/` in the output directory for details

See [docs/troubleshooting.md](troubleshooting.md) for a full catalog.
