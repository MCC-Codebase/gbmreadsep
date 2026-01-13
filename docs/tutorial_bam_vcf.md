# Tutorial: BAM + VCF read assignment

This tutorial is designed for biologists who have a tumor BAM and a somatic VCF.

---

## Step 1: Check prerequisites

```bash
gbmreadsep doctor
```

You should see OK for Python, samtools, and any optional tools you plan to use.

---

## Step 2: Run assignment

```bash
gbmreadsep assign \
  --bam tumor.bam \
  --vcf somatic.vcf.gz \
  --outdir results/
```

**Expected console output (example)**
```
Assigning reads: 100%|██████████| 12345/12345 [00:02<00:00, 4567.8read/s]
results/report.html
```

---

## Step 3: Open the report

```bash
xdg-open results/report.html  # Linux
open results/report.html      # macOS
```

---

## Step 4: Interpret the results

Key sections in the report:
- **Class counts**: number of Tumor / TME / Unassigned reads
- **Posterior histogram**: how confident the assignments are
- **Evidence histogram**: how many anchors each read overlapped

See [docs/interpretation.md](interpretation.md) for deeper guidance.

---

## Troubleshooting tips

- **BAM not indexed** → run `samtools index tumor.bam`
- **VCF not indexed** → run `bgzip -c somatic.vcf > somatic.vcf.gz; tabix -p vcf somatic.vcf.gz`
- **Contig mismatch** → rerun with `--contig-style ucsc` or `--contig-style ensembl`

Full list: [docs/troubleshooting.md](troubleshooting.md)
