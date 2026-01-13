# Tutorial: FASTQ end‑to‑end (ONT)

This tutorial covers the full pipeline: FASTQ → BAM → somatic VCF → read assignment.

---

## Step 1: Prepare inputs

You need:
- Tumor FASTQ files (ONT)
- (Optional) normal FASTQ files
- A reference FASTA (same build you will use for calling)

---

## Step 2: Run the end‑to‑end pipeline

```bash
gbmreadsep nanopore endtoend \
  --tumor-fastq tumor.fq.gz \
  --normal-fastq normal.fq.gz \
  --ref ref.fa \
  --outdir ont_run/
```

**Expected console output (example)**
```
ont_run/report.html
```

---

## Step 3: Open the reports

```bash
xdg-open ont_run/report.html
xdg-open ont_run/readsep/report.html
```

---

## Step 4: Check logs if needed

Logs live in:
- `ont_run/logs/pipeline.log`
- `ont_run/logs/assign.log`
- `ont_run/logs/caller.log`

---

## Tumor‑only note

If you do not have a matched normal, you can omit `--normal-fastq`.
Tumor‑only calling can retain germline variants, so interpret results accordingly.

---

## Troubleshooting tips

- Missing tools → `gbmreadsep doctor`
- Stuck or failed calling → check `ont_run/logs/caller.log`
- Reference mismatch → ensure alignment and calling use the same reference

More details: [docs/troubleshooting.md](troubleshooting.md)
