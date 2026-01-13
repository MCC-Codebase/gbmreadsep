# Troubleshooting catalog

This catalog lists common errors and exact fixes.

---

## BAM not indexed

**Message**
```
BAM is not indexed. Run: samtools index tumor.bam
```

**Fix**
```bash
samtools index tumor.bam
```

---

## VCF not bgzip/tabix indexed

**Message**
```
VCF is not bgzip/tabix indexed. Run: bgzip -c file.vcf > file.vcf.gz; tabix -p vcf file.vcf.gz
```

**Fix**
```bash
bgzip -c file.vcf > file.vcf.gz
tabix -p vcf file.vcf.gz
```

---

## Contig mismatch (chr1 vs 1)

**Message**
```
Contig mismatch between BAM and VCF (e.g., chr1 vs 1). Use --contig-style {ucsc,ensembl,auto} to override.
```

**Fix**
```bash
gbmreadsep assign ... --contig-style ucsc
# or
gbmreadsep assign ... --contig-style ensembl
```

---

## Missing external tools

Run:
```bash
gbmreadsep doctor
```

The output includes exact installation commands for Ubuntu and conda/mamba.

---

## Check logs

Every command writes logs under:
```
<outdir>/logs/
```

Common log files:
- `assign.log`
- `alignment.log`
- `caller.log`
- `pipeline.log`

These logs include stack traces and external tool stderr for debugging.
