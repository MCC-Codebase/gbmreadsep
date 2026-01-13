# CLI reference

This page documents GBMReadSep commands, expected outputs, and common failure modes.

---

## `gbmreadsep quickstart`

**Purpose**: Print three ready‑to‑run recipes and their expected outputs.

**Example**
```bash
gbmreadsep quickstart
```

**Outputs**: Printed to stdout only.

**Errors**: None (no file access).

**Source Anchor**: src/gbmreadsep/cli.py:L492-L521#53c1f533

---

## `gbmreadsep wizard`

**Purpose**: Guided mode that prompts for inputs and prints the final command.

**Key options**
- `--mode {assign,endtoend}`
- `--non-interactive` (require flags; no prompts)
- `--run` (execute the command)
- `--dry-run` (print commands only)

**Example**
```bash
gbmreadsep wizard --mode assign --bam tumor.bam --vcf somatic.vcf.gz \
  --outdir results --purity 0.7 --non-interactive
```

**Outputs**: Printed command. Optionally runs the command.

**Errors**: Missing required flags when `--non-interactive` is set.

**Source Anchor**: src/gbmreadsep/cli.py:L538-L632#268fafe6

---

## `gbmreadsep make-toy-data`

**Purpose**: Generate a tiny reference, BAM, and VCF for demos/tests.

**Example**
```bash
gbmreadsep make-toy-data --outdir toy/
```

**Outputs**
- `toy/toy_ref.fa`
- `toy/tumor.bam`
- `toy/anchors.vcf.gz`
- `toy/toy_summary.json`

**Errors**: Write failures if the outdir is not writable.

**Source Anchor**: src/gbmreadsep/cli.py:L635-L643#a765a060

---

## `gbmreadsep assign`

**Purpose**: Assign reads to tumor vs non‑tumor using a BAM + somatic VCF.

**Example**
```bash
gbmreadsep assign --bam tumor.bam --vcf somatic.vcf.gz --outdir results/
```

**Key outputs**
- `results/report.html`
- `results/assignments.tsv.gz`
- `results/summary.json`
- `results/logs/assign.log`

**Common errors**
- Missing BAM index → `samtools index tumor.bam`
- Missing VCF index → `bgzip -c file.vcf > file.vcf.gz; tabix -p vcf file.vcf.gz`
- Contig mismatch → `--contig-style ucsc` or `--contig-style ensembl`

**Source Anchor**: src/gbmreadsep/cli.py:L646-L758#528c2458

---

## `gbmreadsep nanopore align`

**Purpose**: Align ONT FASTQ to a reference and create sorted + indexed BAM.

**Example**
```bash
gbmreadsep nanopore align --fastq reads.fq.gz --ref ref.fa --out-bam tumor.bam
```

**Outputs**
- `tumor.bam` and `tumor.bam.bai`
- `logs/alignment.log`

**Errors**
- Missing minimap2 or samtools → run `gbmreadsep doctor` for install commands

**Source Anchor**: src/gbmreadsep/cli.py:L796-L834#1bc8ada6

---

## `gbmreadsep nanopore call-vcf`

**Purpose**: Call somatic variants from ONT BAM(s).

**Example**
```bash
gbmreadsep nanopore call-vcf --tumor-bam tumor.bam --ref ref.fa --outdir somatic_vcf/
```

**Outputs**
- `somatic_vcf/snv.vcf.gz` (tumor‑only)
- `somatic_vcf/output.vcf.gz` (tumor‑normal)
- `somatic_vcf/logs/caller.log`

**Errors**
- Missing docker or caller binaries
- Missing BAM index

**Source Anchor**: src/gbmreadsep/cli.py:L837-L883#14b1c4c9

---

## `gbmreadsep nanopore endtoend`

**Purpose**: Full pipeline from FASTQ/BAM to report.

**Example**
```bash
gbmreadsep nanopore endtoend --tumor-fastq tumor.fq.gz --ref ref.fa --outdir ont_run/
```

**Outputs**
- `ont_run/report.html`
- `ont_run/readsep/report.html`
- `ont_run/logs/pipeline.log`

**Errors**: See logs under `ont_run/logs/`.

**Source Anchor**: src/gbmreadsep/cli.py:L886-L972#8803159d
