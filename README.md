# GBMReadSep

GBMReadSep assigns each sequencing read to **tumor**, **TME (non-tumor)**, or **uncertain** using somatic SNV anchors. It is designed for biologists who want clear, repeatable steps from **install → run → interpret → troubleshoot**.

If you are new here, start with:

```
gbmreadsep quickstart
```

That prints three copy‑paste recipes and tells you exactly which output files to expect.

---

## 5‑minute quickstart (BAM + VCF)

**Inputs you need**
- A tumor BAM (sorted + indexed).
- A somatic VCF for the same tumor (bgzip + tabix indexed recommended).

**Commands**

```bash
# 1) (optional) check tools
gbmreadsep doctor

# 2) run read assignment
gbmreadsep assign \
   --bam tumor.bam \
   --vcf somatic.vcf.gz \
   --outdir results/
```

**Outputs (in `results/`)**
- `report.html` — open in a browser; a plain‑English summary and plots.
- `assignments.tsv.gz` — per‑read posterior and class calls.
- `summary.json` — machine‑readable run summary.
- `plots/` — class counts and posterior histograms.

---

## One‑command end‑to‑end (FASTQ → report)

If you have ONT FASTQ files and want everything done for you:

```bash
gbmreadsep nanopore endtoend \
   --tumor-fastq tumor.fq.gz \
   --normal-fastq normal.fq.gz \
   --ref ref.fa \
   --outdir ont_run/
```

This runs:
1) alignment → BAM
2) somatic calling → VCF
3) read assignment → report

Outputs:
- `ont_run/report.html` — top‑level pipeline report.
- `ont_run/readsep/report.html` — read assignment report.

---

## Output folder tour (what each file means)

**For `gbmreadsep assign`:**
- `report.html`: the main report you should share with collaborators.
- `assignments.tsv.gz`: each read with posterior (`p_tumor`), evidence count, and class.
- `summary.json`: totals and parameters used.
- `plots/`: PNGs used in the report.
- `anchor_stats.json`: counts of anchors used and filters applied.
- `logs/assign.log`: detailed logs (for troubleshooting).

**For `gbmreadsep nanopore endtoend`:**
- `alignment/`: aligned BAM(s).
- `somatic_vcf/`: called VCFs.
- `readsep/`: read separation outputs (same layout as `assign`).
- `pipeline_summary.json`: end‑to‑end summary.
- `logs/pipeline.log`: detailed logs.

---

## Safety features (recommended)

- `--dry-run`: validate inputs and print commands without running tools.
- `--resume`: skip steps that already produced outputs.
- Logs live under `<outdir>/logs/` for every command.

---

## Interpretation (in plain language)

- **Posterior** (`p_tumor`): the probability a read came from tumor, given the anchors it overlaps.
- **Uncertainty is normal.** Many reads do **not** overlap anchors, so they are expected to be **UNASSIGNED**.
- **Tumor‑enriched BAM**: if you enable `--split-bams`, reads labeled “T” are enriched for tumor‑origin reads, not 100% pure tumor.

For a deeper guide, see [docs/interpretation.md](docs/interpretation.md).

---

## Troubleshooting (most common fixes)

- **BAM not indexed** → run:
  ```bash
  samtools index tumor.bam
  ```

- **VCF not bgzip/tabix indexed** → run:
  ```bash
  bgzip -c somatic.vcf > somatic.vcf.gz
  tabix -p vcf somatic.vcf.gz
  ```

- **Contig mismatch (chr1 vs 1)** → use:
  ```bash
  gbmreadsep assign ... --contig-style ucsc
  # or: --contig-style ensembl
  ```

- **Missing tools**: run `gbmreadsep doctor` for exact install commands.

More details: [docs/troubleshooting.md](docs/troubleshooting.md).

---

## Toy demo (runs in < 1 minute)

```bash
# Generate tiny test inputs
gbmreadsep make-toy-data --outdir toy/

# Run read assignment
gbmreadsep assign --bam toy/tumor.bam --vcf toy/anchors.vcf.gz --outdir toy_out/

# Open the report
open toy_out/report.html  # macOS
xdg-open toy_out/report.html  # Linux
```

This is the fastest way to sanity‑check installation and understand outputs.

---

## Guided mode (interactive)

If you want step‑by‑step prompts:

```bash
gbmreadsep wizard
```

You’ll be asked about file paths, purity, and tumor‑only vs tumor‑normal. The wizard prints the final command and can run it for you.

---

## Documentation map

- [docs/USER_JOURNEY.md](docs/USER_JOURNEY.md) — end‑to‑end workflows and expectations
- [docs/tutorial_bam_vcf.md](docs/tutorial_bam_vcf.md) — BAM+VCF step‑by‑step
- [docs/tutorial_fastq_end2end.md](docs/tutorial_fastq_end2end.md) — FASTQ end‑to‑end
- [docs/interpretation.md](docs/interpretation.md) — how to read posteriors and plots
- [docs/troubleshooting.md](docs/troubleshooting.md) — error catalog with fixes
- [docs/faq.md](docs/faq.md) — purity, tumor‑only pitfalls, expected rates
- [docs/glossary.md](docs/glossary.md) — plain‑language definitions

---

## Minimal conceptual background

GBMReadSep uses **somatic SNVs** as anchors. If a read overlaps an anchor, it provides evidence for tumor or non‑tumor origin. Many reads will never touch an anchor, especially at typical WGS coverage, so a large **UNASSIGNED** fraction is normal and expected.

---

## How the module works internally (quick tour)

At a high level, GBMReadSep loads somatic SNV anchors from a VCF, builds per‑contig lookup tables, and then evaluates each read in the BAM by aggregating log‑likelihood ratios across any overlapping anchors.

**Core data structures**
- `AnchorVariant`: a single somatic SNV (chrom, 0‑based position, ref/alt, VAF/DP, and tumor‑allele proxy `f_tumor`).
- `ReadAssignment`: per‑read posterior, evidence count, and class (`T`, `N`, `U`).

**Anchor loading + purity estimation**
- `load_anchor_variants` parses the VCF, keeps PASS SNVs, applies optional DP/VAF/QUAL filters, and collects per‑record stats.
- If purity is not provided, `estimate_purity_from_vafs` uses a conservative median heuristic over a VAF window.

**Read assignment model**
- For each read, candidate anchors are selected by genomic overlap.
- The model computes a per‑anchor log‑likelihood ratio using a simple substitution error model that combines base quality and mapping quality.
- The posterior is `sigmoid(logit(prior_p_tumor) + sum(llr))`, and thresholds classify the read as tumor (T), non‑tumor (N), or uncertain (U).

**Outputs**
- `assignments.tsv.gz` with per‑read posteriors and class calls.
- `summary.json` with run parameters and counts.
- `plots/` and `report.html` for a shareable summary.

---

## Development

```bash
# create env
mamba env create -f environment.yml

# dev install
pip install -e .[dev]

# tests
pytest -q
```
