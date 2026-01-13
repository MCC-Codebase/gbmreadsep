# Interpretation guide

This guide explains how to interpret GBMReadSep outputs, including plots and uncertainty.

---

## Key concepts

**Posterior (`p_tumor`)**
- A probability that a read comes from tumor, given its anchor evidence.
- Values near 1.0 → strong tumor evidence.
- Values near 0.0 → strong non‑tumor evidence.

**Classes**
- **T**: tumor
- **N**: non‑tumor / TME
- **U**: uncertain (not enough evidence)

---

## Why UNASSIGNED is normal

Most reads do **not** overlap anchor SNVs. Those reads naturally have little evidence and fall into **U**.
This is expected in WGS coverage and does **not** indicate failure.

---

## Tumor‑enriched BAMs

If you enable `--split-bams` or `--write-tagged-bam`, the “tumor” BAM is **enriched**, not pure.
There will always be some non‑tumor reads due to:
- limited anchors
- sequencing noise
- ambiguous evidence

---

## Thresholds and uncertainty

**Defaults**
- Tumor: `p_tumor ≥ 0.9`
- Normal: `p_tumor ≤ 0.1`
- Otherwise: Unassigned

If you want stricter tumor calls, raise `--tumor-threshold` (e.g., 0.95).
If you want fewer UNASSIGNED reads, lower thresholds (trade‑off: more mis‑classification risk).

---

## Suggested interpretation workflow

1) Start with the report summary and class counts.
2) Look at posterior histogram: do you see clear peaks near 0 and 1?
3) Check evidence histogram: reads with more anchors give more confident calls.
4) Use downstream tools on the tumor‑enriched BAM, but remember it is enriched, not pure.

---

## Common patterns

- **Very few tumor reads**: consider low tumor purity or low anchor count.
- **Posterior ~0.5 for many reads**: can indicate germline anchors or contig mismatch.
- **Low evidence per read**: indicates sparse anchors; consider lowering filters or using more anchors.

---

## Recommended thresholds (starting points)

- WGS (moderate coverage): keep defaults.
- Ultra‑deep targeted: you can tighten thresholds because evidence is richer.
- Tumor‑only calling: consider stricter tumor thresholds due to germline noise.
