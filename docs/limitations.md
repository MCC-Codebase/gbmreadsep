# Limitations and failure modes

## Fundamental limits (bulk DNA-seq)

Bulk DNA reads are not cell-type-labeled. Tumor and TME/immune cells share the same germline genome.
Therefore, read separation is only possible via **tumor-specific somatic events**.

GBMReadSep assigns reads based on whether they overlap **somatic SNV anchors**. Reads that do not overlap
anchors are intrinsically ambiguous and default to the purity prior.

## Practical failure modes

- **Subclonality**: a tumor-derived read can show REF at a somatic locus; absence of ALT is not proof of non-tumor.
- **Copy-number / LOH**: the relationship between mixture VAF and tumor allele fraction is altered.
- **Anchor quality**: tumor-only calling or noisy variants can produce false anchors that mislabel reads.
- **Alignment artifacts**: low MAPQ, paralogs, local assembly issues can bias evidence.

## Recommendations

- Prefer **tumor-normal somatic calling** and strict FILTER=PASS anchors.
- Provide a reliable **purity estimate** when available (pathology/ABSOLUTE/PureCN/FACETS).
- Use split BAMs as *enriched subsets* for downstream tasks (assembly, SV validation), not as truth.

