# FAQ

## What purity should I use?

If you know the tumor purity (from pathology or orthogonal estimates), provide it via `--purity`.
If you don’t, GBMReadSep estimates purity from VAFs and warns you if the estimate is unreliable.

## How many UNASSIGNED reads should I expect?

Often **many**. Most reads will not overlap anchor SNVs. UNASSIGNED is normal and expected.

## Can I use tumor‑only calling?

Yes, but tumor‑only calls can include germline variants. You should interpret tumor‑only results
carefully and apply downstream germline filtering when possible.

## Why are many anchors near VAF ~0.5?

This can indicate germline variants or copy‑number effects. GBMReadSep warns if a large fraction
of anchors are near 0.5 because that can reduce read assignment specificity.

## What is a tumor‑enriched BAM?

Reads labeled “T” are **enriched** for tumor origin, not guaranteed tumor‑only.
Use with caution in downstream analyses.
