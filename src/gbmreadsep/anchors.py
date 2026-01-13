from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple

import pysam

from .models import AnchorVariant
from .utils import clamp

logger = logging.getLogger(__name__)

_MIN_ANCHORS_WARNING = 50
_HIGH_VAF_WINDOW = (0.45, 0.55)
_HIGH_VAF_FRACTION_WARN = 0.3


@dataclass(frozen=True)
class AnchorIndex:
    """Per-contig anchor lookup structure."""

    positions: List[int]  # sorted 0-based positions
    anchors: List[AnchorVariant]  # aligned with positions


def _extract_vaf_and_dp(sample: pysam.libcbcf.VariantRecordSample) -> Tuple[Optional[float], Optional[int]]:
    """Best-effort extraction of (VAF, DP) from a VCF sample field."""
    # Common fields across callers:
    # - AF (float)
    # - AD (ref, alt) and DP
    # - FREQ (% in some callers)
    vaf: Optional[float] = None
    dp: Optional[int] = None

    if "DP" in sample and sample["DP"] is not None:
        try:
            dp = int(sample["DP"])
        except Exception:
            dp = None

    if "AF" in sample and sample["AF"] is not None:
        af = sample["AF"]
        try:
            # AF can be list-like or scalar
            if isinstance(af, (list, tuple)):
                vaf = float(af[0])
            else:
                vaf = float(af)
        except Exception:
            vaf = None

    if vaf is None and "FREQ" in sample and sample["FREQ"] is not None:
        # Sometimes a string like "12.3%"
        try:
            s = sample["FREQ"]
            if isinstance(s, (list, tuple)):
                s = s[0]
            if isinstance(s, str) and s.endswith("%"):
                vaf = float(s[:-1]) / 100.0
            else:
                vaf = float(s)
        except Exception:
            vaf = None

    if vaf is None and "AD" in sample and sample["AD"] is not None:
        try:
            ad = sample["AD"]
            if isinstance(ad, (list, tuple)) and len(ad) >= 2:
                refc = float(ad[0])
                altc = float(ad[1])
                denom = refc + altc
                if denom > 0:
                    vaf = altc / denom
                    if dp is None:
                        dp = int(round(denom))
        except Exception:
            vaf = None

    return vaf, dp


def estimate_purity_from_vafs(
    vafs: Iterable[float],
    *,
    min_vaf: float = 0.05,
    max_vaf: float = 0.45,
) -> Optional[float]:
    """Rough purity estimate from mixture VAFs assuming many clonal het SNVs.

    We use the heuristic purity ≈ median(2 * VAF) over a VAF window.

    This is intentionally conservative and may be biased in high-CN/LOH contexts.
    """
    vals: List[float] = []
    for v in vafs:
        if min_vaf <= v <= max_vaf:
            vals.append(clamp(2.0 * v, 0.0, 1.0))
    if len(vals) < 25:
        return None
    vals.sort()
    mid = len(vals) // 2
    if len(vals) % 2 == 1:
        return vals[mid]
    return 0.5 * (vals[mid - 1] + vals[mid])


def load_anchor_variants(
    vcf_path: str,
    *,
    sample: Optional[str] = None,
    purity: Optional[float] = None,
    min_dp: int = 10,
    min_vaf: Optional[float] = None,
    max_vaf: Optional[float] = None,
    min_qual: Optional[float] = None,
    require_pass: bool = True,
    max_anchors: int = 200_000,
    default_f_tumor: float = 0.5,
) -> Tuple[List[AnchorVariant], float, Dict[str, int]]:
    """Load somatic SNV anchors from VCF.

    Parameters
    ----------
    vcf_path:
        Somatic VCF (optionally bgzip+tabix indexed).
    sample:
        Sample name within VCF. If None, uses first sample.
    purity:
        Tumor purity prior. If None, estimated from VAFs (heuristic).
    min_dp:
        Minimum depth (DP) for anchor inclusion (if DP available).
    min_vaf, max_vaf:
        Optional VAF window for anchors. If set, anchors without VAF are skipped.
        This is useful when you want to focus on high-confidence clonal-ish anchors
        for read separation.
    min_qual:
        Optional minimum VCF QUAL for anchors (if QUAL is present). Records without
        QUAL are kept.
    require_pass:
        If True, require FILTER to be PASS or empty.
    max_anchors:
        Upper bound on the number of anchors (for runtime control). If exceeded, keeps most
        informative anchors based on VAF and DP heuristics.
    default_f_tumor:
        Used when VAF is missing; should be in [0,1].

    Returns
    -------
    anchors:
        List of AnchorVariant objects (sorted by chrom, pos).
    purity_used:
        Purity prior used downstream.
    stats:
        Simple counters about records kept/skipped.
    """
    vcf = pysam.VariantFile(vcf_path)

    if sample is None:
        if len(vcf.header.samples) == 0:
            raise ValueError(
                "VCF has no samples. Provide a VCF with tumor sample fields containing VAF/AD/DP."
            )
        sample = list(vcf.header.samples)[0]
        logger.info("No --sample provided; using first VCF sample: %s", sample)
    if sample not in vcf.header.samples:
        raise ValueError(f"Sample '{sample}' not found in VCF samples: {list(vcf.header.samples)}")

    stats: Dict[str, int] = {
        "records_total": 0,
        "records_pass": 0,
        "records_snv": 0,
        "records_with_vaf": 0,
        "anchors_kept": 0,
        "anchors_skipped_filter": 0,
        "anchors_skipped_non_snv": 0,
        "anchors_skipped_dp": 0,
        "anchors_skipped_vaf": 0,
        "anchors_skipped_qual": 0,
    }

    raw: List[Tuple[str, int, str, str, Optional[float], Optional[int], str]] = []
    vafs_for_purity: List[float] = []

    # Iteration strategy: VCF.fetch() is fast but requires an index for bgzipped VCFs.
    # For portability, fall back to sequential iteration if no index is present.
    try:
        iterator = vcf.fetch()
    except (ValueError, OSError):
        iterator = vcf

    for rec in iterator:
        stats["records_total"] += 1

        if require_pass:
            # In pysam, rec.filter.keys() returns set of filter names; PASS may be absent when empty.
            filt = list(rec.filter.keys())
            if len(filt) > 0 and not (len(filt) == 1 and filt[0] == "PASS"):
                stats["anchors_skipped_filter"] += 1
                continue
        stats["records_pass"] += 1

        ref = rec.ref
        alts = list(rec.alts or [])
        if len(ref) != 1 or len(alts) != 1 or len(alts[0]) != 1:
            stats["anchors_skipped_non_snv"] += 1
            continue
        alt = alts[0]
        if ref.upper() not in {"A", "C", "G", "T"} or alt.upper() not in {"A", "C", "G", "T"}:
            stats["anchors_skipped_non_snv"] += 1
            continue
        stats["records_snv"] += 1

        # Optional QUAL filter
        if min_qual is not None and rec.qual is not None and float(rec.qual) < float(min_qual):
            stats["anchors_skipped_qual"] += 1
            continue

        s = rec.samples[sample]
        vaf, dp = _extract_vaf_and_dp(s)
        if dp is not None and dp < min_dp:
            stats["anchors_skipped_dp"] += 1
            continue

        # Optional VAF window filter
        if min_vaf is not None or max_vaf is not None:
            if vaf is None:
                stats["anchors_skipped_vaf"] += 1
                continue
            lo = float(min_vaf) if min_vaf is not None else 0.0
            hi = float(max_vaf) if max_vaf is not None else 1.0
            if not (lo <= float(vaf) <= hi):
                stats["anchors_skipped_vaf"] += 1
                continue

        chrom = str(rec.contig)
        pos0 = int(rec.pos) - 1  # VCF is 1-based
        rid = rec.id if rec.id is not None else f"{chrom}:{rec.pos}:{ref}:{alt}"

        raw.append((chrom, pos0, ref.upper(), alt.upper(), vaf, dp, rid))
        if vaf is not None:
            stats["records_with_vaf"] += 1
            vafs_for_purity.append(float(vaf))

    if len(raw) == 0:
        raise ValueError("No usable SNV anchors were found in the VCF after filtering.")

    purity_used = purity
    if purity_used is None:
        est = estimate_purity_from_vafs(vafs_for_purity)
        if est is None:
            purity_used = 0.5
            logger.warning(
                "Purity was not provided and could not be estimated reliably from VAFs; "
                "defaulting to purity=0.5. Consider providing --purity."
            )
        else:
            purity_used = clamp(est, 0.05, 0.95)
            logger.info("Estimated purity from VAFs: %.3f", purity_used)
    else:
        purity_used = clamp(float(purity_used), 0.01, 0.99)

    default_f_tumor = clamp(float(default_f_tumor), 0.0, 1.0)

    anchors: List[AnchorVariant] = []
    for chrom, pos0, ref, alt, vaf, dp, rid in raw:
        if vaf is None:
            f_tumor = default_f_tumor
        else:
            # Under a diploid mixture heuristic: VAF ≈ purity * f_tumor
            f_tumor = clamp(float(vaf) / purity_used, 0.0, 1.0)

        anchors.append(
            AnchorVariant(
                chrom=chrom,
                pos0=pos0,
                ref=ref,
                alt=alt,
                vaf=vaf,
                dp=dp,
                f_tumor=f_tumor,
                record_id=rid,
            )
        )

    # Downselect if too many anchors
    if len(anchors) > max_anchors:
        logger.warning(
            "Anchor count (%d) exceeds max_anchors=%d; downselecting for runtime control.",
            len(anchors),
            max_anchors,
        )
        # Heuristic: prefer higher VAF and higher DP.
        def score(a: AnchorVariant) -> float:
            v = a.vaf if a.vaf is not None else 0.0
            d = float(a.dp if a.dp is not None else 0)
            return v * (1.0 + 0.01 * d)

        anchors.sort(key=score, reverse=True)
        anchors = anchors[:max_anchors]

    # Sort by chrom then pos; keep original contig ordering from BAM later.
    anchors.sort(key=lambda a: (a.chrom, a.pos0))
    stats["anchors_kept"] = len(anchors)

    if len(anchors) < _MIN_ANCHORS_WARNING:
        logger.warning(
            "Only %d anchors found. Read-level assignments may be noisy. "
            "Consider lowering filters or providing more anchors.",
            len(anchors),
        )

    vafs_present = [a.vaf for a in anchors if a.vaf is not None]
    if vafs_present:
        lo, hi = _HIGH_VAF_WINDOW
        frac = sum(1 for v in vafs_present if lo <= float(v) <= hi) / float(len(vafs_present))
        if frac >= _HIGH_VAF_FRACTION_WARN:
            logger.warning(
                "Many anchors have VAF near 0.5 (%.0f%%). This can indicate germline variants; "
                "consider tumor-normal calling or germline filtering.",
                frac * 100.0,
            )
    return anchors, purity_used, stats


def build_anchor_index(anchors: List[AnchorVariant]) -> Dict[str, AnchorIndex]:
    """Build a per-contig index for fast read overlap lookup."""
    by_contig: Dict[str, List[AnchorVariant]] = {}
    for a in anchors:
        by_contig.setdefault(a.chrom, []).append(a)

    index: Dict[str, AnchorIndex] = {}
    for chrom, lst in by_contig.items():
        lst_sorted = sorted(lst, key=lambda x: x.pos0)
        positions = [a.pos0 for a in lst_sorted]
        index[chrom] = AnchorIndex(positions=positions, anchors=lst_sorted)
    return index
