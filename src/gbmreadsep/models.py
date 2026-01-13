from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class AnchorVariant:
    """A somatic SNV used as an anchor for read assignment.

    Coordinates are 0-based half-open in internal representation.

    Attributes
    ----------
    chrom:
        Contig name as present in BAM/VCF.
    pos0:
        0-based genomic position of the SNV (reference coordinate).
    ref:
        Reference base (A/C/G/T), uppercase.
    alt:
        Alternate base (A/C/G/T), uppercase.
    vaf:
        Variant allele fraction observed in the mixture sample (if available).
    dp:
        Depth at the locus (if available).
    f_tumor:
        Estimated probability that a *tumor-derived* read carries ALT at this locus.
        This is NOT the mixture VAF; it is a tumor-allele fraction proxy.
    record_id:
        Optional identifier (VCF ID or CHROM:POS:REF:ALT).
    """

    chrom: str
    pos0: int
    ref: str
    alt: str
    vaf: Optional[float]
    dp: Optional[int]
    f_tumor: float
    record_id: str


@dataclass(frozen=True)
class ReadAssignment:
    """Per-read posterior assignment result."""

    qname: str
    chrom: str
    start0: int
    end0: int
    is_read1: bool
    is_read2: bool
    mapq: int
    evidence_count: int
    llr_sum: float
    p_tumor: float
    cls: str  # 'T', 'N', or 'U'
