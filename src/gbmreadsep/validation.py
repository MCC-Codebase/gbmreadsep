from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterable, List

logger = logging.getLogger(__name__)


_UCSC_PREFIX = "chr"


def check_bam_index(bam_path: str | Path) -> None:
    """Ensure a BAM has an index; raise ValueError with fix instructions."""
    bam = Path(bam_path)
    bai1 = bam.with_suffix(bam.suffix + ".bai")
    bai2 = bam.with_suffix(".bai")
    if bai1.exists() or bai2.exists():
        return
    raise ValueError(
        "BAM is not indexed. Run: samtools index " + str(bam)
    )


def check_vcf_index(vcf_path: str | Path) -> None:
    """Ensure a bgzipped VCF has a tabix index; raise ValueError with fix instructions."""
    vcf = Path(vcf_path)
    if vcf.suffixes[-2:] == [".vcf", ".gz"]:
        tbi = vcf.with_suffix(vcf.suffix + ".tbi")
        if not tbi.exists():
            raise ValueError(
                "VCF is not bgzip/tabix indexed. Run: bgzip -c "
                + str(vcf.with_suffix(""))
                + " > "
                + str(vcf)
                + "; tabix -p vcf "
                + str(vcf)
            )
    elif vcf.suffix == ".vcf":
        logger.info(
            "VCF is uncompressed (.vcf). This is supported but slower; "
            "consider bgzip+tabix for large files."
        )


def detect_contig_style(contigs: Iterable[str]) -> str:
    """Infer contig style: 'ucsc' if most contigs start with 'chr', else 'ensembl'."""
    names = [c for c in contigs if c]
    if not names:
        return "unknown"
    chr_like = [c for c in names if c.startswith(_UCSC_PREFIX)]
    if len(chr_like) >= max(1, int(0.5 * len(names))):
        return "ucsc"
    return "ensembl"


def remap_contig(contig: str, style: str) -> str:
    """Remap a contig name to the requested style (ucsc or ensembl)."""
    if style == "ucsc":
        if contig.startswith(_UCSC_PREFIX):
            return contig
        if contig == "MT":
            return "chrM"
        return f"{_UCSC_PREFIX}{contig}"
    if style == "ensembl":
        if contig.startswith(_UCSC_PREFIX):
            core = contig[len(_UCSC_PREFIX) :]
            if core == "M":
                return "MT"
            return core
        return contig
    return contig


def remap_contigs(contigs: Iterable[str], style: str) -> List[str]:
    """Remap a list of contigs to a requested style."""
    return [remap_contig(c, style) for c in contigs]
