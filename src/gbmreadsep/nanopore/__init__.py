"""Nanopore-specific convenience wrappers.

This subpackage is *optional*.

GBMReadSep core functionality only needs a BAM + somatic VCF.
For Oxford Nanopore WGS/WES, this subpackage provides helper commands to:

- align FASTQ -> sorted/indexed BAM (minimap2 + samtools)
- call somatic variants (ClairS / ClairS-TO; typically via Docker)
- run an end-to-end pipeline FASTQ/BAM -> VCF -> read separation report

The code here aims to be pragmatic and robust, not a full workflow engine.
"""

from __future__ import annotations

__all__ = [
    "align_ont_fastq_to_bam",
    "call_somatic_vcf",
    "nanopore_endtoend",
]

from .align import align_ont_fastq_to_bam
from .somatic import call_somatic_vcf
from .pipeline import nanopore_endtoend
