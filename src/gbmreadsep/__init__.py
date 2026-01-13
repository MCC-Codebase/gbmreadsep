"""GBMReadSep: Bayesian read-level tumor vs TME assignment for bulk DNA sequencing.

Public API is intentionally small; most users should use the CLI:

    gbmreadsep assign --bam ... --vcf ... --outdir ...

"""

from __future__ import annotations

__all__ = ["__version__"]

__version__ = "0.2.0"
