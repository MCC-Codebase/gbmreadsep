from __future__ import annotations

import logging
import time
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import pysam

from ..external import ensure_executable_in_path, run_command, run_pipe
from ..utils import ensure_outdir

logger = logging.getLogger(__name__)


_SUPPORTED_PRESETS = {
    # Designed for accurate long reads (e.g., Q20+ ONT, Dorado).
    "lr:hq",
    "lr:hqae",
    # Traditional ONT noisy long reads.
    "map-ont",
}


def _ensure_faidx(ref_fa: Path) -> None:
    """Ensure reference FASTA has a .fai index (required by many callers)."""
    fai = ref_fa.with_suffix(ref_fa.suffix + ".fai")
    if fai.exists():
        return
    logger.info("Creating FASTA index: %s", fai)
    # pysam.faidx wraps samtools faidx via htslib; no external samtools required.
    pysam.faidx(str(ref_fa))


def build_alignment_commands(
    *,
    fastq: Sequence[str | Path],
    ref_fa: str | Path,
    out_bam: str | Path,
    threads: int = 8,
    minimap2_preset: str = "lr:hq",
    sort_mem: str = "2G",
    extra_minimap2_args: Optional[List[str]] = None,
    extra_samtools_sort_args: Optional[List[str]] = None,
) -> Dict[str, List[str]]:
    """Build minimap2/samtools commands for ONT alignment."""
    ref_fa = Path(ref_fa).expanduser().resolve()
    out_bam = Path(out_bam).expanduser().resolve()
    fastq_paths = [Path(x).expanduser().resolve() for x in fastq]

    if not ref_fa.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {ref_fa}")
    if len(fastq_paths) == 0:
        raise ValueError("At least one FASTQ file is required")
    for fq in fastq_paths:
        if not fq.exists():
            raise FileNotFoundError(f"FASTQ not found: {fq}")

    extra_minimap2_args = extra_minimap2_args or []
    extra_samtools_sort_args = extra_samtools_sort_args or []

    producer = [
        "minimap2",
        "-a",
        "-x",
        minimap2_preset,
        "-t",
        str(int(threads)),
        "--secondary=no",
        "--MD",
        str(ref_fa),
    ] + [str(p) for p in fastq_paths] + list(map(str, extra_minimap2_args))

    tmp_prefix = str(out_bam.with_suffix("")) + ".tmp"
    consumer = [
        "samtools",
        "sort",
        "-@",
        str(int(threads)),
        "-m",
        str(sort_mem),
        "-T",
        tmp_prefix,
        "-o",
        str(out_bam),
        "-",
    ] + list(map(str, extra_samtools_sort_args))

    index_cmd = ["samtools", "index", "-@", str(int(threads)), str(out_bam)]

    return {
        "minimap2": producer,
        "samtools_sort": consumer,
        "samtools_index": index_cmd,
    }


def align_ont_fastq_to_bam(
    *,
    fastq: Sequence[str | Path],
    ref_fa: str | Path,
    out_bam: str | Path,
    threads: int = 8,
    minimap2_preset: str = "lr:hq",
    sort_mem: str = "2G",
    extra_minimap2_args: Optional[List[str]] = None,
    extra_samtools_sort_args: Optional[List[str]] = None,
    dry_run: bool = False,
    resume: bool = False,
) -> Dict[str, object]:
    """Align Oxford Nanopore reads (FASTQ) to a reference and create sorted+indexed BAM.

    This function is intended for *convenience* and reproducibility.
    It uses minimap2 + samtools sort (coordinate sort).

    Parameters
    ----------
    fastq:
        One or more FASTQ/FASTQ.GZ files.
    ref_fa:
        Reference FASTA (e.g., GRCh38). Must be the same reference used for VCF calling.
    out_bam:
        Output BAM path (will be coordinate-sorted).
    threads:
        Threads for minimap2 and samtools.
    minimap2_preset:
        Minimap2 preset. Recommended for modern Q20+ ONT WGS is ``lr:hq``.
        Alternatives: ``map-ont`` for noisier reads.
    sort_mem:
        ``samtools sort -m`` memory per thread.
    extra_minimap2_args:
        Extra arguments appended to minimap2.
    extra_samtools_sort_args:
        Extra arguments appended to samtools sort.
    dry_run:
        If True, validate inputs and return planned commands without executing.
    resume:
        If True, skip alignment if out_bam and its index already exist.

    Returns
    -------
    dict
        Summary including runtime and executed commands.
    """
    t0 = time.time()

    ref_fa = Path(ref_fa).expanduser().resolve()
    out_bam = Path(out_bam).expanduser().resolve()
    fastq_paths = [Path(x).expanduser().resolve() for x in fastq]

    if not ref_fa.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {ref_fa}")
    if len(fastq_paths) == 0:
        raise ValueError("At least one FASTQ file is required")
    for fq in fastq_paths:
        if not fq.exists():
            raise FileNotFoundError(f"FASTQ not found: {fq}")

    if minimap2_preset not in _SUPPORTED_PRESETS:
        raise ValueError(
            f"Unsupported minimap2_preset='{minimap2_preset}'. Supported: {sorted(_SUPPORTED_PRESETS)}"
        )

    bai1 = out_bam.with_suffix(out_bam.suffix + ".bai")
    bai2 = out_bam.with_suffix(".bai")
    if resume and out_bam.exists() and (bai1.exists() or bai2.exists()):
        logger.info("Resume enabled: alignment outputs already exist: %s", out_bam)
        return {
            "ref_fa": str(ref_fa),
            "fastq": [str(p) for p in fastq_paths],
            "out_bam": str(out_bam),
            "threads": int(threads),
            "minimap2_preset": minimap2_preset,
            "sort_mem": sort_mem,
            "cmd_minimap2": None,
            "cmd_samtools_sort": None,
            "cmd_samtools_index": None,
            "runtime_seconds": 0.0,
            "skipped": True,
        }

    ensure_executable_in_path(
        "minimap2",
        hint=(
            "Install minimap2. Ubuntu: sudo apt-get install -y minimap2\n"
            "Conda/mamba: mamba install -c bioconda minimap2"
        ),
    )
    ensure_executable_in_path(
        "samtools",
        hint=(
            "Install samtools. Ubuntu: sudo apt-get install -y samtools\n"
            "Conda/mamba: mamba install -c bioconda samtools"
        ),
    )

    cmds = build_alignment_commands(
        fastq=fastq_paths,
        ref_fa=ref_fa,
        out_bam=out_bam,
        threads=threads,
        minimap2_preset=minimap2_preset,
        sort_mem=sort_mem,
        extra_minimap2_args=extra_minimap2_args,
        extra_samtools_sort_args=extra_samtools_sort_args,
    )

    if dry_run:
        return {
            "ref_fa": str(ref_fa),
            "fastq": [str(p) for p in fastq_paths],
            "out_bam": str(out_bam),
            "threads": int(threads),
            "minimap2_preset": minimap2_preset,
            "sort_mem": sort_mem,
            "cmd_minimap2": cmds["minimap2"],
            "cmd_samtools_sort": cmds["samtools_sort"],
            "cmd_samtools_index": cmds["samtools_index"],
            "runtime_seconds": 0.0,
            "dry_run": True,
        }

    ensure_outdir(out_bam.parent)
    _ensure_faidx(ref_fa)

    logger.info("Aligning with minimap2 (%s) and sorting with samtools...", minimap2_preset)
    run_pipe(cmds["minimap2"], cmds["samtools_sort"], check=True)

    logger.info("Indexing BAM: %s", out_bam)
    run_command(cmds["samtools_index"], check=True)

    dt = time.time() - t0

    return {
        "ref_fa": str(ref_fa),
        "fastq": [str(p) for p in fastq_paths],
        "out_bam": str(out_bam),
        "threads": int(threads),
        "minimap2_preset": minimap2_preset,
        "sort_mem": sort_mem,
        "cmd_minimap2": cmds["minimap2"],
        "cmd_samtools_sort": cmds["samtools_sort"],
        "cmd_samtools_index": cmds["samtools_index"],
        "runtime_seconds": float(dt),
    }
