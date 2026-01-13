from __future__ import annotations

import logging
import time
from pathlib import Path
from typing import Dict, List, Mapping, Optional

import pysam

from ..external import DockerMount, build_docker_mounts, docker_run, ensure_executable_in_path, run_command
from ..utils import ensure_outdir

logger = logging.getLogger(__name__)


# Keep these lists aligned with upstream documentation as of late 2025.
# If a new basecalling model appears, we *warn* rather than hard-fail.
_CLAIRS_PLATFORMS = {
    "ont_r10_dorado_sup_4khz",
    "ont_r10_dorado_sup_5khz",
    "ont_r10_dorado_hac_5khz",
    "ont_r10_dorado_hac_4khz",
    "ont_r10_guppy",
    "ont_r9_guppy",
    "ilmn",
    "hifi_sequel2",
    "hifi_revio",
}

_CLAIRSTO_PLATFORMS = {
    "ont_r10_dorado_sup_4khz",
    "ont_r10_dorado_hac_4khz",
    "ont_r10_dorado_sup_5khz",
    "ont_r10_dorado_sup_5khz_ss",
    "ont_r10_dorado_sup_5khz_ssrs",
    "ont_r10_guppy_sup_4khz",
    "ont_r10_guppy_hac_5khz",
    "ilmn",
    "ilmn_ss",
    "ilmn_ssrs",
    "hifi_revio",
    "hifi_revio_ss",
    "hifi_revio_ssrs",
}


def _ensure_bam_index(bam: Path) -> None:
    """Ensure BAM has an index file (.bai)."""
    bai1 = bam.with_suffix(bam.suffix + ".bai")
    bai2 = bam.with_suffix(".bai")
    if bai1.exists() or bai2.exists():
        return
    raise ValueError(f"BAM is not indexed. Run: samtools index {bam}")


def _ensure_faidx(ref_fa: Path) -> None:
    fai = ref_fa.with_suffix(ref_fa.suffix + ".fai")
    if fai.exists():
        return
    logger.info("Creating FASTA index: %s", fai)
    pysam.faidx(str(ref_fa))


def _container_path_for_file(host_file: Path, dir_map: Mapping[Path, str]) -> str:
    d = host_file.parent
    if d not in dir_map:
        raise KeyError(f"Directory {d} was not mounted into the container")
    return f"{dir_map[d]}/{host_file.name}"


def _warn_if_unknown_platform(platform: str, known: set[str], caller_name: str) -> None:
    if platform not in known:
        logger.warning(
            "Platform '%s' is not in the known %s platform list. "
            "Proceeding anyway; if the caller rejects it, re-run with a valid platform.",
            platform,
            caller_name,
        )


def _build_docker_cmd(
    *,
    image: str,
    argv: List[str],
    mounts: List[DockerMount],
) -> List[str]:
    cmd: List[str] = ["docker", "run", "--rm"]
    for m in mounts:
        ro = ":ro" if m.read_only else ""
        cmd.extend(["-v", f"{str(m.host_dir)}:{m.container_dir}{ro}"])
    cmd.append(image)
    cmd.extend(argv)
    return cmd


def build_somatic_commands(
    *,
    tumor_bam: str | Path,
    ref_fa: str | Path,
    outdir: str | Path,
    normal_bam: Optional[str | Path] = None,
    platform: str = "ont_r10_dorado_sup_5khz",
    threads: int = 8,
    sample_name: str = "SAMPLE",
    engine: str = "docker",
    docker_image_clairs: str = "hkubal/clairs:latest",
    docker_image_clairsto: str = "hkubal/clairs-to:latest",
    snv_min_af: Optional[float] = None,
    snv_min_qual: Optional[float] = None,
    indel_min_qual: Optional[float] = None,
    min_coverage: Optional[int] = None,
    bed_fn: Optional[str | Path] = None,
    region: Optional[str] = None,
    ctg_name: Optional[str] = None,
    extra_args: Optional[List[str]] = None,
) -> Dict[str, object]:
    """Build the command(s) for somatic VCF calling."""
    tumor_bam = Path(tumor_bam).expanduser().resolve()
    ref_fa = Path(ref_fa).expanduser().resolve()
    outdir = Path(outdir).expanduser().resolve()
    normal_bam_p: Optional[Path] = None
    if normal_bam is not None:
        normal_bam_p = Path(normal_bam).expanduser().resolve()

    if not tumor_bam.exists():
        raise FileNotFoundError(f"Tumor BAM not found: {tumor_bam}")
    if normal_bam_p is not None and not normal_bam_p.exists():
        raise FileNotFoundError(f"Normal BAM not found: {normal_bam_p}")
    if not ref_fa.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {ref_fa}")

    extra_args = extra_args or []

    opt_args: List[str] = []
    if ctg_name:
        opt_args += ["--ctg_name", str(ctg_name)]
    if region:
        opt_args += ["--region", str(region)]
    if bed_fn:
        bed_p = Path(bed_fn).expanduser().resolve()
        if not bed_p.exists():
            raise FileNotFoundError(f"BED file not found: {bed_p}")
        opt_args += ["--bed_fn", str(bed_p)]
    if snv_min_af is not None:
        opt_args += ["--snv_min_af", str(float(snv_min_af))]
    if snv_min_qual is not None:
        opt_args += ["--snv_min_qual", str(float(snv_min_qual))]
    if indel_min_qual is not None:
        opt_args += ["--indel_min_qual", str(float(indel_min_qual))]
    if min_coverage is not None:
        opt_args += ["--min_coverage", str(int(min_coverage))]

    if normal_bam_p is not None:
        out_vcf = outdir / "output.vcf.gz"
        if engine == "docker":
            mounts, dir_map = build_docker_mounts([tumor_bam, normal_bam_p, ref_fa])
            mounts.append(
                DockerMount(host_dir=outdir.resolve(), container_dir="/mnt/gbmreadsep/out", read_only=False)
            )
            argv = [
                "/opt/bin/run_clairs",
                "--tumor_bam_fn",
                _container_path_for_file(tumor_bam, dir_map),
                "--normal_bam_fn",
                _container_path_for_file(normal_bam_p, dir_map),
                "--ref_fn",
                _container_path_for_file(ref_fa, dir_map),
                "--threads",
                str(int(threads)),
                "--platform",
                str(platform),
                "--output_dir",
                "/mnt/gbmreadsep/out",
                "--sample_name",
                str(sample_name),
            ] + opt_args + list(map(str, extra_args))
            return {
                "caller": "clairs",
                "engine": engine,
                "out_vcf": str(out_vcf),
                "docker_image": docker_image_clairs,
                "docker_cmd": _build_docker_cmd(image=docker_image_clairs, argv=argv, mounts=mounts),
                "docker_mounts": mounts,
                "docker_argv": argv,
            }
        argv = [
            "run_clairs",
            "--tumor_bam_fn",
            str(tumor_bam),
            "--normal_bam_fn",
            str(normal_bam_p),
            "--ref_fn",
            str(ref_fa),
            "--threads",
            str(int(threads)),
            "--platform",
            str(platform),
            "--output_dir",
            str(outdir),
            "--sample_name",
            str(sample_name),
        ] + opt_args + list(map(str, extra_args))
        return {
            "caller": "clairs",
            "engine": engine,
            "out_vcf": str(out_vcf),
            "cmd": argv,
        }

    out_snv_vcf = outdir / "snv.vcf.gz"
    out_indel_vcf = outdir / "indel.vcf.gz"

    if engine == "docker":
        mounts, dir_map = build_docker_mounts([tumor_bam, ref_fa])
        mounts.append(
            DockerMount(host_dir=outdir.resolve(), container_dir="/mnt/gbmreadsep/out", read_only=False)
        )
        argv = [
            "/opt/bin/run_clairs_to",
            "--tumor_bam_fn",
            _container_path_for_file(tumor_bam, dir_map),
            "--ref_fn",
            _container_path_for_file(ref_fa, dir_map),
            "--threads",
            str(int(threads)),
            "--platform",
            str(platform),
            "--output_dir",
            "/mnt/gbmreadsep/out",
            "--sample_name",
            str(sample_name),
        ] + opt_args + list(map(str, extra_args))
        return {
            "caller": "clairs-to",
            "engine": engine,
            "out_snv_vcf": str(out_snv_vcf),
            "out_indel_vcf": str(out_indel_vcf),
            "docker_image": docker_image_clairsto,
            "docker_cmd": _build_docker_cmd(image=docker_image_clairsto, argv=argv, mounts=mounts),
            "docker_mounts": mounts,
            "docker_argv": argv,
        }

    argv = [
        "run_clairs_to",
        "--tumor_bam_fn",
        str(tumor_bam),
        "--ref_fn",
        str(ref_fa),
        "--threads",
        str(int(threads)),
        "--platform",
        str(platform),
        "--output_dir",
        str(outdir),
        "--sample_name",
        str(sample_name),
    ] + opt_args + list(map(str, extra_args))
    return {
        "caller": "clairs-to",
        "engine": engine,
        "out_snv_vcf": str(out_snv_vcf),
        "out_indel_vcf": str(out_indel_vcf),
        "cmd": argv,
    }


def call_somatic_vcf(
    *,
    tumor_bam: str | Path,
    ref_fa: str | Path,
    outdir: str | Path,
    normal_bam: Optional[str | Path] = None,
    platform: str = "ont_r10_dorado_sup_5khz",
    threads: int = 8,
    sample_name: str = "SAMPLE",
    engine: str = "docker",
    docker_image_clairs: str = "hkubal/clairs:latest",
    docker_image_clairsto: str = "hkubal/clairs-to:latest",
    snv_min_af: Optional[float] = None,
    snv_min_qual: Optional[float] = None,
    indel_min_qual: Optional[float] = None,
    min_coverage: Optional[int] = None,
    bed_fn: Optional[str | Path] = None,
    region: Optional[str] = None,
    ctg_name: Optional[str] = None,
    extra_args: Optional[List[str]] = None,
    dry_run: bool = False,
    resume: bool = False,
) -> Dict[str, object]:
    """Call somatic variants from ONT long-read alignments.

    This is a pragmatic wrapper around:

    - **ClairS** for paired tumor/normal somatic calling
    - **ClairS-TO** for tumor-only somatic calling

    Both tools support Docker, Singularity, and conda installs.

    Parameters
    ----------
    tumor_bam:
        Tumor BAM (sorted, indexed).
    normal_bam:
        Optional matched normal BAM.
    ref_fa:
        Reference FASTA (must match alignment).
    outdir:
        Output directory for VCF(s).
    platform:
        Basecalling / platform string expected by ClairS / ClairS-TO.
    engine:
        "docker" (recommended) or "native".
    dry_run:
        If True, validate inputs and return planned commands without executing.
    resume:
        If True, skip calling if the expected output VCF already exists.
    docker_image_clairs, docker_image_clairsto:
        Docker images for ClairS and ClairS-TO.

    Returns
    -------
    dict
        Summary including output VCF paths.
    """
    t0 = time.time()

    tumor_bam = Path(tumor_bam).expanduser().resolve()
    ref_fa = Path(ref_fa).expanduser().resolve()
    outdir = Path(outdir).expanduser().resolve()
    normal_bam_p: Optional[Path] = None
    if normal_bam is not None:
        normal_bam_p = Path(normal_bam).expanduser().resolve()

    if not tumor_bam.exists():
        raise FileNotFoundError(f"Tumor BAM not found: {tumor_bam}")
    if normal_bam_p is not None and not normal_bam_p.exists():
        raise FileNotFoundError(f"Normal BAM not found: {normal_bam_p}")
    if not ref_fa.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {ref_fa}")

    if bed_fn:
        bed_p = Path(bed_fn).expanduser().resolve()
        if not bed_p.exists():
            raise FileNotFoundError(f"BED file not found: {bed_p}")

    if normal_bam_p is None:
        logger.warning(
            "Tumor-only somatic calling without a matched normal or germline filtering can "
            "retain germline variants. Consider providing a normal or filtering downstream."
        )

    if engine == "docker":
        ensure_executable_in_path(
            "docker",
            hint=(
                "Install Docker. Ubuntu: sudo apt-get install -y docker.io\n"
                "Conda: use system Docker or Docker Desktop."
            ),
        )
    else:
        exe = "run_clairs" if normal_bam_p is not None else "run_clairs_to"
        ensure_executable_in_path(
            exe,
            hint=(
                "Install ClairS/ClairS-TO (recommended via Docker) or ensure the caller is in PATH.\n"
                "Upstream: https://github.com/HKU-BAL/ClairS"
            ),
        )

    if resume:
        if normal_bam_p is not None:
            out_vcf = outdir / "output.vcf.gz"
            if out_vcf.exists():
                logger.info("Resume enabled: somatic VCF already exists: %s", out_vcf)
                return {
                    "caller": "clairs",
                    "engine": engine,
                    "tumor_bam": str(tumor_bam),
                    "normal_bam": str(normal_bam_p),
                    "ref_fa": str(ref_fa),
                    "outdir": str(outdir),
                    "platform": platform,
                    "threads": int(threads),
                    "sample_name": sample_name,
                    "out_vcf": str(out_vcf),
                    "runtime_seconds": 0.0,
                    "skipped": True,
                }
        else:
            out_snv_vcf = outdir / "snv.vcf.gz"
            if out_snv_vcf.exists():
                logger.info("Resume enabled: somatic VCF already exists: %s", out_snv_vcf)
                return {
                    "caller": "clairs-to",
                    "engine": engine,
                    "tumor_bam": str(tumor_bam),
                    "normal_bam": None,
                    "ref_fa": str(ref_fa),
                    "outdir": str(outdir),
                    "platform": platform,
                    "threads": int(threads),
                    "sample_name": sample_name,
                    "out_snv_vcf": str(out_snv_vcf),
                    "runtime_seconds": 0.0,
                    "skipped": True,
                }

    if not dry_run:
        outdir = ensure_outdir(outdir)

    cmds = build_somatic_commands(
        tumor_bam=tumor_bam,
        normal_bam=normal_bam_p,
        ref_fa=ref_fa,
        outdir=outdir,
        platform=platform,
        threads=threads,
        sample_name=sample_name,
        engine=engine,
        docker_image_clairs=docker_image_clairs,
        docker_image_clairsto=docker_image_clairsto,
        snv_min_af=snv_min_af,
        snv_min_qual=snv_min_qual,
        indel_min_qual=indel_min_qual,
        min_coverage=min_coverage,
        bed_fn=bed_fn,
        region=region,
        ctg_name=ctg_name,
        extra_args=extra_args,
    )

    if dry_run:
        _ensure_bam_index(tumor_bam)
        if normal_bam_p is not None:
            _ensure_bam_index(normal_bam_p)
        cmds["runtime_seconds"] = 0.0
        cmds["dry_run"] = True
        return cmds

    _ensure_bam_index(tumor_bam)
    if normal_bam_p is not None:
        _ensure_bam_index(normal_bam_p)
    _ensure_faidx(ref_fa)

    if normal_bam_p is not None:
        _warn_if_unknown_platform(platform, _CLAIRS_PLATFORMS, "ClairS")
        if engine == "docker":
            logger.info("Calling somatic variants with ClairS (docker): %s", docker_image_clairs)
            docker_run(
                image=docker_image_clairs,
                argv=cmds["docker_argv"],
                mounts=cmds["docker_mounts"],
            )
        else:
            logger.info("Calling somatic variants with ClairS (native)...")
            run_command(cmds["cmd"], check=True)

        dt = time.time() - t0
        return {
            "caller": "clairs",
            "engine": engine,
            "tumor_bam": str(tumor_bam),
            "normal_bam": str(normal_bam_p),
            "ref_fa": str(ref_fa),
            "outdir": str(outdir),
            "platform": platform,
            "threads": int(threads),
            "sample_name": sample_name,
            "out_vcf": cmds["out_vcf"],
            "runtime_seconds": float(dt),
        }

    _warn_if_unknown_platform(platform, _CLAIRSTO_PLATFORMS, "ClairS-TO")

    if engine == "docker":
        logger.info("Calling somatic variants with ClairS-TO (docker): %s", docker_image_clairsto)
        docker_run(
            image=docker_image_clairsto,
            argv=cmds["docker_argv"],
            mounts=cmds["docker_mounts"],
        )
    else:
        logger.info("Calling somatic variants with ClairS-TO (native)...")
        run_command(cmds["cmd"], check=True)

    dt = time.time() - t0
    return {
        "caller": "clairs-to",
        "engine": engine,
        "tumor_bam": str(tumor_bam),
        "normal_bam": None,
        "ref_fa": str(ref_fa),
        "outdir": str(outdir),
        "platform": platform,
        "threads": int(threads),
        "sample_name": sample_name,
        "out_snv_vcf": cmds["out_snv_vcf"],
        "out_indel_vcf": cmds["out_indel_vcf"],
        "runtime_seconds": float(dt),
    }
