from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Optional

from . import __version__
from .anchors import build_anchor_index, load_anchor_variants
from .assigner import assign_bam
from .doctor import collect_checks
from .plotting import (
    plot_class_counts,
    plot_evidence_hist,
    plot_posterior_hist,
    plot_uncertainty_hist,
)
from .report import render_report
from .utils import ensure_outdir, write_json

# Optional nanopore helpers
from .nanopore import align_ont_fastq_to_bam, call_somatic_vcf, nanopore_endtoend


def _setup_logging(verbosity: int, *, logfile: Optional[Path] = None) -> None:
    level = logging.WARNING
    if verbosity == 1:
        level = logging.INFO
    elif verbosity >= 2:
        level = logging.DEBUG

    log_fmt = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
    logging.basicConfig(level=level, format=log_fmt, stream=sys.stderr)

    if logfile is not None:
        logfile.parent.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(logfile)
        fh.setLevel(level)
        fh.setFormatter(logging.Formatter(log_fmt))
        logging.getLogger().addHandler(fh)


def _path_exists(p: str) -> str:
    if not Path(p).exists():
        raise argparse.ArgumentTypeError(f"Path does not exist: {p}")
    return p


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="gbmreadsep",
        description=(
            "GBMReadSep: Bayesian read-level tumor vs TME assignment for bulk DNA-seq BAMs. "
            "Includes optional Oxford Nanopore WGS helpers (alignment + somatic VCF calling)."
        ),
    )
    p.add_argument("--version", action="version", version=f"gbmreadsep {__version__}")

    sub = p.add_subparsers(dest="cmd", required=True)

    # -----------------
    # assign
    # -----------------
    a = sub.add_parser(
        "assign",
        help="Assign reads in a BAM to tumor vs non-tumor using somatic SNV anchors (VCF).",
    )
    a.add_argument("--bam", required=True, type=_path_exists, help="Input BAM (sorted, indexed).")
    a.add_argument("--vcf", required=True, type=_path_exists, help="Somatic VCF (.vcf/.vcf.gz).")
    a.add_argument(
        "--sample",
        default=None,
        help="VCF sample name to use (default: first sample in VCF).",
    )
    a.add_argument("--outdir", required=True, help="Output directory.")
    a.add_argument(
        "--purity",
        type=float,
        default=None,
        help="Tumor purity prior (0-1). If omitted, estimated heuristically from VAFs.",
    )
    a.add_argument(
        "--min-dp",
        type=int,
        default=10,
        help="Minimum locus depth (DP) in VCF to keep an anchor (if DP available).",
    )
    a.add_argument(
        "--no-require-pass",
        action="store_true",
        help="Do not require FILTER=PASS for anchors.",
    )
    a.add_argument(
        "--max-anchors",
        type=int,
        default=200_000,
        help="Maximum number of anchor SNVs to use (downselects if exceeded).",
    )
    a.add_argument(
        "--default-f-tumor",
        type=float,
        default=0.5,
        help="Fallback tumor ALT fraction for anchors without VAF in VCF.",
    )

    # Model thresholds
    a.add_argument("--tumor-threshold", type=float, default=0.9, help="Posterior >= => class T.")
    a.add_argument("--normal-threshold", type=float, default=0.1, help="Posterior <= => class N.")
    a.add_argument("--min-baseq", type=int, default=20, help="Minimum baseQ for anchor evidence.")

    # Read filters
    a.add_argument("--keep-duplicates", action="store_true", help="Do not skip duplicate reads.")
    a.add_argument("--include-secondary", action="store_true", help="Include secondary alignments.")
    a.add_argument(
        "--include-supplementary", action="store_true", help="Include supplementary alignments."
    )

    # Outputs
    a.add_argument(
        "--write-tagged-bam",
        default=None,
        help="Write a BAM with assignment tags (TP/TC/TE/TL).",
    )
    a.add_argument(
        "--split-bams",
        action="store_true",
        help="Write tumor.bam / tme.bam / uncertain.bam into outdir.",
    )
    a.add_argument(
        "--assignments-tsv",
        default=None,
        help="Optional path for assignments TSV.GZ (default: outdir/assignments.tsv.gz).",
    )

    a.add_argument("-v", "--verbose", action="count", default=0, help="Increase verbosity (-v/-vv).")

    # -----------------
    # doctor
    # -----------------
    d = sub.add_parser(
        "doctor",
        help="Check your environment for required external tools (minimap2/samtools/docker).",
    )
    d.add_argument("-v", "--verbose", action="count", default=0, help="Increase verbosity (-v/-vv).")

    # -----------------
    # nanopore group
    # -----------------
    n = sub.add_parser(
        "nanopore",
        help="Nanopore helpers: FASTQ->BAM alignment, somatic VCF calling, and end-to-end pipeline.",
    )
    nsub = n.add_subparsers(dest="nanocmd", required=True)

    # nanopore align
    na = nsub.add_parser("align", help="Align ONT FASTQ to reference and create sorted+indexed BAM.")
    na.add_argument(
        "--fastq",
        required=True,
        nargs="+",
        type=_path_exists,
        help="One or more FASTQ/FASTQ.GZ files.",
    )
    na.add_argument("--ref", required=True, type=_path_exists, help="Reference FASTA.")
    na.add_argument("--out-bam", required=True, help="Output BAM path.")
    na.add_argument("--threads", type=int, default=8, help="Threads for minimap2/samtools.")
    na.add_argument(
        "--preset",
        type=str,
        default="lr:hq",
        help="Minimap2 preset (recommended for Q20+ ONT: lr:hq; alternative: map-ont).",
    )
    na.add_argument(
        "--sort-mem",
        type=str,
        default="2G",
        help="samtools sort -m (memory per thread, e.g. 2G).",
    )
    na.add_argument("-v", "--verbose", action="count", default=0, help="Increase verbosity (-v/-vv).")

    # nanopore call-vcf
    nv = nsub.add_parser(
        "call-vcf",
        help="Call somatic variants for ONT WGS from BAM(s) using ClairS / ClairS-TO.",
    )
    nv.add_argument("--tumor-bam", required=True, type=_path_exists, help="Tumor BAM (sorted, indexed).")
    nv.add_argument(
        "--normal-bam",
        default=None,
        type=_path_exists,
        help="Optional matched normal BAM (enables ClairS paired calling).",
    )
    nv.add_argument("--ref", required=True, type=_path_exists, help="Reference FASTA (same as used for alignment).")
    nv.add_argument("--outdir", required=True, help="Output directory for VCF(s).")
    nv.add_argument(
        "--platform",
        default="ont_r10_dorado_sup_5khz",
        help=(
            "Sequencing/basecalling platform identifier expected by ClairS/ClairS-TO. "
            "Example: ont_r10_dorado_sup_5khz"
        ),
    )
    nv.add_argument("--threads", type=int, default=8, help="Threads for the caller.")
    nv.add_argument(
        "--engine",
        default="docker",
        choices=["docker", "native"],
        help="Execution engine (docker recommended).",
    )
    nv.add_argument(
        "--sample-name",
        default="SAMPLE",
        help="Sample name to embed in the VCF.",
    )
    nv.add_argument(
        "--snv-min-af",
        type=float,
        default=None,
        help="Optional: ClairS/ClairS-TO --snv_min_af.",
    )
    nv.add_argument(
        "--snv-min-qual",
        type=float,
        default=None,
        help="Optional: ClairS/ClairS-TO --snv_min_qual.",
    )
    nv.add_argument(
        "--indel-min-qual",
        type=float,
        default=None,
        help="Optional: ClairS/ClairS-TO --indel_min_qual.",
    )
    nv.add_argument(
        "--min-coverage",
        type=int,
        default=None,
        help="Optional: ClairS/ClairS-TO --min_coverage.",
    )
    nv.add_argument(
        "--bed",
        default=None,
        type=_path_exists,
        help="Optional BED file restricting calling regions.",
    )
    nv.add_argument(
        "--region",
        default=None,
        help="Optional region string ctg:start-end.",
    )
    nv.add_argument(
        "--ctg-name",
        default=None,
        help="Optional contig(s) to process (comma-separated).",
    )
    nv.add_argument("-v", "--verbose", action="count", default=0, help="Increase verbosity (-v/-vv).")

    # nanopore endtoend
    ne = nsub.add_parser(
        "endtoend",
        help="End-to-end Nanopore workflow: FASTQ/BAM -> somatic VCF -> read assignments + report.",
    )
    group_t = ne.add_mutually_exclusive_group(required=True)
    group_t.add_argument(
        "--tumor-fastq",
        nargs="+",
        type=_path_exists,
        help="Tumor FASTQ/FASTQ.GZ file(s).",
    )
    group_t.add_argument(
        "--tumor-bam",
        type=_path_exists,
        help="Tumor BAM (skip alignment).",
    )

    group_n = ne.add_mutually_exclusive_group(required=False)
    group_n.add_argument(
        "--normal-fastq",
        nargs="+",
        type=_path_exists,
        help="Normal FASTQ/FASTQ.GZ file(s) (enables paired calling).",
    )
    group_n.add_argument(
        "--normal-bam",
        type=_path_exists,
        help="Normal BAM (skip alignment; enables paired calling).",
    )

    ne.add_argument("--ref", required=True, type=_path_exists, help="Reference FASTA.")
    ne.add_argument("--outdir", required=True, help="Output directory.")
    ne.add_argument("--threads", type=int, default=8, help="Threads for alignment/calling.")
    ne.add_argument(
        "--preset",
        type=str,
        default="lr:hq",
        help="Minimap2 preset for alignment (if FASTQ provided).",
    )
    ne.add_argument(
        "--platform",
        default="ont_r10_dorado_sup_5khz",
        help="Platform string for ClairS/ClairS-TO.",
    )
    ne.add_argument(
        "--engine",
        default="docker",
        choices=["docker", "native"],
        help="Variant calling execution engine (docker recommended).",
    )

    # Read separation tuning
    ne.add_argument(
        "--purity",
        type=float,
        default=None,
        help="Tumor purity prior for read separation (if omitted, estimated from VAFs).",
    )
    ne.add_argument("--tumor-threshold", type=float, default=0.9, help="Posterior >= => class T.")
    ne.add_argument("--normal-threshold", type=float, default=0.1, help="Posterior <= => class N.")
    ne.add_argument(
        "--min-baseq",
        type=int,
        default=10,
        help="Minimum baseQ for anchor evidence (ONT default is lower than Illumina).",
    )
    ne.add_argument(
        "--keep-duplicates",
        action="store_true",
        help="Do not skip duplicate reads during assignment.",
    )
    ne.add_argument(
        "--no-split-bams",
        action="store_true",
        help="Do not write tumor/tme/uncertain BAMs.",
    )
    ne.add_argument(
        "--no-tagged-bam",
        action="store_true",
        help="Do not write a tagged BAM.",
    )

    ne.add_argument("-v", "--verbose", action="count", default=0, help="Increase verbosity (-v/-vv).")

    return p


# -----------------
# Command handlers
# -----------------

def cmd_assign(args: argparse.Namespace) -> int:
    outdir = ensure_outdir(args.outdir)
    _setup_logging(args.verbose, logfile=Path(outdir) / "gbmreadsep.log")

    logger = logging.getLogger("gbmreadsep")
    logger.info("gbmreadsep %s", __version__)

    anchors, purity_used, anchor_stats = load_anchor_variants(
        args.vcf,
        sample=args.sample,
        purity=args.purity,
        min_dp=int(args.min_dp),
        require_pass=not bool(args.no_require_pass),
        max_anchors=int(args.max_anchors),
        default_f_tumor=float(args.default_f_tumor),
    )
    sample_used = args.sample or "AUTO"
    anchor_index = build_anchor_index(anchors)

    run = assign_bam(
        bam_path=args.bam,
        anchor_index_by_contig=anchor_index,
        prior_p_tumor=purity_used,
        outdir=outdir,
        tumor_threshold=float(args.tumor_threshold),
        normal_threshold=float(args.normal_threshold),
        min_baseq=int(args.min_baseq),
        skip_duplicates=not bool(args.keep_duplicates),
        include_secondary=bool(args.include_secondary),
        include_supplementary=bool(args.include_supplementary),
        write_tagged_bam=args.write_tagged_bam,
        split_bams=bool(args.split_bams),
        assignments_tsv_gz=args.assignments_tsv,
        progress=True,
    )

    write_json(Path(outdir) / "anchor_stats.json", anchor_stats)

    plots_dir = Path(outdir) / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    class_counts_png = plots_dir / "class_counts.png"
    posterior_png = plots_dir / "posterior_hist.png"
    evidence_png = plots_dir / "evidence_hist.png"
    uncertainty_png = plots_dir / "uncertainty_hist.png"

    plot_class_counts(class_counts=run["counts"], out_png=class_counts_png)
    plot_posterior_hist(
        bin_edges=run["posterior_hist"]["bin_edges"],
        counts=run["posterior_hist"]["counts"],
        out_png=posterior_png,
    )
    plot_evidence_hist(evidence_hist=run["evidence_count_hist"], out_png=evidence_png)
    plot_uncertainty_hist(posterior_hist=run["posterior_hist"], out_png=uncertainty_png)

    plots_rel = {
        "class_counts": str(Path("plots") / class_counts_png.name),
        "posterior_hist": str(Path("plots") / posterior_png.name),
        "evidence_hist": str(Path("plots") / evidence_png.name),
        "uncertainty_hist": str(Path("plots") / uncertainty_png.name),
    }

    report_path = render_report(
        outdir=outdir,
        version=__version__,
        run=run,
        anchor_stats=anchor_stats,
        vcf_path=args.vcf,
        sample=sample_used,
        plots=plots_rel,
    )

    logger.info("Report written: %s", report_path)
    print(str(report_path))
    return 0


def cmd_doctor(args: argparse.Namespace) -> int:
    _setup_logging(args.verbose, logfile=None)

    checks = collect_checks()

    # Human-readable output
    lines = []
    ok_all = True
    for name in ["python", "minimap2", "samtools", "docker"]:
        r = checks[name]
        status = "OK" if r.ok else "MISSING"
        lines.append(f"{name:9s} : {status:7s}  {r.detail}")
        if not r.ok:
            ok_all = False

    print("\n".join(lines))

    # Guidance
    for name in ["minimap2", "samtools", "docker"]:
        r = checks[name]
        if not r.ok and r.howto:
            print("\n---")
            print(f"How to install/fix '{name}':")
            print(r.howto)

    return 0 if ok_all else 1


def cmd_nanopore_align(args: argparse.Namespace) -> int:
    out_bam = Path(args.out_bam).expanduser().resolve()
    _setup_logging(args.verbose, logfile=out_bam.parent / "gbmreadsep_align.log")

    summary = align_ont_fastq_to_bam(
        fastq=args.fastq,
        ref_fa=args.ref,
        out_bam=out_bam,
        threads=int(args.threads),
        minimap2_preset=str(args.preset),
        sort_mem=str(args.sort_mem),
    )

    # Emit a small JSON next to the BAM
    summary_path = out_bam.with_suffix(".align_summary.json")
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(str(out_bam))
    return 0


def cmd_nanopore_call_vcf(args: argparse.Namespace) -> int:
    outdir = ensure_outdir(args.outdir)
    _setup_logging(args.verbose, logfile=Path(outdir) / "gbmreadsep_call_vcf.log")

    summary = call_somatic_vcf(
        tumor_bam=args.tumor_bam,
        normal_bam=args.normal_bam,
        ref_fa=args.ref,
        outdir=outdir,
        platform=args.platform,
        threads=int(args.threads),
        sample_name=str(args.sample_name),
        engine=str(args.engine),
        snv_min_af=args.snv_min_af,
        snv_min_qual=args.snv_min_qual,
        indel_min_qual=args.indel_min_qual,
        min_coverage=args.min_coverage,
        bed_fn=args.bed,
        region=args.region,
        ctg_name=args.ctg_name,
    )

    write_json(Path(outdir) / "call_vcf_summary.json", summary)

    # Print the primary VCF path (snv for tumor-only; output for paired)
    if summary.get("caller") == "clairs":
        print(str(summary.get("out_vcf")))
    else:
        print(str(summary.get("out_snv_vcf")))
    return 0


def cmd_nanopore_endtoend(args: argparse.Namespace) -> int:
    outdir = ensure_outdir(args.outdir)
    _setup_logging(args.verbose, logfile=Path(outdir) / "gbmreadsep_pipeline.log")

    report_path = nanopore_endtoend(
        tumor_fastq=args.tumor_fastq,
        normal_fastq=args.normal_fastq,
        tumor_bam=args.tumor_bam,
        normal_bam=args.normal_bam,
        ref_fa=args.ref,
        outdir=outdir,
        threads=int(args.threads),
        minimap2_preset=str(args.preset),
        platform=str(args.platform),
        caller_engine=str(args.engine),
        purity=args.purity,
        tumor_threshold=float(args.tumor_threshold),
        normal_threshold=float(args.normal_threshold),
        min_baseq=int(args.min_baseq),
        skip_duplicates=not bool(args.keep_duplicates),
        split_bams=not bool(args.no_split_bams),
        write_tagged_bam=not bool(args.no_tagged_bam),
    )

    print(str(report_path))
    return 0


def main(argv: Optional[list[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.cmd == "assign":
        return cmd_assign(args)
    if args.cmd == "doctor":
        return cmd_doctor(args)
    if args.cmd == "nanopore":
        if args.nanocmd == "align":
            return cmd_nanopore_align(args)
        if args.nanocmd == "call-vcf":
            return cmd_nanopore_call_vcf(args)
        if args.nanocmd == "endtoend":
            return cmd_nanopore_endtoend(args)
        parser.error(f"Unknown nanopore subcommand: {args.nanocmd}")

    parser.error(f"Unknown command: {args.cmd}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
