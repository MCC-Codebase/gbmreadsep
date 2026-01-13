from __future__ import annotations

import datetime as _dt
import json
import logging
import time
from pathlib import Path
from typing import Any, Dict, Optional, Sequence

from contextlib import contextmanager

from jinja2 import Template

from .. import __version__
from ..anchors import build_anchor_index, load_anchor_variants
from ..assigner import assign_bam
from ..plotting import plot_class_counts, plot_evidence_hist, plot_posterior_hist, plot_uncertainty_hist
from ..report import render_report
from ..utils import ensure_outdir, write_json
from .align import align_ont_fastq_to_bam
from .somatic import call_somatic_vcf

logger = logging.getLogger(__name__)


_LOG_FMT = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"


@contextmanager
def _step_log(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    handler = logging.FileHandler(path)
    handler.setFormatter(logging.Formatter(_LOG_FMT))
    root_logger = logging.getLogger()
    root_logger.addHandler(handler)
    try:
        yield
    finally:
        root_logger.removeHandler(handler)
        handler.close()


_PIPELINE_TEMPLATE = Template(
    """<!doctype html>
<html lang=\"en\">
<head>
  <meta charset=\"utf-8\">
  <title>GBMReadSep Nanopore Pipeline Report</title>
  <style>
    body { font-family: Arial, Helvetica, sans-serif; margin: 24px; }
    code, pre { background: #f6f8fa; padding: 2px 4px; border-radius: 4px; }
    pre { padding: 12px; overflow-x: auto; }
    h1, h2, h3 { margin-top: 1.2em; }
    table { border-collapse: collapse; margin-top: 0.6em; }
    th, td { border: 1px solid #ddd; padding: 8px; }
    th { background: #f2f2f2; text-align: left; }
    .grid { display: grid; grid-template-columns: 1fr 1fr; gap: 16px; }
    .card { border: 1px solid #ddd; border-radius: 8px; padding: 12px; }
    .small { color: #666; font-size: 0.9em; }
  </style>
</head>
<body>

<h1>GBMReadSep Nanopore Pipeline Report</h1>
<p class=\"small\">Generated: {{ generated_at }}</p>

<h2>Overview</h2>
<ul>
  <li>Alignment (FASTQ → BAM): <code>{{ alignment_dir }}</code></li>
  <li>Somatic variant calling (BAM → VCF): <code>{{ vcf_dir }}</code></li>
  <li>Read separation (BAM + VCF → tumor/TME assignments): <code>{{ readsep_dir }}</code></li>
</ul>

<h2>Key outputs</h2>
<ul>
  <li><a href=\"{{ readsep_report_rel }}\">Read separation report</a></li>
  <li><code>{{ pipeline_summary }}</code> (machine-readable pipeline summary)</li>
</ul>

<h2>Alignment</h2>
<div class=\"card\">
  <table>
    <tr><th>Tumor BAM</th><td><code>{{ tumor_bam }}</code></td></tr>
    {% if normal_bam %}
    <tr><th>Normal BAM</th><td><code>{{ normal_bam }}</code></td></tr>
    {% endif %}
    <tr><th>Reference</th><td><code>{{ ref_fa }}</code></td></tr>
    <tr><th>Preset</th><td><code>{{ minimap2_preset }}</code></td></tr>
    <tr><th>Threads</th><td>{{ threads }}</td></tr>
    <tr><th>Runtime</th><td>{{ alignment_runtime }} s</td></tr>
  </table>
</div>

<h2>Somatic variant calling</h2>
<div class=\"card\">
  <table>
    <tr><th>Caller</th><td><code>{{ caller }}</code></td></tr>
    <tr><th>Engine</th><td><code>{{ engine }}</code></td></tr>
    <tr><th>Platform</th><td><code>{{ platform }}</code></td></tr>
    <tr><th>VCF</th><td><code>{{ vcf_path }}</code></td></tr>
    <tr><th>Runtime</th><td>{{ vcf_runtime }} s</td></tr>
  </table>
</div>

<h2>Next step</h2>
<p>
Open the <a href=\"{{ readsep_report_rel }}\">read separation report</a> for uncertainty statistics,
read counts, and per-read assignments.
</p>

<hr>
<p class=\"small\">GBMReadSep {{ version }}</p>
</body>
</html>"""
)


def nanopore_endtoend(
    *,
    tumor_fastq: Optional[Sequence[str | Path]],
    normal_fastq: Optional[Sequence[str | Path]],
    tumor_bam: Optional[str | Path],
    normal_bam: Optional[str | Path],
    ref_fa: str | Path,
    outdir: str | Path,
    threads: int = 8,
    minimap2_preset: str = "lr:hq",
    platform: str = "ont_r10_dorado_sup_5khz",
    caller_engine: str = "docker",
    purity: Optional[float] = None,
    tumor_threshold: float = 0.9,
    normal_threshold: float = 0.1,
    min_baseq: int = 10,
    skip_duplicates: bool = True,
    split_bams: bool = True,
    write_tagged_bam: bool = True,
    dry_run: bool = False,
    resume: bool = False,
) -> Path:
    """Run FASTQ/BAM → somatic VCF → read separation.

    Returns
    -------
    Path
        Path to the top-level pipeline HTML report.
    """
    t0 = time.time()
    outdir_p = Path(outdir).expanduser().resolve()

    alignment_dir = outdir_p / "alignment"
    vcf_dir = outdir_p / "somatic_vcf"
    readsep_dir = outdir_p / "readsep"

    ref_fa_p = Path(ref_fa).expanduser().resolve()

    if dry_run:
        logger.info("Dry-run: building planned commands for alignment and calling.")

    # 1) Alignment (optional)
    aln_summary: Dict[str, Any] = {}

    if tumor_bam is None:
        if not tumor_fastq:
            raise ValueError("Provide either tumor_bam or tumor_fastq")
        tumor_bam_p = alignment_dir / "tumor.bam"
        with _step_log(outdir_p / "logs" / "alignment.log"):
            aln_summary["tumor"] = align_ont_fastq_to_bam(
                fastq=tumor_fastq,
                ref_fa=ref_fa_p,
                out_bam=tumor_bam_p,
                threads=threads,
                minimap2_preset=minimap2_preset,
                dry_run=dry_run,
                resume=resume,
            )
    else:
        tumor_bam_p = Path(tumor_bam).expanduser().resolve()

    normal_bam_p: Optional[Path] = None
    if normal_bam is not None:
        normal_bam_p = Path(normal_bam).expanduser().resolve()
    elif normal_fastq:
        normal_bam_p = alignment_dir / "normal.bam"
        with _step_log(outdir_p / "logs" / "alignment.log"):
            aln_summary["normal"] = align_ont_fastq_to_bam(
                fastq=normal_fastq,
                ref_fa=ref_fa_p,
                out_bam=normal_bam_p,
                threads=threads,
                minimap2_preset=minimap2_preset,
                dry_run=dry_run,
                resume=resume,
            )

    # 2) Somatic VCF calling
    with _step_log(outdir_p / "logs" / "caller.log"):
        vcf_summary = call_somatic_vcf(
            tumor_bam=tumor_bam_p,
            normal_bam=normal_bam_p,
            ref_fa=ref_fa_p,
            outdir=vcf_dir,
            threads=threads,
            platform=platform,
            engine=caller_engine,
            sample_name="TUMOR" if normal_bam_p is not None else "SAMPLE",
            dry_run=dry_run,
            resume=resume,
        )

    if vcf_summary.get("caller") == "clairs":
        vcf_path = Path(str(vcf_summary.get("out_vcf")))
    else:
        vcf_path = Path(str(vcf_summary.get("out_snv_vcf")))

    if dry_run:
        logger.info("Dry-run: would assign reads using %s and %s", tumor_bam_p, vcf_path)
        return outdir_p / "report.html"

    outdir_p = ensure_outdir(outdir_p)
    alignment_dir = ensure_outdir(alignment_dir)
    vcf_dir = ensure_outdir(vcf_dir)
    readsep_dir = ensure_outdir(readsep_dir)

    # 3) Read separation
    if resume and (readsep_dir / "summary.json").exists():
        logger.info("Resume enabled: read separation summary already exists.")
        return outdir_p / "report.html"

    anchors, purity_used, anchor_stats = load_anchor_variants(
        str(vcf_path),
        sample=None,
        purity=purity,
        min_dp=10,
        require_pass=True,
    )
    anchor_index = build_anchor_index(anchors)

    tagged_bam_path: Optional[Path] = None
    if write_tagged_bam:
        tagged_bam_path = readsep_dir / "tagged.bam"

    with _step_log(outdir_p / "logs" / "assign.log"):
        run_summary = assign_bam(
            bam_path=str(tumor_bam_p),
            anchor_index_by_contig=anchor_index,
            prior_p_tumor=float(purity_used),
            outdir=readsep_dir,
            tumor_threshold=tumor_threshold,
            normal_threshold=normal_threshold,
            min_baseq=min_baseq,
            skip_duplicates=skip_duplicates,
            include_secondary=False,
            include_supplementary=False,
            write_tagged_bam=str(tagged_bam_path) if tagged_bam_path else None,
            split_bams=split_bams,
            assignments_tsv_gz=str(readsep_dir / "assignments.tsv.gz"),
            progress=True,
        )

    # Plots + report
    plots_dir = ensure_outdir(readsep_dir / "plots")

    plot_class_counts(class_counts=run_summary["counts"], out_png=plots_dir / "class_counts.png")
    plot_posterior_hist(
        bin_edges=run_summary["posterior_hist"]["bin_edges"],
        counts=run_summary["posterior_hist"]["counts"],
        out_png=plots_dir / "posterior_hist.png",
    )
    plot_evidence_hist(
        evidence_hist=run_summary["evidence_count_hist"],
        out_png=plots_dir / "evidence_hist.png",
    )
    plot_uncertainty_hist(
        posterior_hist=run_summary["posterior_hist"],
        out_png=plots_dir / "uncertainty_hist.png",
    )

    (readsep_dir / "anchor_stats.json").write_text(json.dumps(anchor_stats, indent=2), encoding="utf-8")

    readsep_report = render_report(
        outdir=readsep_dir,
        version=__version__,
        run=run_summary,
        anchor_stats=anchor_stats,
        vcf_path=str(vcf_path),
        sample="AUTO",
        plots={
            "class_counts": "plots/class_counts.png",
            "posterior_hist": "plots/posterior_hist.png",
            "evidence_hist": "plots/evidence_hist.png",
            "uncertainty_hist": "plots/uncertainty_hist.png",
        },
    )

    # 4) Top-level pipeline report + summary.json
    dt = time.time() - t0

    pipeline_summary = {
        "ref_fa": str(ref_fa_p),
        "tumor_bam": str(tumor_bam_p),
        "normal_bam": str(normal_bam_p) if normal_bam_p is not None else None,
        "alignment": aln_summary,
        "somatic_vcf": vcf_summary,
        "readsep": {
            "outdir": str(readsep_dir),
            "report_html": str(readsep_report),
            "purity_used": float(purity_used),
            "anchor_stats": anchor_stats,
            "counts": run_summary.get("counts", {}),
        },
        "runtime_seconds": float(dt),
    }

    pipeline_summary_path = outdir_p / "pipeline_summary.json"
    write_json(pipeline_summary_path, pipeline_summary)

    html = _PIPELINE_TEMPLATE.render(
        generated_at=_dt.datetime.now().isoformat(timespec="seconds"),
        version=__version__,
        alignment_dir=str(alignment_dir),
        vcf_dir=str(vcf_dir),
        readsep_dir=str(readsep_dir),
        readsep_report_rel=str(readsep_report.relative_to(outdir_p)),
        pipeline_summary=str(pipeline_summary_path.name),
        tumor_bam=str(tumor_bam_p),
        normal_bam=str(normal_bam_p) if normal_bam_p is not None else None,
        ref_fa=str(ref_fa_p),
        minimap2_preset=minimap2_preset,
        threads=int(threads),
        alignment_runtime=float(aln_summary.get("tumor", {}).get("runtime_seconds", 0.0)),
        caller=str(vcf_summary.get("caller")),
        engine=str(vcf_summary.get("engine")),
        platform=str(vcf_summary.get("platform")),
        vcf_path=str(vcf_path),
        vcf_runtime=float(vcf_summary.get("runtime_seconds", 0.0)),
    )

    report_path = outdir_p / "report.html"
    report_path.write_text(html, encoding="utf-8")

    logger.info("Pipeline complete. Report: %s", report_path)
    return report_path
