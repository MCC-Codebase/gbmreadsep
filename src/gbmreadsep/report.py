from __future__ import annotations

import datetime as _dt
import logging
from pathlib import Path
from typing import Any, Dict

from jinja2 import Template

logger = logging.getLogger(__name__)


_REPORT_TEMPLATE = Template(
    """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>GBMReadSep Report</title>
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
    img { max-width: 100%; height: auto; border: 1px solid #eee; border-radius: 6px; }
  </style>
</head>
<body>

<h1>GBMReadSep Report</h1>
<p class="small">Generated: {{ generated_at }}</p>

<h2>Run summary</h2>
<div class="grid">
  <div class="card">
    <h3>Inputs</h3>
    <table>
      <tr><th>BAM</th><td><code>{{ bam_path }}</code></td></tr>
      <tr><th>VCF</th><td><code>{{ vcf_path }}</code></td></tr>
      <tr><th>Sample</th><td><code>{{ sample }}</code></td></tr>
    </table>
  </div>
  <div class="card">
    <h3>Model</h3>
    <table>
      <tr><th>Prior tumor fraction (purity)</th><td>{{ prior_p_tumor }}</td></tr>
      <tr><th>Tumor threshold</th><td>{{ tumor_threshold }}</td></tr>
      <tr><th>Normal threshold</th><td>{{ normal_threshold }}</td></tr>
      <tr><th>Min baseQ</th><td>{{ min_baseq }}</td></tr>
    </table>
  </div>
</div>

<h2>Anchors</h2>
<table>
  <tr><th>Anchors kept</th><td>{{ anchor_stats.anchors_kept }}</td></tr>
  <tr><th>VCF records total</th><td>{{ anchor_stats.records_total }}</td></tr>
  <tr><th>SNV records</th><td>{{ anchor_stats.records_snv }}</td></tr>
  <tr><th>Records with VAF</th><td>{{ anchor_stats.records_with_vaf }}</td></tr>
  <tr><th>Skipped non-SNV</th><td>{{ anchor_stats.anchors_skipped_non_snv }}</td></tr>
  <tr><th>Skipped DP</th><td>{{ anchor_stats.anchors_skipped_dp }}</td></tr>
  <tr><th>Skipped FILTER</th><td>{{ anchor_stats.anchors_skipped_filter }}</td></tr>
</table>

<h2>Read assignments</h2>
<table>
  <tr><th>Total reads seen</th><td>{{ counts.reads_total }}</td></tr>
  <tr><th>Unmapped skipped</th><td>{{ counts.reads_unmapped }}</td></tr>
  <tr><th>Duplicates skipped</th><td>{{ counts.reads_skipped_duplicates }}</td></tr>
  <tr><th>Secondary skipped</th><td>{{ counts.reads_skipped_secondary }}</td></tr>
  <tr><th>Supplementary skipped</th><td>{{ counts.reads_skipped_supplementary }}</td></tr>
  <tr><th>Reads with evidence</th><td>{{ counts.reads_with_evidence }}</td></tr>
  <tr><th>Tumor (T)</th><td>{{ counts.class_T }}</td></tr>
  <tr><th>TME/Normal (N)</th><td>{{ counts.class_N }}</td></tr>
  <tr><th>Uncertain (U)</th><td>{{ counts.class_U }}</td></tr>
</table>

<h2>Plots</h2>

<div class="grid">
  <div class="card">
    <h3>Assignment counts</h3>
    <img src="{{ plots.class_counts }}" alt="class counts">
  </div>
  <div class="card">
    <h3>Posterior distribution</h3>
    <img src="{{ plots.posterior_hist }}" alt="posterior histogram">
  </div>
</div>

<div class="grid" style="margin-top:16px;">
  <div class="card">
    <h3>Evidence per read</h3>
    <img src="{{ plots.evidence_hist }}" alt="evidence histogram">
  </div>
  <div class="card">
    <h3>Uncertainty</h3>
    <img src="{{ plots.uncertainty_hist }}" alt="uncertainty histogram">
  </div>
</div>

<h2>Outputs</h2>
<ul>
  <li><code>{{ assignments_tsv_gz }}</code> (per-read assignments)</li>
  {% if write_tagged_bam %}
  <li><code>{{ write_tagged_bam }}</code> (tagged BAM)</li>
  {% endif %}
  {% if split_bams %}
  <li><code>tumor.bam</code>, <code>tme.bam</code>, <code>uncertain.bam</code> (split BAMs)</li>
  {% endif %}
  <li><code>summary.json</code> (machine-readable summary)</li>
</ul>

<h2>Interpretation notes</h2>
<ul>
  <li>Reads not overlapping somatic anchors are intrinsically ambiguous and default to the prior.</li>
  <li>"Normal" here means "non-tumor" (TME/normal contamination) and is not immune-specific.</li>
  <li>Purity estimation is heuristic; provide <code>--purity</code> for best control.</li>
</ul>

<hr>
<p class="small">GBMReadSep {{ version }}</p>
</body>
</html>"""
)


def render_report(
    *,
    outdir: str | Path,
    version: str,
    run: Dict[str, Any],
    anchor_stats: Dict[str, Any],
    vcf_path: str,
    sample: str,
    plots: Dict[str, str],
) -> Path:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Backward compatibility: if uncertainty plot missing, re-use posterior.
    plots = dict(plots)
    plots.setdefault("uncertainty_hist", plots.get("posterior_hist", ""))

    html = _REPORT_TEMPLATE.render(
        generated_at=_dt.datetime.now().isoformat(timespec="seconds"),
        version=version,
        bam_path=run.get("bam_path"),
        vcf_path=vcf_path,
        sample=sample,
        prior_p_tumor=run.get("prior_p_tumor"),
        tumor_threshold=run.get("tumor_threshold"),
        normal_threshold=run.get("normal_threshold"),
        min_baseq=run.get("min_baseq"),
        assignments_tsv_gz=run.get("assignments_tsv_gz"),
        write_tagged_bam=run.get("write_tagged_bam"),
        split_bams=run.get("split_bams"),
        counts=run.get("counts", {}),
        anchor_stats=anchor_stats,
        plots=plots,
    )

    out_path = outdir / "report.html"
    out_path.write_text(html, encoding="utf-8")
    return out_path
