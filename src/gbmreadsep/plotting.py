from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)


def plot_posterior_hist(
    *,
    bin_edges: List[float],
    counts: List[int],
    out_png: str | Path,
    title: str = "Posterior P(tumor) distribution",
) -> None:
    out_png = Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)

    widths = [bin_edges[i + 1] - bin_edges[i] for i in range(len(counts))]
    centers = [bin_edges[i] + widths[i] / 2.0 for i in range(len(counts))]

    plt.figure()
    plt.bar(centers, counts, width=widths, align="center")
    plt.xlabel("Posterior P(tumor)")
    plt.ylabel("Read count")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=160)
    plt.close()


def plot_class_counts(
    *,
    class_counts: Dict[str, int],
    out_png: str | Path,
    title: str = "Read assignments",
) -> None:
    out_png = Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)

    labels = ["Tumor (T)", "TME/Normal (N)", "Uncertain (U)"]
    values = [
        int(class_counts.get("class_T", 0)),
        int(class_counts.get("class_N", 0)),
        int(class_counts.get("class_U", 0)),
    ]

    plt.figure()
    plt.bar(labels, values)
    plt.ylabel("Read count")
    plt.title(title)
    plt.xticks(rotation=15, ha="right")
    plt.tight_layout()
    plt.savefig(out_png, dpi=160)
    plt.close()


def plot_evidence_hist(
    *,
    evidence_hist: Dict[int, int],
    out_png: str | Path,
    title: str = "Evidence count per read",
    max_bin: int = 20,
) -> None:
    out_png = Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)

    # Collapse tail into max_bin+
    xs = list(range(0, max_bin + 1))
    ys = [0] * len(xs)
    tail = 0
    for k, v in evidence_hist.items():
        if k <= max_bin:
            ys[k] += int(v)
        else:
            tail += int(v)

    xticklabels = [str(x) for x in range(0, max_bin + 1)]
    if tail > 0:
        xs.append(max_bin + 1)
        ys.append(tail)
        xticklabels.append(f"{max_bin + 1}+")

    plt.figure()
    plt.bar(range(len(xs)), ys)
    plt.xlabel("Number of anchor observations in read")
    plt.ylabel("Read count")
    plt.title(title)
    plt.xticks(range(len(xs)), xticklabels, rotation=0)
    plt.tight_layout()
    plt.savefig(out_png, dpi=160)
    plt.close()


def plot_uncertainty_hist(
    *,
    posterior_hist: Dict[str, object],
    out_png: str | Path,
    title: str = "Uncertainty distribution",
    nbins: int = 25,
) -> None:
    """Plot an uncertainty histogram from a posterior histogram.

    We define uncertainty as ``u = min(p, 1-p)``, which ranges from 0 (very certain)
    to 0.5 (maximally uncertain).

    Parameters
    ----------
    posterior_hist:
        Dict with keys ``bin_edges`` and ``counts`` (as produced by assign_bam).
    nbins:
        Number of uncertainty bins between [0, 0.5].

    Notes
    -----
    We only have a posterior histogram (not per-read posteriors) at this stage, so
    the plot is an approximation based on bin centers.
    """
    out_png = Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)

    bin_edges = list(map(float, posterior_hist.get("bin_edges", [])))
    counts = list(map(int, posterior_hist.get("counts", [])))
    if len(bin_edges) != len(counts) + 1:
        raise ValueError("posterior_hist must contain bin_edges of length len(counts)+1")

    # Map posterior bins -> uncertainty bins.
    u_edges = [i * 0.5 / nbins for i in range(nbins + 1)]
    u_counts = [0] * nbins

    for i, c in enumerate(counts):
        left = bin_edges[i]
        right = bin_edges[i + 1]
        p = 0.5 * (left + right)
        u = min(p, 1.0 - p)

        # Find uncertainty bin
        j = min(int(u / 0.5 * nbins), nbins - 1)
        if j < 0:
            j = 0
        u_counts[j] += int(c)

    widths = [u_edges[i + 1] - u_edges[i] for i in range(nbins)]
    centers = [u_edges[i] + widths[i] / 2.0 for i in range(nbins)]

    plt.figure()
    plt.bar(centers, u_counts, width=widths, align="center")
    plt.xlabel("Uncertainty u = min(p, 1-p)")
    plt.ylabel("Read count")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=160)
    plt.close()
