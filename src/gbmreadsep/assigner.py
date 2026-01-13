from __future__ import annotations

import bisect
import logging
import time
from dataclasses import asdict
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Tuple

import numpy as np
import pysam
from tqdm import tqdm

from .anchors import AnchorIndex
from .models import AnchorVariant, ReadAssignment
from .utils import clamp, ensure_outdir, logit, open_textmaybe_gzip, phred_to_error_prob, sigmoid, write_json

logger = logging.getLogger(__name__)


def _base_prob(obs_base: str, true_base: str, e: float) -> float:
    """Probability of observing obs_base given true_base under a simple substitution model."""
    if obs_base == true_base:
        return 1.0 - e
    # distribute sequencing error uniformly among the other 3 bases
    return e / 3.0


def _effective_error_prob(baseq: int, mapq: int, *, cap: float = 0.25) -> float:
    """Combine base quality and mapping quality into an effective substitution error.

    This is a heuristic: mapping errors are not strictly substitution errors, but low MAPQ increases
    the risk that the observed base is uninformative for the intended locus.

    cap:
        Upper cap on error to avoid numerical pathologies on very low-quality reads.
    """
    e_seq = phred_to_error_prob(baseq)
    e_map = phred_to_error_prob(mapq)
    # independence-style union
    e = 1.0 - (1.0 - e_seq) * (1.0 - e_map)
    return clamp(e, 1e-6, cap)


def llr_for_observation(
    anchor: AnchorVariant,
    obs_base: str,
    baseq: int,
    mapq: int,
) -> float:
    """Log-likelihood ratio log P(obs | tumor) / P(obs | normal) for a single SNV observation."""
    obs = obs_base.upper()
    if obs not in {"A", "C", "G", "T"}:
        # Non-ACGT bases carry minimal information; treat as uninformative.
        return 0.0

    e = _effective_error_prob(baseq, mapq)

    p_obs_given_ref = _base_prob(obs, anchor.ref, e)
    p_obs_given_alt = _base_prob(obs, anchor.alt, e)

    # Tumor: mixture of ALT and REF alleles per tumor read with f_tumor.
    f = clamp(anchor.f_tumor, 0.0, 1.0)
    p_tumor = f * p_obs_given_alt + (1.0 - f) * p_obs_given_ref

    # Normal: assumed REF at somatic anchor
    p_normal = p_obs_given_ref

    # Numerical guards
    p_tumor = max(p_tumor, 1e-12)
    p_normal = max(p_normal, 1e-12)
    return float(np.log(p_tumor / p_normal))


def extract_bases_at_positions(
    read: pysam.AlignedSegment, positions0: List[int]
) -> Dict[int, Tuple[str, int]]:
    """Extract (base, baseq) at a sorted list of reference positions (0-based) for one read.

    This function walks the CIGAR once and extracts bases for candidate positions without iterating
    over all aligned pairs.

    Returns a mapping {pos0: (base, base_quality)} only for positions that are aligned as a base
    in the read (i.e. not deletions / ref skips at that locus).
    """
    if read.is_unmapped or read.cigartuples is None or len(positions0) == 0:
        return {}

    seq = read.query_sequence
    if seq is None:
        return {}
    quals = read.query_qualities  # can be None

    out: Dict[int, Tuple[str, int]] = {}

    pos_idx = 0
    ref_pos = read.reference_start
    query_pos = 0

    # Fast-forward if positions are before the read
    while pos_idx < len(positions0) and positions0[pos_idx] < ref_pos:
        pos_idx += 1

    for op, length in read.cigartuples:
        if pos_idx >= len(positions0):
            break

        if op in (0, 7, 8):  # M, =, X: consumes query and ref
            ref_end = ref_pos + length

            # process all candidate positions within this match block
            while pos_idx < len(positions0) and positions0[pos_idx] < ref_end:
                p0 = positions0[pos_idx]
                if p0 >= ref_pos:
                    qpos = query_pos + (p0 - ref_pos)
                    if 0 <= qpos < len(seq):
                        base = seq[qpos]
                        bq = int(quals[qpos]) if quals is not None else 0
                        out[p0] = (base, bq)
                pos_idx += 1

            ref_pos = ref_end
            query_pos += length

        elif op == 1:  # I: consumes query only
            query_pos += length
        elif op in (2, 3):  # D, N: consumes ref only
            ref_pos += length
        elif op == 4:  # S: consumes query only
            query_pos += length
        elif op in (5, 6):  # H, P: consumes neither
            continue
        else:
            # Unknown/rare CIGAR op, treat conservatively
            continue

    return out


def iter_candidate_anchors_for_read(
    anchor_index: AnchorIndex,
    start0: int,
    end0: int,
) -> Tuple[List[int], List[AnchorVariant]]:
    """Return the anchor positions and objects within [start0, end0)."""
    left = bisect.bisect_left(anchor_index.positions, start0)
    right = bisect.bisect_left(anchor_index.positions, end0)
    return anchor_index.positions[left:right], anchor_index.anchors[left:right]


def assign_read(
    read: pysam.AlignedSegment,
    *,
    anchor_index_by_contig: Dict[str, AnchorIndex],
    prior_p_tumor: float,
    tumor_threshold: float,
    normal_threshold: float,
    min_baseq: int,
) -> ReadAssignment:
    """Assign a single read and return posterior + class."""
    chrom = read.reference_name if read.reference_name is not None else "*"
    start0 = int(read.reference_start) if read.reference_start is not None else -1
    end0 = int(read.reference_end) if read.reference_end is not None else -1

    # Prior
    prior = clamp(prior_p_tumor, 1e-6, 1.0 - 1e-6)
    llr_sum = 0.0
    evidence_count = 0

    if chrom in anchor_index_by_contig and start0 >= 0 and end0 > start0:
        idx = anchor_index_by_contig[chrom]
        pos_list, anchors = iter_candidate_anchors_for_read(idx, start0, end0)
        if pos_list:
            bases = extract_bases_at_positions(read, pos_list)
            mapq = int(read.mapping_quality)
            for pos0, anchor in zip(pos_list, anchors):
                if pos0 not in bases:
                    continue
                base, bq = bases[pos0]
                if bq < min_baseq:
                    continue
                llr_sum += llr_for_observation(anchor, base, bq, mapq)
                evidence_count += 1

    post = sigmoid(logit(prior) + llr_sum)

    if post >= tumor_threshold:
        cls = "T"
    elif post <= normal_threshold:
        cls = "N"
    else:
        cls = "U"

    return ReadAssignment(
        qname=str(read.query_name),
        chrom=chrom,
        start0=start0,
        end0=end0,
        is_read1=bool(read.is_read1),
        is_read2=bool(read.is_read2),
        mapq=int(read.mapping_quality),
        evidence_count=evidence_count,
        llr_sum=float(llr_sum),
        p_tumor=float(post),
        cls=cls,
    )


def assign_bam(
    *,
    bam_path: str,
    anchor_index_by_contig: Dict[str, AnchorIndex],
    prior_p_tumor: float,
    outdir: str | Path,
    tumor_threshold: float = 0.9,
    normal_threshold: float = 0.1,
    min_baseq: int = 20,
    skip_duplicates: bool = True,
    include_secondary: bool = False,
    include_supplementary: bool = False,
    write_tagged_bam: Optional[str] = None,
    split_bams: bool = False,
    assignments_tsv_gz: Optional[str] = None,
    progress: bool = True,
) -> Dict[str, object]:
    """Main workhorse: iterate BAM, assign reads, write outputs, and return summary dict."""
    t0 = time.time()
    outdir_path = ensure_outdir(outdir)

    if tumor_threshold <= normal_threshold:
        raise ValueError("tumor_threshold must be > normal_threshold")

    bam = pysam.AlignmentFile(bam_path, "rb")

    # Outputs
    tagged_bam_fh: Optional[pysam.AlignmentFile] = None
    tumor_bam_fh: Optional[pysam.AlignmentFile] = None
    tme_bam_fh: Optional[pysam.AlignmentFile] = None
    uncertain_bam_fh: Optional[pysam.AlignmentFile] = None

    if write_tagged_bam is not None:
        tagged_bam_fh = pysam.AlignmentFile(str(write_tagged_bam), "wb", template=bam)

    if split_bams:
        tumor_bam_fh = pysam.AlignmentFile(str(outdir_path / "tumor.bam"), "wb", template=bam)
        tme_bam_fh = pysam.AlignmentFile(str(outdir_path / "tme.bam"), "wb", template=bam)
        uncertain_bam_fh = pysam.AlignmentFile(
            str(outdir_path / "uncertain.bam"), "wb", template=bam
        )

    if assignments_tsv_gz is None:
        assignments_tsv_gz = str(outdir_path / "assignments.tsv.gz")

    tsv_fh = open_textmaybe_gzip(assignments_tsv_gz, "wt")
    tsv_fh.write(
        "\t".join(
            [
                "qname",
                "chrom",
                "start0",
                "end0",
                "is_read1",
                "is_read2",
                "mapq",
                "evidence_count",
                "llr_sum",
                "p_tumor",
                "class",
            ]
        )
        + "\n"
    )

    # Streaming histograms
    posterior_bins = np.linspace(0.0, 1.0, 101)
    posterior_counts = np.zeros(len(posterior_bins) - 1, dtype=np.int64)

    evidence_count_hist: Dict[int, int] = {}

    counts = {
        "reads_total": 0,
        "reads_written": 0,
        "reads_unmapped": 0,
        "reads_skipped_secondary": 0,
        "reads_skipped_supplementary": 0,
        "reads_skipped_duplicates": 0,
        "reads_with_evidence": 0,
        "class_T": 0,
        "class_N": 0,
        "class_U": 0,
    }

    # Iterate reads (unfiltered fetch uses BAM order)
    it: Iterable[pysam.AlignedSegment] = bam.fetch(until_eof=True)
    if progress:
        it = tqdm(it, unit="read", desc="Assigning reads")

    for read in it:
        counts["reads_total"] += 1

        if read.is_unmapped:
            counts["reads_unmapped"] += 1
            continue
        if read.is_secondary and not include_secondary:
            counts["reads_skipped_secondary"] += 1
            continue
        if read.is_supplementary and not include_supplementary:
            counts["reads_skipped_supplementary"] += 1
            continue
        if skip_duplicates and read.is_duplicate:
            counts["reads_skipped_duplicates"] += 1
            continue

        res = assign_read(
            read,
            anchor_index_by_contig=anchor_index_by_contig,
            prior_p_tumor=prior_p_tumor,
            tumor_threshold=tumor_threshold,
            normal_threshold=normal_threshold,
            min_baseq=min_baseq,
        )

        if res.evidence_count > 0:
            counts["reads_with_evidence"] += 1

        posterior_counts += np.histogram([res.p_tumor], bins=posterior_bins)[0]
        evidence_count_hist[res.evidence_count] = evidence_count_hist.get(res.evidence_count, 0) + 1

        if res.cls == "T":
            counts["class_T"] += 1
        elif res.cls == "N":
            counts["class_N"] += 1
        else:
            counts["class_U"] += 1

        # TSV
        tsv_fh.write(
            f"{res.qname}\t{res.chrom}\t{res.start0}\t{res.end0}\t"
            f"{int(res.is_read1)}\t{int(res.is_read2)}\t{res.mapq}\t"
            f"{res.evidence_count}\t{res.llr_sum:.6f}\t{res.p_tumor:.6f}\t{res.cls}\n"
        )

        # BAM tagging / splitting
        if tagged_bam_fh is not None or split_bams:
            # add tags to the read object (will be written to whichever outputs are enabled)
            read.set_tag("TP", float(res.p_tumor), value_type="f")
            read.set_tag("TC", res.cls, value_type="Z")
            read.set_tag("TE", int(res.evidence_count), value_type="i")
            read.set_tag("TL", float(res.llr_sum), value_type="f")

        if tagged_bam_fh is not None:
            tagged_bam_fh.write(read)
            counts["reads_written"] += 1

        if split_bams:
            assert tumor_bam_fh is not None and tme_bam_fh is not None and uncertain_bam_fh is not None
            if res.cls == "T":
                tumor_bam_fh.write(read)
            elif res.cls == "N":
                tme_bam_fh.write(read)
            else:
                uncertain_bam_fh.write(read)

    tsv_fh.close()

    # Close BAM outputs
    bam.close()
    if tagged_bam_fh is not None:
        tagged_bam_fh.close()
    if split_bams:
        assert tumor_bam_fh is not None and tme_bam_fh is not None and uncertain_bam_fh is not None
        tumor_bam_fh.close()
        tme_bam_fh.close()
        uncertain_bam_fh.close()

    dt = time.time() - t0

    summary = {
        "bam_path": bam_path,
        "prior_p_tumor": float(prior_p_tumor),
        "tumor_threshold": float(tumor_threshold),
        "normal_threshold": float(normal_threshold),
        "min_baseq": int(min_baseq),
        "skip_duplicates": bool(skip_duplicates),
        "include_secondary": bool(include_secondary),
        "include_supplementary": bool(include_supplementary),
        "assignments_tsv_gz": str(assignments_tsv_gz),
        "write_tagged_bam": str(write_tagged_bam) if write_tagged_bam is not None else None,
        "split_bams": bool(split_bams),
        "counts": counts,
        "posterior_hist": {
            "bin_edges": posterior_bins.tolist(),
            "counts": posterior_counts.tolist(),
        },
        "evidence_count_hist": evidence_count_hist,
        "runtime_seconds": float(dt),
    }

    write_json(outdir_path / "summary.json", summary)
    return summary
