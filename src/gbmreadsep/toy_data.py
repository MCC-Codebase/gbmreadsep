from __future__ import annotations

import random
from pathlib import Path
from typing import Dict, List, Tuple

import pysam

from .utils import ensure_outdir, write_json


def _write_fasta(path: Path, contig: str, seq: str) -> None:
    lines = [f">{contig}"]
    for i in range(0, len(seq), 60):
        lines.append(seq[i : i + 60])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _mutate_base(base: str) -> str:
    for alt in ["A", "C", "G", "T"]:
        if alt != base:
            return alt
    return "A"


def _make_read(
    name: str,
    contig: str,
    start0: int,
    seq: str,
    mapq: int = 60,
) -> pysam.AlignedSegment:
    a = pysam.AlignedSegment()
    a.query_name = name
    a.query_sequence = seq
    a.flag = 0
    a.reference_id = 0
    a.reference_start = start0
    a.mapping_quality = mapq
    a.cigartuples = [(0, len(seq))]
    a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
    return a


def make_toy_data(*, outdir: str | Path) -> Dict[str, str]:
    """Create a tiny reference, BAM, and VCF suitable for quick demos/tests.

    The outputs include:
    - toy_ref.fa (+ .fai)
    - tumor.bam (+ .bai)
    - anchors.vcf.gz (+ .tbi)

    Returns
    -------
    dict
        Paths to the generated files.
    """
    outdir_p = ensure_outdir(outdir)

    contig = "chr1"
    ref_seq = ("ACGT" * 50)[:200]
    ref_fa = outdir_p / "toy_ref.fa"
    _write_fasta(ref_fa, contig, ref_seq)
    pysam.faidx(str(ref_fa))

    anchor_positions = [50, 120]  # 0-based positions
    anchors: List[Tuple[int, str, str]] = []
    for pos0 in anchor_positions:
        ref_base = ref_seq[pos0]
        alt_base = _mutate_base(ref_base)
        anchors.append((pos0, ref_base, alt_base))

    # Build BAM
    bam_path = outdir_p / "tumor.bam"
    header = {
        "HD": {"VN": "1.6"},
        "SQ": [{"SN": contig, "LN": len(ref_seq)}],
    }

    reads: List[pysam.AlignedSegment] = []
    rng = random.Random(7)

    for i in range(10):
        start0 = 30 + i
        seq = list(ref_seq[start0 : start0 + 50])
        # mutate anchor 0 for half the reads
        if i % 2 == 0:
            rel = anchor_positions[0] - start0
            if 0 <= rel < len(seq):
                seq[rel] = anchors[0][2]
        reads.append(_make_read(f"r1_{i}", contig, start0, "".join(seq)))

    for i in range(10):
        start0 = 100 + i
        seq = list(ref_seq[start0 : start0 + 50])
        if i % 2 == 1:
            rel = anchor_positions[1] - start0
            if 0 <= rel < len(seq):
                seq[rel] = anchors[1][2]
        # shuffle slightly to add minor variability
        if rng.random() < 0.2:
            seq[-1] = _mutate_base(seq[-1])
        reads.append(_make_read(f"r2_{i}", contig, start0, "".join(seq)))

    reads.sort(key=lambda r: r.reference_start)

    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
        for r in reads:
            bam.write(r)

    pysam.index(str(bam_path))

    # Build VCF
    vcf_path = outdir_p / "anchors.vcf"
    header = pysam.VariantHeader()
    header.add_meta("fileformat", "VCFv4.2")
    header.add_sample("TUMOR")
    header.contigs.add(contig, length=len(ref_seq))
    header.formats.add("AF", number=1, type="Float", description="Alt allele fraction")
    header.formats.add("DP", number=1, type="Integer", description="Depth")

    with pysam.VariantFile(str(vcf_path), "w", header=header) as vcf:
        for pos0, ref_base, alt_base in anchors:
            rec = vcf.new_record(
                contig=contig,
                start=pos0,
                stop=pos0 + 1,
                alleles=(ref_base, alt_base),
                id=f"{contig}:{pos0+1}:{ref_base}:{alt_base}",
                qual=60,
                filter="PASS",
            )
            rec.samples[0]["AF"] = 0.3
            rec.samples[0]["DP"] = 50
            vcf.write(rec)

    vcf_gz = outdir_p / "anchors.vcf.gz"
    pysam.tabix_compress(str(vcf_path), str(vcf_gz), force=True)
    pysam.tabix_index(str(vcf_gz), preset="vcf", force=True)

    summary = {
        "ref_fa": str(ref_fa),
        "tumor_bam": str(bam_path),
        "anchors_vcf": str(vcf_gz),
        "outdir": str(outdir_p),
    }

    write_json(outdir_p / "toy_summary.json", summary)
    return summary
