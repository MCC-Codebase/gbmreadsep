import pysam

from gbmreadsep.assigner import extract_bases_at_positions
from gbmreadsep.models import AnchorVariant
from gbmreadsep.assigner import llr_for_observation


def make_read(seq: str, start: int = 100) -> pysam.AlignedSegment:
    a = pysam.AlignedSegment()
    a.query_name = "r1"
    a.query_sequence = seq
    a.flag = 0
    a.reference_start = start
    a.mapping_quality = 60
    a.cigartuples = [(0, len(seq))]  # M
    # high qualities
    a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
    return a


def test_extract_bases_simple_match():
    read = make_read("ACGTACGTAA", start=100)
    pos = [100, 103, 109]  # A, T, A
    out = extract_bases_at_positions(read, pos)
    assert out[100][0] == "A"
    assert out[103][0] == "T"
    assert out[109][0] == "A"


def test_llr_signs():
    anchor = AnchorVariant(
        chrom="1",
        pos0=100,
        ref="A",
        alt="G",
        vaf=0.3,
        dp=100,
        f_tumor=0.5,
        record_id="x",
    )
    llr_alt = llr_for_observation(anchor, "G", baseq=40, mapq=60)
    llr_ref = llr_for_observation(anchor, "A", baseq=40, mapq=60)
    assert llr_alt > 0
    assert llr_ref < 0
