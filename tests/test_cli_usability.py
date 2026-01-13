import subprocess
import sys
from pathlib import Path

import pysam

from gbmreadsep.toy_data import make_toy_data


def _run_cli(args: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, "-m", "gbmreadsep"] + args,
        check=False,
        capture_output=True,
        text=True,
    )


def _make_small_vcf(path: Path, contig: str) -> Path:
    header = pysam.VariantHeader()
    header.add_meta("fileformat", "VCFv4.2")
    header.add_sample("TUMOR")
    header.contigs.add(contig, length=200)
    header.formats.add("AF", number=1, type="Float", description="Alt allele fraction")
    header.formats.add("DP", number=1, type="Integer", description="Depth")

    vcf_path = path / "anchors.vcf"
    with pysam.VariantFile(str(vcf_path), "w", header=header) as vcf:
        rec = vcf.new_record(
            contig=contig,
            start=49,
            stop=50,
            alleles=("A", "G"),
            id=f"{contig}:50:A:G",
            qual=60,
            filter="PASS",
        )
        rec.samples[0]["AF"] = 0.3
        rec.samples[0]["DP"] = 50
        vcf.write(rec)

    vcf_gz = path / "anchors.vcf.gz"
    pysam.tabix_compress(str(vcf_path), str(vcf_gz), force=True)
    pysam.tabix_index(str(vcf_gz), preset="vcf", force=True)
    return vcf_gz


def test_quickstart_output() -> None:
    cp = _run_cli(["quickstart"])
    assert cp.returncode == 0
    assert "gbmreadsep assign" in cp.stdout
    assert "gbmreadsep nanopore endtoend" in cp.stdout


def test_wizard_non_interactive_command(tmp_path: Path) -> None:
    toy = make_toy_data(outdir=tmp_path / "toy")
    cp = _run_cli(
        [
            "wizard",
            "--mode",
            "assign",
            "--bam",
            toy["tumor_bam"],
            "--vcf",
            toy["anchors_vcf"],
            "--outdir",
            str(tmp_path / "out"),
            "--purity",
            "0.7",
            "--non-interactive",
        ]
    )
    assert cp.returncode == 0
    assert "gbmreadsep assign" in cp.stdout
    assert "--purity 0.7" in cp.stdout


def test_assign_dry_run_does_not_write_outputs(tmp_path: Path) -> None:
    toy = make_toy_data(outdir=tmp_path / "toy")
    outdir = tmp_path / "assign"
    cp = _run_cli(
        [
            "assign",
            "--bam",
            toy["tumor_bam"],
            "--vcf",
            toy["anchors_vcf"],
            "--outdir",
            str(outdir),
            "--dry-run",
        ]
    )
    assert cp.returncode == 0
    assert "Dry-run" in cp.stdout
    assert not (outdir / "summary.json").exists()


def test_make_toy_data_and_assign(tmp_path: Path) -> None:
    toy_dir = tmp_path / "toy"
    cp = _run_cli(["make-toy-data", "--outdir", str(toy_dir)])
    assert cp.returncode == 0

    outdir = tmp_path / "out"
    cp = _run_cli(
        [
            "assign",
            "--bam",
            str(toy_dir / "tumor.bam"),
            "--vcf",
            str(toy_dir / "anchors.vcf.gz"),
            "--outdir",
            str(outdir),
        ]
    )
    assert cp.returncode == 0
    assert (outdir / "report.html").exists()


def test_contig_mismatch_message(tmp_path: Path) -> None:
    toy = make_toy_data(outdir=tmp_path / "toy")
    vcf_dir = tmp_path / "vcf"
    vcf_dir.mkdir()
    vcf_gz = _make_small_vcf(vcf_dir, contig="1")

    cp = _run_cli(
        [
            "assign",
            "--bam",
            toy["tumor_bam"],
            "--vcf",
            str(vcf_gz),
            "--outdir",
            str(tmp_path / "out"),
            "--contig-style",
            "ensembl",
        ]
    )
    assert cp.returncode != 0
    assert "Contig mismatch" in cp.stderr
