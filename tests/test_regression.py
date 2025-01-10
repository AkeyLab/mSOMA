import pytest
from click.testing import CliRunner
import pysam
from msoma.command_line import msoma
import random
from pathlib import Path

runner = CliRunner()


@pytest.mark.docker
def test_filter_passing_and_non_primary_alignment(tmp_path):
    """
    This is a simplified regression test for the following bug/feature:
        A paired-end alignment has an R1 and R2 proper alignment on Chr1
        but R1 ALSO has a non-primary alignment on Chr3

    The flags correspond to
        99: read paired, read mapped in proper pair, mate reverse strand, first in pair (R1)
        147: read paired, read mapped in proper pair, read reverse strand, second in pair (R2)
        329: read paired, mate unmapped, first in pair, not primary alignment (R1)

    In this situation we want to count the R1 and R2 alignments on Chr1
    but not the R1 alignment on Chr3
    """
    SAM_str = "\n".join(
        (
            "@HD	VN:1.6	SO:coordinate",
            "@SQ	SN:1	LN:100",
            "@SQ	SN:2	LN:100",
            "@SQ	SN:3	LN:100",
            "QNAME1\t99\t1\t10\t60\t10M\t*\t0\t0\tAAAAAAAAAA\t==========\tNM:i:0",
            "QNAME1\t147\t1\t50\t60\t10M\t*\t0\t0\tAAAAAAAAAA\t==========\tNM:i:0",
            "QNAME1\t329\t3\t30\t0\t10M\t*\t0\t0\t*\t*\tNM:i:0",
        )
    )
    sam_path = tmp_path / "test.sam"
    sam_path.write_text(SAM_str)
    bam_path = tmp_path / "test.bam"
    pysam.tabix_compress(str(sam_path), str(bam_path))
    pysam.index(str(bam_path))

    fasta_str = "\n".join(
        (
            ">1",
            "".join(random.choices(["A", "C", "G", "T"], k=100)),
            ">2",
            "".join(random.choices(["A", "C", "G", "T"], k=100)),
            ">3",
            "".join(random.choices(["A", "C", "G", "T"], k=100)),
        )
    )
    fasta_path = tmp_path / "test.fasta"
    fasta_path.write_text(fasta_str)

    bed_str = "\n".join(
        (
            "1\t0\t100",
            "2\t0\t100",
            "3\t0\t100",
        )
    )
    bed_path = tmp_path / "test.bed"
    bed_path.write_text(bed_str)

    count_path = str(tmp_path / "test.count.gz")
    mpileup_path = str(tmp_path / "test.mpileup")

    result = runner.invoke(
        msoma,
        [
            "count",
            "--fasta",
            str(fasta_path),
            "--bed",
            str(bed_path),
            "--seq-length",
            "10",
            "--ntrim",
            "0",
            "--mpileup-intermediate",
            mpileup_path,
            "--output",
            count_path,
            str(bam_path),
        ],
    )

    assert result.exit_code == 0

    # QNAME1 should be present in the mpileup intermediate file
    # this is QNAME coming from the R1 alignment on Chr1
    with open(mpileup_path, "r") as f_in:
        assert "QNAME1" in f_in.read()
