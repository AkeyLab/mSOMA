import pytest
from click.testing import CliRunner

from msoma.command_line import msoma
from msoma import __version__

runner = CliRunner()


def test_version():
    result = runner.invoke(msoma, ["--version"])
    assert result.exit_code == 0
    assert result.output == f"msoma, version {__version__}\n"


@pytest.mark.docker
def test_count():
    result = runner.invoke(
        msoma,
        [
            "count",
            "--fasta",
            "tests/data/test.fa",
            "--bed",
            "tests/data/test.bed",
            "--seq-length",
            "75",
            "--output",
            "test.count",
            "tests/data/test.bam",
        ],
    )
    assert result.exit_code == 0


@pytest.mark.docker
def test_count_mpileup_intermediate():
    result = runner.invoke(
        msoma,
        [
            "count",
            "--fasta",
            "tests/data/test.fa",
            "--bed",
            "tests/data/test.bed",
            "--seq-length",
            "75",
            "--output",
            "test.count",
            "--mpileup-intermediate",
            "test.mpileup",
            "tests/data/test.bam",
        ],
    )
    assert result.exit_code == 0


@pytest.mark.docker
def test_qname_whitelist():
    result = runner.invoke(
        msoma,
        [
            "qname-whitelist",
            "-o",
            "qnames.txt",
            "-q",
            "20",
            "-f",
            "2",
            "-F",
            "3844",
            "-I",
            "0",
            "-m",
            "0.1",
            "-s",
            "0.5",
            "tests/data/test.bam",
        ],
    )
    assert result.exit_code == 0


def test_fail_count_no_input_bam():
    result = runner.invoke(
        msoma,
        [
            "count",
            "--fasta",
            "tests/data/test.fa",
            "--seq-length",
            "75",
            "--output",
            "test.count",
        ],
    )
    assert result.exit_code != 0
    assert "Error: Missing argument 'INPUT_BAM'" in result.output
