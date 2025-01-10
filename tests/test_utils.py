from collections import Counter

import pytest

from msoma import utils as bb_utils
import gzip
import io


@pytest.mark.parametrize(
    "filename,mode",
    [
        ("-", ""),
        ("-", "rb"),
        ("-", "r+"),
        ("-", "wb"),
        ("-", "w+"),
    ],
)
def test_smart_open_invalid_mode(filename, mode):
    with pytest.raises(ValueError):
        with bb_utils.smart_open(filename=filename, mode=mode) as fh:
            pass


@pytest.mark.parametrize("filename", ["-", None])
def test_smart_open_read_stdin(monkeypatch, filename):
    test_lines = ["test", "test2"]
    monkeypatch.setattr("sys.stdin", io.StringIO("\n".join(test_lines)))
    with bb_utils.smart_open(filename=filename, mode="r") as fh:
        for line in fh:
            assert line.strip() == test_lines.pop(0)


@pytest.mark.parametrize("filename", ["-", None])
def test_smart_open_write_stdout(capfd, filename):
    test_lines = ["test", "test2"]
    with bb_utils.smart_open(filename=filename, mode="w") as fh:
        for line in test_lines:
            fh.write(line + "\n")

    out, err = capfd.readouterr()
    assert out == "\n".join(test_lines) + "\n"


def test_smart_open_read_flat(tmp_path):
    f_path = tmp_path / "test.txt"
    test_lines = ["test", "test2"]
    f_path.write_text("\n".join(test_lines))
    with bb_utils.smart_open(filename=str(f_path), mode="r") as fh:
        for line in fh:
            assert line.strip() == test_lines.pop(0)


def test_smart_open_write_flat(tmp_path):
    f_path = tmp_path / "test.txt"
    test_lines = ["test", "test2"]
    with bb_utils.smart_open(filename=str(f_path), mode="w") as fh:
        for line in test_lines:
            fh.write(line + "\n")

    assert f_path.read_text() == "\n".join(test_lines) + "\n"


def test_smart_open_read_gzip(tmp_path):
    f_path = tmp_path / "test.txt.gz"
    test_lines = ["test", "test2"]
    with gzip.open(filename=str(f_path), mode="wt") as fh:
        for line in test_lines:
            fh.write(line + "\n")

    with bb_utils.smart_open(filename=str(f_path), mode="r") as fh:
        for line in fh:
            assert line.strip() == test_lines.pop(0)


def test_smart_open_write(tmp_path):
    f_path = tmp_path / "test.txt.gz"
    test_lines = ["test", "test2"]
    with bb_utils.smart_open(filename=str(f_path), mode="w") as fh:
        for line in test_lines:
            fh.write(line + "\n")

    with gzip.open(filename=str(f_path), mode="rt") as fh:
        for line in fh:
            assert line.strip() == test_lines.pop(0)


@pytest.mark.parametrize(
    "locus, context",
    [
        (("chr1", 2), ("A", "C")),
        (("chr1", 3), ("A", "G")),
        (("chr2", 5), ("G", "A")),
    ],
)
def test_get_context(fasta, locus, context):
    chrom, pos = locus
    bef, aft = context
    BEF, REF, AFT = bb_utils.get_context(fasta, chrom, pos)
    assert BEF == bef
    assert AFT == aft


@pytest.mark.parametrize(
    "locus, context, error",
    [
        (("Q", 2), ("A", "T"), KeyError),  # chromosome does not exist
        (("1", 1), ("C", "T"), ValueError),  # position is out of bounds
    ],
)
def test_get_context_exceptions(fasta, locus, context, error):
    chrom, pos = locus
    bef, aft = context
    with pytest.raises(error):
        bb_utils.get_context(fasta, chrom, pos) == (bef, aft)


@pytest.mark.parametrize(
    "cigarstring,exp_count",
    [
        ("10M", {"M": 10}),
        ("10M1D", {"M": 10, "D": 1}),
        ("76S130M", {"S": 76, "M": 130}),
        ("10M1D1M", {"D": 1, "M": 11}),
        (
            "10M1D1M1M",
            {"D": 1, "M": 12},
        ),  # technically invalid cigarstring, but count works
    ],
)
def test_cigarstring_counts(cigarstring, exp_count):
    counts = bb_utils.cigarstring_counts(cigarstring)
    assert counts == exp_count


@pytest.mark.parametrize(
    "cigarstring",
    [
        "10Q",  # Q is not a valid cigar operation
        "10M1D1M1",  # missing operation character
    ],
)
def test_cigarstring_counts_invalid(cigarstring):
    with pytest.raises(ValueError):
        bb_utils.cigarstring_counts(cigarstring)
