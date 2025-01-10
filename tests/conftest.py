import gzip
import io

import pysam
import pytest

from msoma import utils


def pytest_addoption(parser):
    parser.addoption(
        "--docker",
        action="store_true",
        default=False,
        help="run tests requiring docker container",
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "docker: mark test as requiring docker to run")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--docker"):
        # --docker given in cli: do not skip docker requiring tests
        return
    skip_docker = pytest.mark.skip(reason="need --docker option to run")
    for item in items:
        if "docker" in item.keywords:
            item.add_marker(skip_docker)


@pytest.fixture(scope="session")
def fasta_file(tmp_path_factory):
    fasta_data = "\n".join(
        (">chr1", "AACGTACATAGTGTATTGA", ">chr2", "ATCGGATTTTATTATTATG")
    )
    p = tmp_path_factory.mktemp("data") / "test.fasta"
    p.write_text(fasta_data + "\n")
    return str(p)


@pytest.fixture(scope="function")
def fasta(fasta_file):
    return pysam.FastaFile(fasta_file)


@pytest.fixture(scope="session")
def bam_file(tmp_path_factory):
    bam_out = tmp_path_factory.mktemp("data") / "test.bam"
    bam_out = str(bam_out)

    header = {
        "HD": {"VN": "1.0"},
        "SQ": [{"LN": 1575, "SN": "chr1"}, {"LN": 1584, "SN": "chr2"}],
    }

    reads = [
        utils._test_pysam_alignment(ref_id=0, mapq=10, flag=0),  # chr1 fwd
        utils._test_pysam_alignment(ref_id=1, mapq=20, flag=16),  # chr2 rev
        utils._test_pysam_alignment(
            ref_id=0, mapq=30, flag=81
        ),  # chr1 paired, rev, 1st in pair
        utils._test_pysam_alignment(
            ref_id=0, mapq=40, flag=145
        ),  # chr1 paired, rev, 2nd in pair
        utils._test_pysam_alignment(
            ref_id=1, mapq=50, flag=129
        ),  # chr2 paired, fwd, 2nd in pair
        utils._test_pysam_alignment(
            ref_id=0, mapq=255, flag=65
        ),  # chr1 paired, fwd, 1st in pair
        utils._test_pysam_alignment(
            ref_id=0, mapq=255, flag=321
        ),  # chr1 paired, fwd, 1st in pair, not primary align
    ]
    with pysam.AlignmentFile(bam_out, "wb", header=header) as outf:
        for r in reads:
            outf.write(r)

    pysam.sort("-o", bam_out, bam_out)
    pysam.index(bam_out)

    return bam_out


@pytest.fixture(scope="function")
def bam(bam_file):
    bam = pysam.AlignmentFile(bam_file, "rb")
    yield bam
    bam.close()


@pytest.fixture(scope="session")
def bed_file(tmp_path_factory):
    bed_out = tmp_path_factory.mktemp("data") / "test.bed"
    bed_out = str(bed_out)

    bed_data = "\t".join(("chr1", "0", "100")) + "\n"
    with open(bed_out, "w") as outf:
        outf.write(bed_data)

    return bed_out


# chromosome,position,reference_base,read_count,read_bases,base_qualities,map_qualities,base_positions
pileup_data = (
    "\n".join(
        (
            "chr1	4 G	195	.$.$............,,......................,...........,,,..,......,...,.,,..,..,...,..............,.......,,........,..,.,......,......,...,......,,...,..,..,.........,......,............,.,.....,,^].^],	ooJKtKKqKJJpwu?GvIwIMuuLvLqLLvqLLvKLLLKLLLtuLtILLLKKKIIMLfLLLLKLILKLLHLrKLMKKLKJLLLLLLLLLLLLLwKLLKqLLLJKLLLLLLLLKwLKLLKLMLLL9LLLLLLKKLLLLLLLLKLKKJBKLLJLLLKLKLKLILKKLILKIHKIJILJHIKKIJJIIHGGIJJDIDE	]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]	75,75,74,74,73,73,73,71,70,70,70,69,69,69,69,68,67,67,67,66,66,64,64,64,63,62,60,59,59,59,59,58,56,56,56,56,56,56,56,55,55,55,54,54,53,53,53,52,52,52,51,50,50,49,49,49,48,48,47,47,47,47,47,43,43,43,43,42,42,42,41,41,41,40,40,40,39,39,38,38,37,37,37,36,36,36,36,35,35,35,33,33,33,33,33,31,31,31,31,30,30,30,29,29,28,28,28,27,27,27,27,27,27,26,26,26,25,25,24,24,24,23,23,23,23,22,22,21,21,21,21,21,20,19,19,19,18,18,18,18,18,18,17,17,16,16,16,16,15,15,15,14,13,13,12,12,12,12,12,11,11,11,11,10,9,9,9,9,9,9,9,8,7,7,6,6,6,5,5,5,4,4,4,4,3,3,2,2,2,2,2,2,2,1,1",
            "chr1	5 T	198	.$.$..........,.,......................,,...........,,,..,......,...,.,,..,..,...,..............,.......,,........,..,,.,......,......,...,......,,...,..,..,.........,......,............,.,.....,,.,^].^].	DCkEEmFDDjjjEBEjDjDDiiDjDiDDiiDDjDDDDGGDDDjiDjDDDDGGFDDGD_DDDDEDDDGDGFDjFDDGDDDGDDDDDDDDDDDDDiFDDDjDDDFFCCDDhDDDFjDBFDGDDDDDD7DDDDDDFDDCEDDDDDDHGCCDEDDEDDGDDDDCCDDCDDCCCCDACCCCCCBCBCCDGBHCCBCCFHAFDD	]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]	75,75,74,74,74,72,71,71,71,70,70,70,70,69,69,68,68,68,67,67,65,65,65,64,63,61,60,60,60,60,59,57,57,57,57,57,57,57,57,56,56,56,55,55,54,54,54,53,53,53,52,51,51,50,50,50,49,49,48,48,48,48,48,44,44,44,44,43,43,43,42,42,42,41,41,41,40,40,39,39,38,38,38,37,37,37,37,36,36,36,34,34,34,34,34,32,32,32,32,31,31,31,30,30,29,29,29,28,28,28,28,28,28,27,27,27,27,26,26,25,25,25,24,24,24,24,23,23,22,22,22,22,22,21,20,20,20,19,19,19,19,19,19,18,18,17,17,17,17,16,16,16,15,14,14,13,13,13,13,13,12,12,12,12,11,10,10,10,10,10,10,10,9,8,8,7,7,7,6,6,6,5,5,5,5,4,4,3,3,3,3,3,3,3,2,2,1,1",
            "chr1	6 A	197	.$.$.$.......,.,......................,,...........,,,..,......,...,.,,..,..,...,..............,.......,,........,..,,.,......,......,...,......,,...,..,..,.........,......,............,.,.....,,.,..^].	nEEmIIItsvJ=IvHvKKvpLvJtLKvvKHxMJJKKKKKLvvKtIKKKLLLLKLLgLKKKLKKKLLL@LvLLLKKKKLJKLLLLLLKKLLKvLKKKtKKKLLKJKKvLKKKwLLLKMKLLLLK<KKJLKLLKKLJKLJKKKKLKHHLHKLKKLJKJKKLJKKMLJKKJHIKIJLJHJKJJLFJG?FFHIEJKEJ?HA	]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]	75,75,75,73,72,72,72,71,71,71,71,70,70,69,69,69,68,68,66,66,66,65,64,62,61,61,61,61,60,58,58,58,58,58,58,58,58,57,57,57,56,56,55,55,55,54,54,54,53,52,52,51,51,51,50,50,49,49,49,49,49,45,45,45,45,44,44,44,43,43,43,42,42,42,41,41,40,40,39,39,39,38,38,38,38,37,37,37,35,35,35,35,35,33,33,33,33,32,32,32,31,31,30,30,30,29,29,29,29,29,29,28,28,28,28,27,27,26,26,26,25,25,25,25,24,24,23,23,23,23,23,22,21,21,21,20,20,20,20,20,20,19,19,18,18,18,18,17,17,17,16,15,15,14,14,14,14,14,13,13,13,13,12,11,11,11,11,11,11,11,10,9,9,8,8,8,7,7,7,6,6,6,6,5,5,4,4,4,4,4,4,4,3,3,2,2,1",
            "chr1	7 C	198	.......,.,......................,,...........,,,..,......,...,.,,..,..,...,..............,.......,,,........,..,,.,......,......,...,......,,...,..,..,.........,......,............,.,.....,,.,...^].^].^].	nJJIsqpGHFpJoJKppJoJoJIqqKJqJJJKGGJJJpqJpJJJJFGGJJGIkJJJJHJIJGnGBIqGJJGJAJGJJJJJJJJJKJJJqGJIJpJJJGEGJJJJtKJJGpJGGJGJJJJJJ7JJJJIJGJJJGHHJJJJFFJHHGIJGJJGJJJJJJIJIGJIJIIIGIIIJIGIJJIIIGHGIIJJIHIJEJIIDAA	]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]	74,73,73,73,72,72,72,72,71,71,70,70,70,69,69,67,67,67,66,65,63,62,62,62,62,61,59,59,59,59,59,59,59,59,58,58,58,57,57,56,56,56,55,55,55,54,53,53,52,52,52,51,51,50,50,50,50,50,46,46,46,46,45,45,45,44,44,44,43,43,43,42,42,41,41,40,40,40,39,39,39,39,38,38,38,36,36,36,36,36,34,34,34,34,33,33,33,33,32,32,31,31,31,30,30,30,30,30,30,29,29,29,29,28,28,27,27,27,26,26,26,26,25,25,24,24,24,24,24,23,22,22,22,21,21,21,21,21,21,20,20,19,19,19,19,18,18,18,17,16,16,15,15,15,15,15,14,14,14,14,13,12,12,12,12,12,12,12,11,10,10,9,9,9,8,8,8,7,7,7,7,6,6,5,5,5,5,5,5,5,4,4,3,3,2,1,1,1",
        )
    )
    + "\n"
)


@pytest.fixture(scope="session")
def pileup_file(tmp_path_factory):
    p = tmp_path_factory.mktemp("data") / "test.pileup"
    p.write_text(pileup_data)
    return str(p)


@pytest.fixture(scope="session")
def gzipped_pileup_file(tmp_path_factory):
    p = tmp_path_factory.mktemp("data") / "test.pileup.gz"
    with gzip.open(p, "wt") as fp:
        fp.write(pileup_data)

    return str(p)


# not sure if I can use a larger scope w/out messing up sys.stdin
@pytest.fixture(scope="function")
def stdin_pileup_file(monkeypatch):
    monkeypatch.setattr("sys.stdin", io.StringIO(pileup_data))


@pytest.fixture(scope="function")
def stdin_bam_file(monkeypatch, bam_file):
    mock_stream = io.BytesIO(pysam.view("-b", bam_file))
    monkeypatch.setattr("sys.stdin", mock_stream)
