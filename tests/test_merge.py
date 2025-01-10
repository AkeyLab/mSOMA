import pytest

from msoma import merge


def test_merge_two_identical_files_unzipped(tmp_path):
    count_data_1 = (
        "CHR	POS	REF	ALT	BEFORE	AFTER	REF_COUNT	ALT_COUNT	REF_INDEL	ALT_INDEL	REF_FWD	REF_REV	ALT_FWD	ALT_REV	P_STRAND_BIAS	P_BP_BIAS	P_BQ_BIAS	P_MQ_BIAS	MEAN_REF_BQ	MEAN_ALT_BQ	MEAN_REF_MQ	MEAN_ALT_MQ	REF_END5BP	ALT_END5BP",
        "chr1	4	G	.	C	T	195	0	0	.	160	35	0	0	.	.	.	.	.	.	.	.	.	.",
        "chr1	5	T	.	G	A	198	0	0	.	161	37	0	0	.	.	.	.	.	.	.	.	.	.",
        "chr1	6	A	.	T	C	197	0	0	.	160	37	0	0	.	.	.	.	.	.	.	.	.	.",
        "chr1	7	C	.	A	A	198	0	0	.	160	38	0	0	.	.	.	.	.	.	.	.	.	.",
    )
    count_path_1 = tmp_path / "counts1.tsv"
    count_path_1.write_text("\n".join(count_data_1))

    count_data_2 = (
        "CHR	POS	REF	ALT	BEFORE	AFTER	REF_COUNT	ALT_COUNT	REF_INDEL	ALT_INDEL	REF_FWD	REF_REV	ALT_FWD	ALT_REV	P_STRAND_BIAS	P_BP_BIAS	P_BQ_BIAS	P_MQ_BIAS	MEAN_REF_BQ	MEAN_ALT_BQ	MEAN_REF_MQ	MEAN_ALT_MQ	REF_END5BP	ALT_END5BP",
        "chr1	4	G	.	C	T	195	0	0	.	160	35	0	0	.	.	.	.	.	.	.	.	.	.",
        "chr1	5	T	.	G	A	198	0	0	.	161	37	0	0	.	.	.	.	.	.	.	.	.	.",
        "chr1	6	A	.	T	C	197	0	0	.	160	37	0	0	.	.	.	.	.	.	.	.	.	.",
        "chr1	7	C	.	A	A	198	0	0	.	160	38	0	0	.	.	.	.	.	.	.	.	.	.",
    )
    count_path_2 = tmp_path / "counts2.tsv"
    count_path_2.write_text("\n".join(count_data_2))

    count_paths = [str(count_path_1), str(count_path_2)]

    out_path = str(tmp_path / "counts_merged.tsv")

    # run the merge
    merge.merge(count_paths, out_path)

    # check the output by reading the merged file
    with open(out_path) as f:
        lines = f.readlines()
        assert len(lines) == 5
        assert lines[0].startswith("CHR")
        assert lines[1].startswith("chr1\t4")
        assert lines[2].startswith("chr1\t5")
        assert lines[3].startswith("chr1\t6")
        assert lines[4].startswith("chr1\t7")

        assert lines[1].split("\t")[6] == "390"
        assert lines[2].split("\t")[6] == "396"
        assert lines[3].split("\t")[6] == "394"
        assert lines[4].split("\t")[6] == "396"
