import logging
import pysam

from . import __version__
from . import utils


def whitelist(input_path: str, output=None, **kwargs):
    """
    Generate a whitelist of QNAMEs from a BAM or CRAM file

    :param input_path: Path to input BAM/CRAM file
    :type input_path: str

    :param output: Path to output whitelist file, default is to write to stdout if not provided
    :type output: str or None

    Keyword Args:
        - reference (str): Path to reference FASTA file, required iff input is CRAM
        - seq_length (int): Read length, if not provided, will be calculated from each read
        - min_MQ (int): Minimum mapping quality
        - require_flags (int): Require bitwise SAM flags
        - exclude_flags (int): Exclude bitwise SAM flags
        - max_indel (int): Maximum number of indel bases allowed in a read
        - mismatch_frac (float): Maximum fraction of mismatch bases allowed in a read
        - softclip_frac (float): Maximum fraction of softclipped bases allowed in a read
    """
    options = {
        "reference": None,
        "seq_length": None,
        "min_MQ": 0,
        "require_flags": 0,
        "exclude_flags": 0,
        "max_indel": None,
        "mismatch_frac": 1,
        "softclip_frac": 1,
    }
    options.update(kwargs)

    failing_qnames = set()

    if input_path.endswith(".bam"):
        fhandle = pysam.AlignmentFile(input_path, "rb")
    elif input_path.endswith(".cram") and options["reference"] is not None:
        fhandle = pysam.AlignmentFile(
            input_path, "rc", reference_filename=options["reference"]
        )
    elif input_path.endswith(".cram") and options["reference"] is None:
        raise ValueError("CRAM file requires a reference FASTA file")
    else:
        raise ValueError(f"Input file must be a BAM or CRAM file: {input_path}")

    for r in fhandle:
        qname = r.query_name

        # If the QNAME has previously failed, skip
        if qname in failing_qnames:
            continue

        # NOTE this is slow to do for every read
        read_length = options["seq_length"] or r.query_length

        # NOTE this is slow to do for every read especially if not filtering on mismatch_frac, softclip_frac, etc
        # NOTE compose a list of all the filters that need to be applied with functools?
        # Check if read is passing all the filters
        # will stop early if any filter fails
        cigar_counts, _ = r.get_cigar_stats()
        passes_all_filters = (
            r.mapping_quality >= options["min_MQ"]
            and r.flag & options["require_flags"] == options["require_flags"]
            and r.flag & options["exclude_flags"] == 0
            and (
                options["max_indel"] == None
                or cigar_counts[1] + cigar_counts[2] <= options["max_indel"]
            )
            and r.get_tag("NM") / read_length <= options["mismatch_frac"]
            and cigar_counts[4] / read_length <= options["softclip_frac"]
        )

        # If the read failed any of the filters
        # then add to the failing QNAMEs set
        # and remove from the passing QNAMEs set if present
        if not passes_all_filters:
            failing_qnames.add(qname)

    # Close the file handle
    fhandle.close()

    # Read through the file again, writing out QNAMEs not in failing_qnames to the output file
    if input_path.endswith(".bam"):
        fhandle = pysam.AlignmentFile(input_path, "rb")
    elif input_path.endswith(".cram") and options["reference"] is not None:
        fhandle = pysam.AlignmentFile(
            input_path, "rc", reference_filename=options["reference"]
        )
    elif input_path.endswith(".cram") and options["reference"] is None:
        raise ValueError("CRAM file requires a reference FASTA file")
    else:
        raise ValueError(f"Input file must be a BAM or CRAM file: {input_path}")

    # NOTE, this can result in duplicate QNAMEs in the output file
    # but this is not a problem for the downstream count step and
    # it is much faster than checking if the QNAME has already been written
    with utils.smart_open(output, "w") as out_f:
        for r in fhandle:
            qname = r.query_name
            if qname not in failing_qnames:
                out_f.write(f"{qname}\n")

    # Close the file handle
    fhandle.close()
