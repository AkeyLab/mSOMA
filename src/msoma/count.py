import logging

import subprocess
from . import __version__


def count(
    seq_path: str,
    output_path: str,
    fasta_path: str,
    bed_path: str,
    seq_len: int,
    **kwargs,
):
    """Generate REF and ALT counts for a given BAM or CRAM file.
    Uses samtools and bamUtils subprocess piped implementation.

    :param seq_path: path to BAM or CRAM file
    :type seq_path: str

    :param output_path: path to output file
    :type output_path: str

    :param fasta_path: path to reference FASTA file
    :type fasta_path: str

    :param bed_path: path to BED file
    :type bed_path: str

    :param seq_len: length of sequence
    :type seq_len: int

    Keyword Arguments:
        - min_mapq: (default 20)
        - require_flags: (default 2)
        - exclude_flags: (default 3844)
        - ntrim: (default 5)
        - max_indel: (default 0)
        - mismatch_frac: (default 0.1)
        - softclip_frac: (default 0.5)
        - min_bq: (default 20)
        - min_depth: (default 10)
        - max_alt_allele: (default 1)
        - qname_whitelist: (default None)
        - mpileup_intermediate: (default None)
    """
    # Create basic logger
    log_path = output_path + ".log"
    logging.basicConfig(
        format="%(asctime)s %(message)s",
        filename=log_path,
        encoding="utf-8",
        level=logging.DEBUG,
    )
    logging.info(f"MtSOMA version: {__version__}")
    logging.info(f"Starting count function")

    # NOTE maybe don't provide any default values?
    options = {
        "reference": None,
        "min_mapq": 20,
        "require_flags": 2,
        "exclude_flags": 3844,
        "ntrim": 5,
        "max_indel": 0,
        "mismatch_frac": 0.1,
        "softclip_frac": 0.5,
        "min_bq": 20,
        "min_depth": 10,
        "max_alt_allele": 1,
        "qname_whitelist": None,
        "mpileup_intermediate": None,
    }
    # remove None-valued kwargs
    kwargs = {k: v for k, v in kwargs.items() if v is not None}

    # check for unexpected kwargs
    unexpected = set(kwargs.keys()) - set(options.keys())
    if unexpected:
        raise TypeError(f"Unexpected keyword arguments: {unexpected}")

    options.update(kwargs)

    logging.info("Options: " + str(options))
    logging.info("BAM: " + seq_path)
    logging.info("FASTA: " + fasta_path)
    logging.info("BED: " + bed_path)
    logging.info("Output: " + output_path)
    logging.info("Seq Length: " + str(seq_len))

    if options["qname_whitelist"] is not None:
        logging.info("Using qname whitelist: " + options["qname_whitelist"])
    else:
        logging.info("Not using qname whitelist")

    if options["mpileup_intermediate"] is not None:
        logging.info(
            "Using intermediate mpileup file: " + options["mpileup_intermediate"]
        )
    else:
        logging.info("Not using intermediate mpileup file")

    SCLEN = int(options["softclip_frac"] * seq_len)
    MMLEN = int(options["mismatch_frac"] * seq_len)

    view_cmd = [
        "samtools",
        "view",
        "-h",
        "--region-file",
        bed_path,
        "-f",
        str(options["require_flags"]),
        "-F",
        str(options["exclude_flags"]),
        "-q",
        str(options["min_mapq"]),
        "-e",
        "(sclen<={SCLEN})&(([NM]<={MMLEN})||(!exists([NM])))".format(
            SCLEN=SCLEN, MMLEN=MMLEN
        ),
    ]

    if options["qname_whitelist"]:
        view_cmd.append("-N")
        view_cmd.append(options["qname_whitelist"])

    if seq_path.endswith(".cram"):
        view_cmd.append("--reference")
        view_cmd.append(fasta_path)

    view_cmd.append(seq_path)

    logging.debug(" ".join(view_cmd))

    indel_cmd = ["awk", '$1~"^@"||$6!~"I|D"']
    logging.debug(" ".join(indel_cmd))

    trim_cmd = ["bam", "trimBam", "-", "-", str(options["ntrim"])]
    logging.debug(" ".join(trim_cmd))

    # this is a debugging option to output intermediate mpileups info
    # don't need this for most cases
    if options["mpileup_intermediate"] is not None:
        mpileup_cmd = [
            "samtools",
            "mpileup",
            "-l",
            bed_path,
            "--min-BQ",
            str(options["min_bq"]),
            "--ignore-RG",
            "--fasta-ref",
            fasta_path,
            "--no-BAQ",
            "--output-BP",
            "--output-MQ",
            "--output-QNAME",
            "--output",
            options["mpileup_intermediate"],
            "-",
        ]
    else:
        mpileup_cmd = [
            "samtools",
            "mpileup",
            "-l",
            bed_path,
            "--min-BQ",
            str(options["min_bq"]),
            "--ignore-RG",
            "--fasta-ref",
            fasta_path,
            "--no-BAQ",
            "--output-BP",
            "--output-MQ",
            "--output-QNAME",
            "-",
        ]
    logging.debug(" ".join(mpileup_cmd))

    # if debugging with mpileup intermediate, then use that as input
    # otherwise, read from pipe
    if options["mpileup_intermediate"] is not None:
        p2c_cmd = [
            "msoma",
            "pileup2counts",
            "--output",
            output_path,
            "--fasta",
            fasta_path,
            "--seq-length",
            str(seq_len),
            "--input",
            options["mpileup_intermediate"],
        ]
    else:
        p2c_cmd = [
            "msoma",
            "pileup2counts",
            "--output",
            output_path,
            "--fasta",
            fasta_path,
            "--seq-length",
            str(seq_len),
            "--input",
            "-",
        ]
    logging.debug(" ".join(p2c_cmd))

    # Pipe each command to the next and then communicate to force execution
    p1 = subprocess.Popen(view_cmd, stdout=subprocess.PIPE)
    p2 = subprocess.Popen(indel_cmd, stdin=p1.stdout, stdout=subprocess.PIPE)
    p3 = subprocess.Popen(trim_cmd, stdin=p2.stdout, stdout=subprocess.PIPE)

    # if debugging with mpileup intermediate p4 outputs to file, not pipe
    if options["mpileup_intermediate"] is not None:
        p4 = subprocess.Popen(mpileup_cmd, stdin=p3.stdout)
        out, err = p4.communicate()
        p5 = subprocess.check_output(p2c_cmd)

    else:
        p4 = subprocess.Popen(mpileup_cmd, stdin=p3.stdout, stdout=subprocess.PIPE)
        p5 = subprocess.Popen(p2c_cmd, stdin=p4.stdout)
        out, err = p5.communicate()

    logging.info("Finished count function")
