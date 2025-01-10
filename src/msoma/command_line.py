import click
import subprocess
import importlib.resources
from . import count as bb_count
from . import pileup2counts as p2c
from . import qnames
from . import utils
from . import merge
from . import __version__


CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=__version__)
def msoma():
    """
    cli for mSOMA functionality
    """
    pass


@msoma.command()
@click.argument(
    "input_bam",
)
@click.option(
    "--fasta",
    type=click.Path(exists=True),
    help="Path to reference fasta file",
    required=True,
)
@click.option(
    "--output",
    "-o",
    help="Path to write counts file",
    default="-",
)
@click.option(
    "--bed",
    "-L",
    type=click.Path(exists=True),
    help="Path to bed file for callable regions",
    required=True,
)
@click.option(
    "--seq-length",
    "-l",
    type=int,
    help="Read length",
    required=True,
)
@click.option(
    "--min-MQ",
    "-q",
    type=int,
    help="Minimum mapping quality",
)
@click.option(
    "--require-flags",
    "-f",
    type=int,
    help="Require bitwise SAM flags",
)
@click.option(
    "--exclude-flags",
    "-F",
    type=int,
    help="Exclude bitwise SAM flags",
)
@click.option(
    "--ntrim",
    "-n",
    type=int,
    help="Number of bases to trim from both ends of the reads",
)
@click.option(
    "--max-indel",
    "-I",
    type=int,
    help="Maximum number of indel bases allowed in a read",
)
@click.option(
    "--mismatch-frac",
    "-m",
    type=float,
    help="Maximum fraction of mismatch bases allowed in a read",
)
@click.option(
    "--softclip-frac",
    "-s",
    type=float,
    help="Maximum fraction of softclipped bases allowed in a read",
)
@click.option(
    "--min-BQ",
    "-b",
    type=int,
    help="Minimum base quality for pileup",
)
@click.option(
    "--min-depth",
    "-d",
    type=int,
    help="Minimum depth for pileup",
)
@click.option(
    "--max-alt-allele",
    "-a",
    type=int,
    default=1,
    help="Maximum number of alternate alleles allowed per locus",
)
@click.option(
    "--qname-whitelist",
    "-N",
    type=click.Path(),
    help="Path to QNAME whitelist file, only reads with QNAMEs in this file will be counted",
    default=None,
)
@click.option(
    "--mpileup-intermediate",
    "-m",
    type=click.Path(),
    help="Path to write intermediate mpileup file for debugging, default is to not write",
    default=None,
)
def count(
    input_bam,
    fasta,
    output,
    bed,
    seq_length,
    min_mq,
    require_flags,
    exclude_flags,
    ntrim,
    max_indel,
    mismatch_frac,
    softclip_frac,
    min_bq,
    min_depth,
    max_alt_allele,
    qname_whitelist,
    mpileup_intermediate,
):
    """
    Create counts file from BAM file for somatic variant calling
    Use samtools and bamUtils piped implementation
    """
    bb_count.count(
        seq_path=input_bam,
        fasta_path=fasta,
        output_path=output,
        bed_path=bed,
        seq_len=seq_length,
        min_mapq=min_mq,
        require_flags=require_flags,
        exclude_flags=exclude_flags,
        ntrim=ntrim,
        max_indel=max_indel,
        mismatch_frac=mismatch_frac,
        softclip_frac=softclip_frac,
        min_bq=min_bq,
        min_depth=min_depth,
        max_alt_allele=max_alt_allele,
        qname_whitelist=qname_whitelist,
        mpileup_intermediate=mpileup_intermediate,
    )


@msoma.command(hidden=True)
@click.option(
    "--input",
    "-i",
    default="-",
    help="Path to pileup file",
)
@click.option(
    "--output",
    "-o",
    default="-",
    help="Path to write counts file",
)
@click.option(
    "--fasta",
    "-f",
    type=click.Path(exists=True),
    help="Path to reference fasta file",
    required=True,
)
@click.option(
    "--seq-length",
    "-l",
    type=int,
    help="Length of the sequence",
    required=True,
)
def pileup2counts(input, output, fasta, seq_length):
    """
    Convert from mpileup format to counts format
    """
    p2c.count_pileup(fasta=fasta, seqlen=seq_length, output=output, mpileup=input)


@msoma.command()
@click.argument(
    "input_counts",
)
@click.option(
    "--output",
    "-o",
    type=click.Path(),
    help="Path to write p-value output file",
    required=True,
)
@click.option(
    "--min-depth",
    "-d",
    type=int,
    help="Minimum read depth to consider locus for p-value calculation",
    required=True,
)
@click.option(
    "--ab",
    "-a",
    help="Output file to write alpha and beta parameter estimates",
    required=True,
)
def mle(input_counts, output, min_depth, ab):
    """
    Calculate p-values for each locus using maximum likelihood estimation from counts file

    NOTE uses an Rscript included as package-data in msoma
    """
    R_script_path = importlib.resources.files("msoma").joinpath("get_betabin.R")
    subprocess.check_output(
        [
            "Rscript",
            R_script_path,
            "--input",
            input_counts,
            "--output",
            output,
            "--mindp",
            str(min_depth),
            "--ab_file",
            ab,
        ]
    )


@msoma.command(hidden=True)
@click.argument(
    "input_bam",
)
@click.option(
    "--output",
    "-o",
    type=click.Path(),
    help="Path to write whitelisted QNAMEs file, default is to write to stdout",
)
@click.option(
    "--reference",
    "-r",
    type=click.Path(exists=True),
    help="Path to reference FASTA file, required iff input is CRAM",
)
@click.option(
    "--seq-length",
    "-l",
    type=int,
    help="Read length",
)
@click.option(
    "--min-MQ",
    "-q",
    type=int,
    help="Minimum mapping quality",
)
@click.option(
    "--require-flags",
    "-f",
    type=int,
    help="Require bitwise SAM flags",
)
@click.option(
    "--exclude-flags",
    "-F",
    type=int,
    help="Exclude bitwise SAM flags",
)
@click.option(
    "--max-indel",
    "-I",
    type=int,
    help="Maximum number of indel bases allowed in a read",
)
@click.option(
    "--mismatch-frac",
    "-m",
    type=float,
    help="Maximum fraction of mismatch bases allowed in a read",
)
@click.option(
    "--softclip-frac",
    "-s",
    type=float,
    help="Maximum fraction of softclipped bases allowed in a read",
)
def qname_whitelist(
    input_bam,
    output,
    reference,
    seq_length,
    min_mq,
    require_flags,
    exclude_flags,
    max_indel,
    mismatch_frac,
    softclip_frac,
):
    """
    Create a whitelist of QNAMEs from BAM file.
    This file can then be used to filter out reads from a BAM file during the count step.
    """
    qnames.whitelist(
        input_path=input_bam,
        output=output,
        reference=reference,
        seq_length=seq_length,
        min_MQ=min_mq,
        require_flags=require_flags,
        exclude_flags=exclude_flags,
        max_indel=max_indel,
        mismatch_frac=mismatch_frac,
        softclip_frac=softclip_frac,
    )


@msoma.command()
@click.argument(
    "count_paths",
    nargs=-1,
    required=True,
)
@click.option(
    "--output",
    "-o",
    type=click.Path(),
    help="Path to write merged counts file. Default is to write to stdout",
    default="-",
)
def merge_counts(count_paths, output):
    """
    Merge count files into a single count file
    """
    merge.merge(count_paths, output)

@msoma.command()
def check_dependencies():
    """
    Check if dependencies are installed and return version and path
    """
    executable_dependencies = ["Rscript", "samtools", "bam", "awk"]
    for dependency in executable_dependencies:
        utils.check_executable_dependency(dependency)

    # check R library dependencies
    r_dependencies = ["survcomp", "data.table", "VGAM", "argparse", "bbmle", "tidyverse", "Biostrings", "dplyr", "qvalue"]
    for dependency in r_dependencies:
        utils.check_R_library_dependency(dependency)
