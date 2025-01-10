import collections
import gzip
import re
import sys
import click
from contextlib import contextmanager
from typing import Dict, Iterable, Optional, Tuple
import subprocess

import pandas as pd
import pysam
import scipy as scp


BASES = {"A", "C", "G", "T"}

# regex to match cigar string operations
# i.e. 10M1I25M -> [(10, M), (1, I), (25, M)]
cigar_re = re.compile(r"(\d+)([MIDNSHP=X])")


@contextmanager
def smart_open(filename: Optional[str] = None, mode: str = "r"):
    """A wrapper for stdin/open/gzip.open logic as a context manager

    :param filename: filename to open or None for stdin/stdout
    :type filename: str or None

    :param mode: mode to open file in, either "r" or "w"
    :type mode: str
    """
    if mode not in {"r", "w"}:
        raise ValueError('ERROR: smart_open mode must be "r" or "w"')

    # stdin/stdout case
    if filename is None or filename == "-":
        if mode == "r":
            yield sys.stdin
        else:
            yield sys.stdout

    # file case
    else:
        if filename.endswith(".gz"):
            open_file = gzip.open(filename, mode + "t")
        else:
            open_file = open(filename, mode)

        yield open_file
        open_file.close()


def get_context(fasta, chrom: str, position: int) -> Tuple[str, str, str]:
    """Return reference base before and after locus provided by chrom and position

    :param fasta: reference fasta file object
    :type fasta: pysam.FastaFile
    :param chrom: chromosome name
    :type chrom: str
    :param position: position of locus (1-based)
    :type position: int

    :return: 3-tuple of base before locus, base at locus, base after locus
    """
    # want surrounding two bases, fasta is zero indexed while position is unit
    # bases are returned as upper case since in ref to plus strand
    BEF, REF, AFT = fasta.fetch(chrom, position - 2, position + 1).upper()
    return BEF, REF, AFT


def sim_betabinom(ns: Iterable[int], a: float, b: float) -> pd.DataFrame:
    """Simulate successes from events with distribution
    following a beta-binomial k_{i} ~ BetaBinom(n_{i}, a, b)

    :param ns: is a list/array/tuple of total trials (ints)
    :type ns: Iterable[int]
    :param a: is a float for the alpha param to the betabinom
    :type a: float
    :param b: is a float for the beta param to the betabinom
    :type b: float

    :return: df is a pandas dataframe with columns n (total trials), k (successes), j (failures)
    :rtype: pd.DataFrame
    """
    df = pd.DataFrame(
        {
            "n": ns,
            "k": scp.stats.betabinom.rvs(n=ns, a=a, b=b),
        }
    )
    df["j"] = df["n"] - df["k"]
    return df


def cigarstring_counts(cigarstring: str) -> Dict[str, int]:
    """Parse a cigarstring and return counts of each cigar operation
    Sums the number of each operation in the cigarstring

    :param cigarstring: cigarstring to parse such as "10M1I25M"
    :type cigarstring: str

    :return: counts of cigar operation counts
    :rtype: dict[str, int]

    Example:
    cigarstring_counts('10M1I25M') -> {'M': 35, 'I': 1}

    Allowed operations are:

        - M: match
        - I: insertion to the reference
        - D: deletion from the reference
        - N: skipped region from the reference
        - S: soft clipping (clipped sequences present in SEQ)
        - H: hard clipping (clipped sequences NOT present in SEQ)
        - P: padding (silent deletion from padded reference)
        - =: sequence match
        - X: sequence mismatch
    """
    matches = cigar_re.findall(cigarstring)
    counts = collections.defaultdict(int)  # type: Dict[str, int]
    parsed_cigar_len = 0

    for count, kind in matches:
        counts[kind] += int(count)
        parsed_cigar_len += len(count) + 1  # +1 for the operation character

    # Every character in the cigarstring should be parsed, any missing means
    # that the cigarstring was not parsed correctly
    if len(cigarstring) != parsed_cigar_len:
        raise ValueError(f"ERROR: cigarstring {cigarstring} not parsed correctly")

    return counts

def check_executable_dependency(dependency):
    yes_found = click.style('Dependency found   external : ', fg='green')
    not_found = click.style('Dependency missing external : ', fg='red')
    try:
        dependency_path_result = subprocess.check_output(["which", dependency], stderr=subprocess.STDOUT)
        dependency_path = dependency_path_result.decode().strip()
        click.echo(yes_found + f"{dependency}: path: {dependency_path}")
    except subprocess.CalledProcessError:
        click.echo(not_found + f"{dependency} not found in PATH")

def check_R_library_dependency(dependency):
    yes_found = click.style('Dependency found   R-library: ', fg='green')
    not_found = click.style('Dependency missing R-library: ', fg='red')
    try:
        subprocess.check_output(["Rscript", "-e", f"library({dependency})"], stderr=subprocess.STDOUT)
        click.echo(yes_found + f"{dependency}")
    except subprocess.CalledProcessError:
        click.echo(not_found + f"{dependency}")