import argparse
import os
import re
from collections import Counter, namedtuple
from string import Template

import numpy as np
import pysam
from scipy.stats import fisher_exact, mannwhitneyu

from typing import Dict, Iterable, Optional
import gzip
import sys

import warnings

import signal

from . import utils

signal.signal(signal.SIGPIPE, signal.SIG_DFL)

# Template to adhere to for output lines and header
out_t = Template(
    "$CHR\t$POS\t$REF\t$ALT\t$BEFORE\t$AFTER\t"
    "$REF_COUNT\t$ALT_COUNT\t$REF_INDEL\t$ALT_INDEL\t"
    "$REF_FWD\t$REF_REV\t$ALT_FWD\t$ALT_REV\t"
    "$P_STRAND_BIAS\t$P_BP_BIAS\t$P_BQ_BIAS\t"
    "$P_MQ_BIAS\t$MEAN_REF_BQ\t$MEAN_ALT_BQ\t"
    "$MEAN_REF_MQ\t$MEAN_ALT_MQ\t$REF_END5BP\t"
    "$ALT_END5BP\n"
)

BASES = {"A", "C", "G", "T"}

_PileupLine = namedtuple(
    "_PileupLine",
    ["chrom", "pos", "ref", "read_count", "bases", "bqs", "mqs", "bps", "qnames"],
)

indel_base_pattern = re.compile(r".[+-]\d+")
indel_number_pattern = re.compile(r"[+-](?P<length>\d+)")
start_pattern = re.compile(r"\^.")
end_pattern = re.compile(r"\$")


def count_pileup(fasta: str, seqlen: int, output, mpileup=None):
    """
    Process an mpileup file or input from stdin to create a bgzip counts file

    :param fasta: path to reference fasta file
    :type fasta: str
    :param seqlen: length of the sequencing reads
    :type seqlen: int
    :param output: path to write the output file
    :type output: str
    :param mpileup: path to mpileup file (if null, read from stdin instead)
    :type mpileup: str

    Outputs:

      - bgzip counts file (i.e. sample.gz)
      - tabix index file (i.e. sample.gz.tbi)
      - bed file of loci with alt allele count > 0 (i.e. sample.alt.txt)
    """

    # ensure that the output file ends with .gz
    if not output.endswith(".gz"):
        output += ".gz"

    unzipped_output = output[:-3]

    # create a position text file of loci with alt allele count > 0
    # replace the .gz at the end of the output file name with .alt.txt
    pos_path = unzipped_output + ".alt.txt"

    fasta = pysam.Fastafile(fasta)

    # Using open() for output instead of gzip.open() because we
    # compress the output file after writing it with bgzip  for tabix
    with utils.smart_open(mpileup, "r") as pipupfile, open(
        unzipped_output, "w"
    ) as of, open(pos_path, "w") as posf:
        # use the template as is for the header (just remove the $'s)
        output_header = out_t.template.replace("$", "")
        of.write(output_header)
        for line in pipupfile:
            # skip parsing mpileup header lines
            if line.startswith("#"):
                continue

            pileup = _PileupLine(*line.split())

            # Produce warning, but skip locus if reference base is not A, C, G, or T
            if pileup.ref not in BASES:
                warnings.warn(
                    f'WARNING: Unknown reference base "{pileup.ref}", must be in {BASES}, skipping locus {pileup.chrom}:{pileup.pos}'
                )
                continue

            # count all indel bases
            indel_counts = Counter(
                m[0] for m in indel_base_pattern.findall(pileup.bases)
            )

            base_string = parse_base_string(pileup.bases)
            # RB NOTE, test len(base_string) == pileup.read_count

            base_counts = Counter(base_string)

            # get index of ref_base (either period (plus strand) or comma (minus)
            ref_index = [i for i, b in enumerate(base_string) if b in {".", ","}]

            # skip locus without any REF reads
            # RB NOTE check against pileup.read_count instead
            if len(ref_index) == 0:
                continue

            # set of alt alleles
            alt_alleles = BASES.intersection(base_string.upper())
            n_alt = len(alt_alleles)

            if n_alt > 1:
                # skip loci with multiple alternative alleles
                continue

            if n_alt == 0:
                # Default values for invariant sites
                ALT = "."
                alt_fwd, alt_rev = 0, 0
                alt_indel = "."
                mean_ref_bq, mean_alt_bq, p_bq_bias_test = ".", ".", "."
                mean_ref_mq, mean_alt_mq, p_mq_bias_test = ".", ".", "."
                ref_end5bp, alt_end5bp, p_bp_bias_test = ".", ".", "."
                p_strand_bias_test = "."
            else:
                # Single alternative allele
                ALT = alt_alleles.pop()
                alt = ALT.lower()
                alt_index = [i for i, B in enumerate(base_string.upper()) if B == ALT]

                alt_fwd = base_counts[ALT]
                alt_rev = base_counts[alt]

                alt_indel = str(indel_counts[ALT] + indel_counts[alt])

                # Base quality, mapping quality, and bp bias stats
                mean_ref_bq, mean_alt_bq, p_bq_bias_test = get_qual_stat(
                    pileup.bqs, ref_index, alt_index
                )
                mean_ref_mq, mean_alt_mq, p_mq_bias_test = get_qual_stat(
                    pileup.mqs, ref_index, alt_index
                )
                ref_end5bp, alt_end5bp, p_bp_bias_test = get_bp_bias_stat(
                    pileup.bps, seqlen, ref_index, alt_index
                )

                # strand-bias test
                _, p_strand_bias_test = fisher_exact(
                    table=[
                        [base_counts["."], base_counts[","]],
                        [alt_fwd, alt_rev],
                    ],
                    alternative="two-sided",
                )

                # write txt file of loci with alt allele count > 0
                # NOTE consecutive loci with alt allele count > 0 are not merged into a single entry
                # I don't think this happens frequently, so the increase in file size is prob negligible
                posf.write(f"{pileup.chrom}\t{int(pileup.pos)}\n")

            # Get surrounding context
            before, after = get_context(fasta, pileup.chrom, int(pileup.pos))

            # Aggregate info for return
            locus = dict(
                CHR=pileup.chrom,
                POS=int(pileup.pos),
                REF=pileup.ref,
                ALT=ALT,
                BEFORE=before,
                AFTER=after,
                REF_COUNT=base_counts[","] + base_counts["."],
                REF_INDEL=indel_counts[","] + indel_counts["."],
                ALT_COUNT=alt_fwd + alt_rev,
                ALT_INDEL=alt_indel,
                REF_FWD=base_counts["."],
                REF_REV=base_counts[","],
                ALT_FWD=alt_fwd,
                ALT_REV=alt_rev,
                P_STRAND_BIAS=p_strand_bias_test,
                P_BP_BIAS=p_bp_bias_test,
                P_BQ_BIAS=p_bq_bias_test,
                P_MQ_BIAS=p_mq_bias_test,
                MEAN_REF_BQ=mean_ref_bq,
                MEAN_ALT_BQ=mean_alt_bq,
                MEAN_REF_MQ=mean_ref_mq,
                MEAN_ALT_MQ=mean_alt_mq,
                REF_END5BP=ref_end5bp,
                ALT_END5BP=alt_end5bp,
            )
            of.write(out_t.substitute(locus))

    # compress output file and create tabix index
    # then remove the uncompressed file
    # (would be better to directly write out in bgzip format, future todo)
    pysam.tabix_compress(unzipped_output, output, force=True)
    os.remove(unzipped_output)

    pysam.tabix_index(
        output,
        seq_col=0,
        start_col=1,
        end_col=1,
        line_skip=1,
        force=True,
    )


def parse_base_string(base_string: str) -> str:
    """Process a single pileup base string, removing indels and start/end sites.

    :param base_string: pileup base string
    :type base_string: str

    :return: pileup base string with indels and start/end sites removed
    """
    # remove start sites and score, leaving bases
    base_string, _ = start_pattern.subn("", base_string)

    # remove end sites
    base_string, _ = end_pattern.subn("", base_string)

    # remove indels
    base_string = remove_matches(
        base_string, indel_number_pattern.finditer(base_string)
    )

    return base_string


def remove_matches(string, matches):
    """Given string a finditer result, remove all matching substrings for indels
    assumes a group `length` used for removing the indel sequence

    :param string: pileup base string
    :type string: str
    :param matches: finditer result
    :type matches: re.Match

    :return: pileup base string with regex matches removed
    """
    result = []
    i = 0
    for match in matches:
        # add from last match to start of this match
        result.append(string[i : match.start()])
        length = int(match.group("length"))
        end = length + match.end()
        i = end
    result.append(string[i:])

    return "".join(result)


def get_qual_stat(qual_string, ref_ind, alt_ind):
    """compute P-values of base/mapping-quality-bias tests (two-sided Mann–Whitney U test), along with mean quality of ref and alt bases.

    :param qual_string: string of base/mapping qualities
    :type qual_string: str
    :param ref_ind: index of reference base
    :type ref_ind: list
    :param alt_ind: index of alternative base
    :type alt_ind: list

    :return: 3-tuple of mean reference base quality, mean alternative base quality, and P-value of base/mapping-quality-bias test
    """
    if not alt_ind:
        return ".", ".", "."

    ref_qual = get_qual(qual_string, ref_ind)
    mean_ref_qual = sum(ref_qual) / len(ref_qual)
    alt_qual = get_qual(qual_string, alt_ind)
    mean_alt_qual = sum(alt_qual) / len(alt_qual)
    p_qual_bias_test = get_p_mannwhitneyu(ref_qual, alt_qual)
    return mean_ref_qual, mean_alt_qual, p_qual_bias_test


def get_qual(qual_string, index):
    """Return substring (specified by index) of PHRED-scale value of a quality string (specified by qual_string).

    :param qual_string: PHRED-scale quality string
    :type qual_string: str
    :param index: index of substring
    :type index: list

    :return: list of PHRED-scale quality values
    """
    qual = [ord(qual_string[i]) - 33 for i in index]
    return qual


def get_bp_bias_stat(bp_string, seqlen, ref_ind, alt_ind):
    """compute P-values of base-positiion-bias tests (two-sided Mann–Whitney U test).

    :param bp_string: string of base positions
    :type bp_string: str
    :param seqlen: length of the sequencing read
    :type seqlen: int
    :param ref_ind: index of reference base
    :type ref_ind: list
    :param alt_ind: index of alternative base
    :type alt_ind: list

    :return: 3-tuple of number of reference bases within 5bp of read ends, number of alternative bases within 5bp of read ends, and P-value of base-position-bias test
    """
    if not alt_ind:
        return ".", ".", "."

    bp_string = bp_string.split(",")
    ref_dist = [max(seqlen - int(bp_string[i]), int(bp_string[i])) for i in ref_ind]
    alt_dist = [max(seqlen - int(bp_string[i]), int(bp_string[i])) for i in alt_ind]
    ref_end5bp = len([i for i in ref_dist if int(i) < 6])
    alt_end5bp = len([i for i in alt_dist if int(i) < 6])
    p_bp_bias_test = get_p_mannwhitneyu(ref_dist, alt_dist)
    return ref_end5bp, alt_end5bp, p_bp_bias_test


def get_p_mannwhitneyu(ref_string, alt_string):
    """compute P-values of two-sided Mann–Whitney U tests.

    :param ref_string: string of reference bases
    :type ref_string: str
    :param alt_string: string of alternative bases
    :type alt_string: str

    :return: P-value of two-sided Mann–Whitney U test
    """
    p_mannwhitneyu = "."
    if not (ref_string == alt_string) and not set(ref_string) == set(alt_string):
        _, p_mannwhitneyu = mannwhitneyu(ref_string, alt_string)
    return p_mannwhitneyu


def get_context(fasta, chrom: str, position: int):
    """Return reference base before and after locus provided by chrom and position.

    :param fasta: reference fasta file object
    :type fasta: pysam.FastaFile
    :param chrom: chromosome name
    :type chrom: str
    :param position: position of locus (1-based)
    :type position: int

    :return: 2-tuple of base before locus, and base after locus
    """
    # want surrounding two bases, fasta is zero indexed while position is unit
    before, _, after = fasta.fetch(chrom, position - 2, position + 1).upper()

    return before, after
