import heapq
from typing import List

from . import utils

import gzip


class CountFile:
    """Class to manage a count file."""

    def __init__(self, path):
        self.path = path
        self.is_gzipped = path.endswith(".gz")
        self.f = gzip.open(path, "rb") if self.is_gzipped else open(path, "r")
        self.f.readline()  # skip the header

        self.chrom = None
        self.pos = None
        self.ref = None
        self.alt = None
        self.bef = None
        self.aft = None
        self.ref_count = None
        self.alt_count = None
        self.ref_indel = None
        self.alt_index = None
        self.ref_fwd = None
        self.ref_rev = None
        self.alt_fwd = None
        self.alt_rev = None

        self.exhausted = False

        self.advance()

    def advance(self):
        """
        Advance to the next position
        Updates the CountFile object and returns True if succeeded
        Returns False if the file is exhausted
        """
        if self.exhausted:
            return False

        line = self.f.readline().decode() if self.is_gzipped else self.f.readline()

        # check for end of file
        if not line:
            self.exhausted = True
            self.f.close()
            return False

        split_line = line.strip().split("\t")

        self.chrom = split_line[0]
        self.pos = split_line[1]
        self.ref = split_line[2]
        self.alt = split_line[3]
        self.bef = split_line[4]
        self.aft = split_line[5]
        self.ref_count = CountFile.to_numeric(split_line[6])
        self.alt_count = CountFile.to_numeric(split_line[7])
        self.ref_indel = CountFile.to_numeric(split_line[8])
        self.alt_indel = CountFile.to_numeric(split_line[9])
        self.ref_fwd = CountFile.to_numeric(split_line[10])
        self.ref_rev = CountFile.to_numeric(split_line[11])
        self.alt_fwd = CountFile.to_numeric(split_line[12])
        self.alt_rev = CountFile.to_numeric(split_line[13])

        return True

    def pos_fields(self):
        return [
            self.chrom,
            self.pos,
            self.ref,
            self.alt,
            self.bef,
            self.aft,
        ]

    def count_fields(self):
        return [
            self.ref_count,
            self.alt_count,
            self.ref_indel,
            self.alt_indel,
            self.ref_fwd,
            self.ref_rev,
            self.alt_fwd,
            self.alt_rev,
        ]

    def __eq__(self, other):
        return self.chrom == other.chrom and self.pos == other.pos

    def __lt__(self, other):
        self_int_chrom = int(self.chrom.replace("chr", ""))
        other_int_chrom = int(other.chrom.replace("chr", ""))
        self_int_pos = int(self.pos)
        other_int_pos = int(other.pos)

        if self_int_chrom < other_int_chrom:
            return True

        elif self_int_chrom > other_int_chrom:
            return False

        # chroms must be the same here, so just compare positions
        return self_int_pos < other_int_pos

    def __repr__(self):
        return (
            f"{self.chrom}:{self.pos} "
            f"{self.bef}[{self.ref}>{self.alt}]{self.aft} "
            f"REF_COUNT={self.ref_count} ALT_COUNT={self.alt_count} "
            f"{self.path}"
        )

    @staticmethod
    def to_numeric(s):
        if s == ".":
            return 0
        else:
            return int(s)

    @staticmethod
    def from_numeric(n):
        if n == 0:
            return "."
        else:
            return str(n)


def merge(count_paths: List[str], out_path: str):
    """
    Merge count files into a single count file

    :param count_paths: paths of count files to merge
    :type count_paths: list[str]

    :param out_path: path to output file
    :type out_path: str
    """
    with utils.smart_open(out_path, "w") as out_f:

        header_fields = [
            "CHR",
            "POS",
            "REF",
            "ALT",
            "BEFORE",
            "AFTER",
            "REF_COUNT",
            "ALT_COUNT",
            "REF_INDEL",
            "ALT_INDEL",
            "REF_FWD",
            "REF_REV",
            "ALT_FWD",
            "ALT_REV",
        ]
        out_f.write("\t".join(header_fields) + "\n")

        # create a min-heap to keep track of the count files
        count_heap = []
        for p in count_paths:
            c = CountFile(p)
            if not c.exhausted:
                heapq.heappush(count_heap, c)

        # as long as there are positions in heap, keep going
        while count_heap:

            # create list of all count_files that are at the current position
            matched_count_files = [heapq.heappop(count_heap)]
            while count_heap and count_heap[0] == matched_count_files[0]:
                matched_count_files.append(heapq.heappop(count_heap))

            # get the set of the alternative alleles
            alt_bases = {c.alt for c in matched_count_files}

            # intersection with BASES={'A','C','G','T'} must be 0 or 1
            # will be zero if alt_bases={'.'} and one if alt_bases={'.','A'} for example
            if len(alt_bases & utils.BASES) < 2:
                # if the alt bases are {'.','A'}, we want the alt to be 'A'
                # if the alt bases are just {'.'} then we'll have it be '.'
                if len(alt_bases) > 1:
                    alt_bases -= {"."}

                alt_base = alt_bases.pop()

                # this is a bit ugly
                locus_info = matched_count_files[0].pos_fields()
                locus_info[3] = alt_base
                locus_info = "\t".join(locus_info)

                # adding up all the numerical fields from the count files
                count_fields = [c.count_fields() for c in matched_count_files]
                count_values = [sum(values) for values in zip(*count_fields)]
                count_strs = ["." if v == 0 else str(v) for v in count_values]
                count_info = "\t".join(count_strs)

                # write out to file
                out_f.write(locus_info + "\t" + count_info + "\n")

            # advance the count_files that are at this position
            # and add them back to the heap if they are not exhausted
            for c in matched_count_files:
                if c.advance():
                    heapq.heappush(count_heap, c)
