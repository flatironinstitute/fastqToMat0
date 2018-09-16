from __future__ import print_function

import re
import gzip

from fastqTomat0.lib.degenerate_tools import make_regex_str, rc_str, convert_pattern
from fastqTomat0.lib.fastq import fastqProcessor

# These are the minimum quality scores required to consider a sequence as a valid barcode
BC_MIN_QUAL = 25
UMI_MIN_QUAL = 15
RC_TABLE = dict(A="T", G="C", T="A", C="G")


class ReadFromSeed:
    """
    Find a barcode using flanking sequences to identify location (common for amplicons)
    """

    seed = None  # str
    seed_len = None  # int
    bc_len = None  # int
    min_quality = None  # int

    allow_reverse = False  # bool

    regex = None  # Compiled re

    def __init__(self, seed, bc_length, min_quality=BC_MIN_QUAL, allow_reverse=False):

        self.bc_len = bc_length

        self.allow_reverse = allow_reverse
        self.regex = re.compile(make_regex_str(seed), re.IGNORECASE)
        self.seed = str(seed).upper().strip()
        self.seed_len = len(seed)
        self.min_quality = min_quality

    def scan_sequence(self, seq, qual):
        """
        For a sequence and quality score record, scan for the seed sequence and then extract the barcode if found
        :param seq: str
            Sequence
        :param qual: list(int)
            Quality scores
        :return bc:
            Return None if not found
        """
        try:
            return self._extract_barcode_regex(seq, qual)
        except (SeedNotFoundError, InsufficientSequence):
            pass
        except QualityError:
            return None

        if self.allow_reverse:
            try:
                return self._extract_barcode_regex(*self._rc(seq, qual))
            except (SeedNotFoundError, InsufficientSequence, QualityError):
                return None
        else:
            return None

    def _extract_barcode_regex(self, seq, qual):

        # Search the string with the regular expression
        seq = seq.upper()
        re_match = self.regex.search(seq)

        # There is no hit to the seed sequence
        if re_match is None:
            raise SeedNotFoundError

        # Get locations for the seed and barcode
        ss, se = re_match.start(), re_match.end()
        bcs, bce = se, se + self.bc_len

        # The sequence terminates before the BC can be fully extracted
        if bce > len(seq):
            raise InsufficientSequence

        # The quality score is too low
        if min(qual[bcs:bce]) < self.min_quality:
            raise QualityError

        return seq[bcs:bce]

    @staticmethod
    def _rc(seq, qual):
        """
        Reverse and complement a string and qual pairing
        """
        return rc_str(seq), qual[::-1]


class Link10xBCCounts:
    """
    Link the 10x genomics triple-read in order to connect 10x individual-cell barcodes with transcript-level barcode
    present in the genome
    """
    gz = False  # bool

    # First Barcode
    bc1_l = None  # int
    bc1_r = None  # int
    bc1_min_qual = None

    # UMI
    umi_l = None  # int
    umi_r = None  # int
    umi_min_qual = None  # int

    # Second Barcode
    scanner = None  # ReadFromSeed object
    bc2_min_qual = None  # int

    # Index
    index_min_qual = None  # int

    def __init__(self, bc1_pattern, bc2_len, bc2_seed, is_zipped=False, bc1_min_qual=BC_MIN_QUAL,
                 bc2_min_qual=BC_MIN_QUAL, umi_min_qual=UMI_MIN_QUAL, index_min_qual=UMI_MIN_QUAL):

        self.scanner = ReadFromSeed(bc2_seed, bc2_len, min_quality=bc2_min_qual)
        self.read_pattern(bc1_pattern)
        self.gz = is_zipped

        self.bc1_min_qual = bc1_min_qual
        self.bc2_min_qual = bc2_min_qual
        self.umi_min_qual = umi_min_qual
        self.index_min_qual = index_min_qual

    def parse_fastq(self, fastq1, fastq2, fastq3):
        totes, pf, uniques = 0, 0, 0
        bc_seen = {}

        fastq1_fh, fastq2_fh, fastq3_fh = self.open_file(fastq1), self.open_file(fastq2), self.open_file(fastq3)
        for fastq_records in fastqProcessor().fastq_gen(fastq1_fh, fastq2_fh, fastq3_fh):

            _, seq1, qual1 = fastq_records[0]
            _, seq2, qual2 = fastq_records[1]
            _, seq3, qual3 = fastq_records[2]

            if totes % 100000 == 0:
                print("Parsed {totes} sequence reads".format(totes=totes))

            totes += 1
            try:
                bc1, umi = self.parse_10x_bc(seq1, qual1)
                bc2 = self.parse_transcript_bc(seq2, qual2)
                idx = self.parse_index(seq3, qual3)
                pf += 1
            except (QualityError, BCNotFoundError, IndexError):
                continue

            try:
                bc_seen[idx][bc1][bc2].add(umi)
            except KeyError:
                nest_dict(bc_seen, idx, bc1, bc2)
                bc_seen[idx][bc1][bc2] = set(umi)

        fastq1_fh.close(), fastq2_fh.close(), fastq3_fh.close()
        return bc_seen

    def parse_fastq_mp(self, fq_tuple):
        fastq1, fastq2, fastq3 = fq_tuple
        return self.parse_fastq(fastq1, fastq2, fastq3)

    def parse_10x_bc(self, seq, qual):
        if min(qual[self.bc1_l:self.bc1_r]) < self.bc1_min_qual or min(qual[self.umi_l:self.umi_r]) < self.umi_min_qual:
            raise QualityError

        bc = seq[self.bc1_l:self.bc1_r]
        umi = seq[self.umi_l:self.umi_r]

        return bc, umi

    def parse_transcript_bc(self, seq, qual):
        bc = self.scanner.scan_sequence(seq, qual)
        if bc is None:
            raise BCNotFoundError
        else:
            return bc

    def parse_index(self, seq, qual):
        if min(qual) < self.index_min_qual:
            raise QualityError
        else:
            return seq

    def open_file(self, file_name):
        if self.gz:
            fh = gzip.open(file_name, mode="rt")
        else:
            fh = open(file_name, mode="rt")
        return fh

    def read_pattern(self, bc_pattern):
        pattern_idx = convert_pattern(bc_pattern)
        self.bc1_l, self.bc1_r = pattern_idx["B"]
        self.umi_l, self.umi_r = pattern_idx["U"]

        if self.bc1_l > self.umi_r and self.umi_l > self.bc1_r:
            raise ValueError("Incompatible pattern")


def nest_dict(d, *keys):
    dref = d
    for k in keys:
        try:
            dref[k]
        except KeyError:
            dref[k] = dict()
        dref = dref[k]


class QualityError(Exception):
    pass


class BCNotFoundError(Exception):
    pass


class SeedNotFoundError(Exception):
    pass


class InsufficientSequence(Exception):
    pass
