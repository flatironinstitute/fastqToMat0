from __future__ import print_function

import argparse
import gzip
import os
import multiprocessing

from fastqTomat0.lib.fastq import fastq_gen
from fastqTomat0.lib.degenerate_tools import convert_pattern

BC_MIN_QUAL = 25
UMI_MIN_QUAL = 15
PATTERN = "BBBBBBBBBBBBBBBBUUUUUUUUUU"


def main():
    ap = argparse.ArgumentParser(description="Quick UMI/BC Counts from Single-Cell FastQs")
    ap.add_argument("-f", "--fastq", dest="fastq", help="Paired-End Read Barcode FASTQ", nargs="+", metavar="FILE",
                    required=True)
    ap.add_argument("-o", "--output", dest="out", help="Output TSV FILE", metavar="FILE", required=True)
    ap.add_argument("--gzip", dest="gzip", help="Unzip FASTQ", action='store_const', const=True, default=False)
    ap.add_argument("-c", "--cpu", dest="cores", help="NUMBER of cores to use", metavar="NUMBER", type=int, default=1)
    args = ap.parse_args()

    for fastq in args.fastq:
        if not os.path.exists(fastq):
            raise FileNotFoundError("{fq} File Not Found".format(fq=fastq))

    umi_count_per_bc(args.fastq, args.out, is_zipped=args.gzip, cores=args.cores)


def umi_count_per_bc(fastq_list, out_file_path, is_zipped=False, cores=1):
    pattern_idx = convert_pattern(PATTERN)
    bc_start, bc_stop = pattern_idx["B"]
    u_start, u_stop = pattern_idx["U"]

    if bc_start > u_stop and u_start > bc_stop:
        raise ValueError("Incompatible pattern")

    counter = FindUniqueCounts(bc_start, bc_stop, u_start, u_stop, is_zipped=is_zipped)
    pool = multiprocessing.Pool(processes=cores, maxtasksperchild=100)

    bc_pileup = {}

    for seq in pool.imap_unordered(counter.parse_fastq, fastq_list):
        for bc, umis in seq.items():
            try:
                bc_pileup[bc] = bc_pileup[bc].union(umis)
            except KeyError:
                bc_pileup[bc] = umis

    with open(out_file_path, mode="wt") as ofh:
        for bc, umi in bc_pileup.items():
            print("{bc}\t{num}".format(bc=bc, num=len(umi)), file=ofh)


class FindUniqueCounts:

    def __init__(self, bc_start, bc_stop, u_start, u_stop, is_zipped=False):
        self.bc_1, self.bc_2, self.u_1, self.u_2 = bc_start, bc_stop, u_start, u_stop
        self.gz = is_zipped

    def parse_fastq(self, fastq):
        totes, pf, uniques = 0, 0, 0
        bc_seen = {}

        fastq_fh = self.open_file(fastq)
        for cont, seq, qual in fastq_gen(fastq_fh):
            if totes % 100000 == 0:
                print("Parsed {totes} sequence reads ({pf} Pass Filter, {uni} Unique)".format(totes=totes,
                                                                                              uni=uniques,
                                                                                              pf=pf))
            totes += 1
            if min(qual[self.bc_1:self.bc_2]) < BC_MIN_QUAL or min(qual[self.u_1:self.u_2]) < UMI_MIN_QUAL:
                continue
            pf += 1

            bc = seq[self.bc_1:self.bc_2]
            umi = seq[self.u_1:self.u_2]

            try:
                bc_seen[bc].add(umi)
            except KeyError:
                bc_seen[bc] = set(umi)

        fastq_fh.close()
        return bc_seen

    def open_file(self, file_name):
        if self.gz:
            fh = gzip.open(file_name, mode="rt")
        else:
            fh = open(file_name, mode="rt")
        return fh


if __name__ == '__main__':
    main()
