from __future__ import unicode_literals, print_function

import argparse
import glob
import gzip
import multiprocessing
import os
import subprocess

import pandas as pd

from fastqTomat0.processor.degenerate_tools import search_list_degenerate, convert_pattern

PATTERN = "BBBBBBBBBBBBBBBBUUUUUUUUUU"


def main():
    ap = argparse.ArgumentParser(description="Coerce non-10x barcodes to look like 10x barcodes (to test cellranger)")
    ap.add_argument("-f", "--fastq", dest="fastq", help="Paired-End Read Barcode FASTQ", nargs="+", metavar="FILE",
                    required=True)
    ap.add_argument("-b", "--barcodes", dest="barcode", help="Translation table (Current BC -> 10x BC)",
                    metavar="FILE", required=True)
    ap.add_argument("-c", "--cpu", dest="cores", help="NUMBER of cores to use", metavar="NUMBER", type=int, default=1)
    args = ap.parse_args()

    file_names = []
    for fastq in args.fastq:
        file_names.extend(glob.glob(fastq))

    for fn in file_names:
        if not os.path.exists(fn):
            raise FileNotFoundError("{fq} File Not Found".format(fq=fastq))

    translate_barcodes(file_names, args.barcode, cores=args.cores)


def translate_barcodes(fastq_files, barcode_file, cores=1):
    pattern_idx = convert_pattern(PATTERN)
    bc_start, bc_stop = pattern_idx["B"]
    u_start, u_stop = pattern_idx["U"]

    bc_table = pd.read_csv(barcode_file, sep="\t", index_col=0)
    bc_dict = bc_table.iloc[:, 0].to_dict()

    fixer = ReplaceBarcodes(bc_dict, bc_start, bc_stop, u_start, u_stop)
    pool = multiprocessing.Pool(processes=cores, maxtasksperchild=100)

    output = []
    for parsed, changed, filename in pool.imap_unordered(fixer.read_and_replace, fastq_files):
        output.append("Parsed {l} lines (Changed {ch}) in file {f}".format(l=parsed, ch=changed, f=filename))

    print("\n".join(output))


class ReplaceBarcodes:

    def __init__(self, bc_dict, bc_start, bc_stop, u_start, u_stop):
        self.bc_dict = bc_dict
        self.bc_1, self.bc_2, self.u_1, self.u_2 = bc_start, bc_stop, u_start, u_stop

    def read_and_replace(self, fastq):
        of_path, of_file_name = os.path.split(fastq)
        os.makedirs(os.path.join(of_path, "BCFIX"), exist_ok=True)
        of_file_name = of_file_name[:-3] if of_file_name.endswith(".gz") else of_file_name
        out_file_name = os.path.join(of_path, "BCFIX", of_file_name)
        allowed_barcodes = list(self.bc_dict.keys())
        with self.open_file(fastq) as infh, self.open_file(out_file_name, mode="wt") as outfh:
            sequence_line = False
            sequences_printed = 0
            sequences_changed = 0
            print("Loading data from file {fl}".format(fl=fastq))
            for line in infh.read().splitlines():
                line = line.strip()
                if len(line) == 0:
                    print(line, file=outfh)
                elif line[0] == "+":
                    print(line, file=outfh)
                elif line[0] == "@":
                    sequence_line = True  # Next line is a sequence
                    print(line, file=outfh)
                elif sequence_line:
                    sequence_line = False  # Last line was an ID

                    if sequences_printed % 1000000 == 0:
                        print("Parsed {totes} sequence reads".format(totes=sequences_printed))

                    sequences_printed += 1
                    bc = line[self.bc_1:self.bc_2]
                    umi = line[self.u_1:self.u_2]

                    try:
                        new_bc = self.bc_dict[bc]
                        sequences_changed += 1
                        print(new_bc + umi, file=outfh)
                    except KeyError:
                        print(line, file=outfh)
                else:
                    print(line, file=outfh)

        print("Zipping output file {fi}".format(fi=out_file_name))
        subprocess.call(["gzip", os.path.abspath(os.path.expanduser(os.path.expandvars(out_file_name)))])
        return sequences_printed, sequences_changed, fastq


    def open_file(self, file_name, mode="rt"):
        file_name = os.path.abspath(os.path.expanduser(os.path.expandvars(file_name)))

        if file_name.endswith(".gz"):
            fh = gzip.open(file_name, mode=mode)
        else:
            fh = open(file_name, mode=mode)
        return fh


if __name__ == '__main__':
    main()
