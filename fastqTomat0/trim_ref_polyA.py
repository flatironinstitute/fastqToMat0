import argparse
import gzip
import os
import pathlib
import multiprocessing
import itertools
import re
import textwrap

from fastqTomat0.processor.fasta import fasta_gen


def main():
    ap = argparse.ArgumentParser(description="Remove polyA tracts from reference cDNA FASTA")
    ap.add_argument("-f", "--fasta", dest="fasta", help="FASTA file to modify", metavar="FILE", required=True)
    ap.add_argument("-o", "--out", dest="out", help="Output FASTA", metavar="FILE", required=True)
    ap.add_argument("--polyA_len", dest="pam", help="NUMBER As to consider polyA", metavar="NUMBER", type=int,
                    default=20)
    ap.add_argument("--polyA_shrink_len", dest="shrink", help="Shrink polyA tract to NUMBER bases", metavar="NUMBER",
                    type=int, default=5)
    args = ap.parse_args()

    poly_a_trimmer(fasta=args.fasta, output_fasta=args.out, tract_len=args.pam, shrink_len=args.shrink)


def poly_a_trimmer(fasta, output_fasta, tract_len, shrink_len, width=70):

    def _opener(x):
        return gzip.open if x.endswith(".gz") else open

    replace_re = re.compile("A{{{tl},}}".format(tl=tract_len), re.IGNORECASE)
    replace_str = "A" * shrink_len

    changes = []

    with _opener(fasta)(fasta) as fasta_fh, _opener(output_fasta)(output_fasta, mode="w") as output_fh:
        for record_id, record_seq in fasta_gen(fasta_fh):
            new_seq = replace_re.sub(replace_str, record_seq)

            if record_seq != new_seq:
                changes.append((record_id, len(record_seq) - len(new_seq)))

            print(">" + record_id, file=output_fh)
            print(textwrap.fill(new_seq, width=width), file=output_fh)

    print("Trimmed {n} sequences:".format(n=len(changes)))
    for record_id, bases_changed in changes:
        print("{g}: Removed {a} As".format(g=record_id, a=bases_changed))


if __name__ == '__main__':
    main()
