from __future__ import print_function

from fastqTomat0.lib.fasta import fasta_to_dict
from fastqTomat0.lib.barcode_indexer import Link10xBCCounts, create_10x_genotype_df

import argparse
import multiprocessing

import pandas as pd

# These are the minimum quality scores required to consider a sequence as a valid barcode
BC_MIN_QUAL = 25
UMI_MIN_QUAL = 15

# This is the pattern for 10x indexes
PATTERN = "BBBBBBBBBBBBBBBBUUUUUUUUUU"

# This is the pattern for transcript level genomic barcodes
BC2_LENGTH = 19
BC2_SEED = 'ttttctaaggatcca'

UNKNOWN_STR = "UNK"

# Pandas column names
IDX, BARCODE, GENOTYPE, UMI_COUNT = 'Library_Index', 'Cell_Barcode', 'Genotype', 'Umi_Count'


def main():
    ap = argparse.ArgumentParser(description="Extract Secondary Barcodes From 10x Data")
    ap.add_argument("-1", "--fastq1", dest="fastq1", help="Cell Index (26bp) FASTQ", nargs="+", metavar="FILE",
                    required=True)
    ap.add_argument("-2", "--fastq2", dest="fastq2", help="Transcript Index (98bp) FASTQ", nargs="+", metavar="FILE",
                    required=True)
    ap.add_argument("-3", "--fastq3", dest="fastq3", help="Library Index (8bp) FASTQ", nargs="+", metavar="FILE",
                    required=True)
    ap.add_argument("-i", "--indexes", dest="indexes", help="Indexes to extract", nargs="+", metavar="SEQ",
                    required=True)
    ap.add_argument("-m", "--mismatch", dest="idx_mismatch", help="Number of allowed mismatches in index",
                    metavar="NUM", type=int, default=1)
    ap.add_argument("-o", "--output", dest="out", help="Output TSV FILE", metavar="FILE", required=True)
    ap.add_argument("--gzip", dest="gzip", help="Unzip FASTQ", action='store_const', const=True, default=False)
    ap.add_argument("-c", "--cpu", dest="cores", help="NUMBER of cores to use", metavar="NUMBER", type=int, default=1)
    ap.add_argument("--bc_len", dest="bc_len", help="LEN of transcript barcode", metavar="LEN", type=int,
                    default=BC2_LENGTH)
    ap.add_argument("--bc_seed", dest="bc_seed", help="SEQ flanking transcript barcode (5' side)", metavar="SEQ",
                    default=BC2_SEED)

    args = ap.parse_args()

    link_barcodes(args.fastq1, args.fastq2, args.fastq3, args.out, is_zipped=args.gzip, allowed_indexes=args.indexes,
                  bc2_seed=args.bc_seed, bc2_len=args.bc_len)


def link_barcodes(bc_fastq_1, bc_fastq_2, bc_fastq_3, out_file_path=None, is_zipped=False, allowed_indexes=None,
                  max_index_mismatch=1, bc1_pattern=PATTERN, bc2_seed=BC2_SEED, bc2_len=BC2_LENGTH, bc2_mapfile=None):
    assert len(bc_fastq_1) == len(bc_fastq_2)

    if bc2_mapfile is not None:
        with open(bc2_mapfile, mode="r") as mapfh:
            bc2_map = fasta_to_dict(mapfh, key_type='seq')
    else:
        bc2_map = None

    linker = Link10xBCCounts(bc1_pattern, bc2_len, bc2_seed, is_zipped=is_zipped)
    bcs = linker.parse_fastq_mp(bc_fastq_1, bc_fastq_2, bc_fastq_3)
    bc_df = create_10x_genotype_df(bcs, allowed_indexes=allowed_indexes, max_index_mismatch=max_index_mismatch,
                                   bc2_map=bc2_map, include_unknowns=False)
    bc_df.drop_duplicates(subset=BARCODE, keep=False)

    if out_file_path is not None:
        bc_df.to_csv(out_file_path, sep="\t")
    return bc_df


if __name__ == '__main__':
    main()
