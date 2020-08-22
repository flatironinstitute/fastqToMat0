from __future__ import print_function

from fastqTomat0.processor.barcode_indexer import Link10xv31, create_10x_map_df

import argparse
import pandas as pd

# These are the minimum quality scores required to consider a sequence as a valid barcode
BC_MIN_QUAL = 25
UMI_MIN_QUAL = 15

# This is the pattern for 10x indexes
PATTERN = "BBBBBBBBBBBBBBBBUUUUUUUUUUUU"

# This is the pattern for transcript level genomic barcodes
BC2_LENGTH = 27
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
    ap.add_argument("-o", "--output", dest="out", help="Output TSV FILE", metavar="FILE", required=True)
    ap.add_argument("--gzip", dest="gzip", help="Unzip FASTQ", action='store_const', const=True, default=False)
    ap.add_argument("-c", "--cpu", dest="cores", help="NUMBER of cores to use", metavar="NUMBER", type=int, default=1)
    ap.add_argument("--bc_len", dest="bc_len", help="LEN of transcript barcode", metavar="LEN", type=int,
                    default=BC2_LENGTH)
    ap.add_argument("--bc_seed", dest="bc_seed", help="SEQ flanking transcript barcode (5' side)", metavar="SEQ",
                    default=BC2_SEED)
    ap.add_argument("--bc1_whitelist", dest="bc_white", help="Single-cell barcode whitelist", metavar="FILE",
                    default=BC2_SEED)
    ap.add_argument("--bc2_fasta", dest="bc_fasta", help="Transcript barcode whitelist", metavar="FILE",
                    default=BC2_SEED)

    args = ap.parse_args()

    link_barcodes(args.fastq1, args.fastq2, args.out, is_zipped=args.gzip,
                  bc2_seed=args.bc_seed, bc2_len=args.bc_len, bc1_whitelist=args.bc_white, bc2_mapfile=args.bc_fasta)


def link_barcodes(bc_fastq_1, bc_fastq_2, out_file_path=None, is_zipped=False,
                  max_index_mismatch=1, bc1_pattern=PATTERN, bc2_seed=BC2_SEED, bc2_len=BC2_LENGTH, bc2_mapfile=None,
                  include_unknowns=False, bc1_whitelist=None):

    assert len(bc_fastq_1) == len(bc_fastq_2)

    print("Loading whitelists")
    # Load whitelists
    bc1_whitelist = pd.read_csv(bc1_whitelist, header=None).iloc[:, 0].tolist() if bc1_whitelist is not None else None
    bc2_map = pd.read_csv(bc2_mapfile, sep="\t", index_col=0) if bc2_mapfile is not None else None
    bc2_whitelist = bc2_map.index.tolist() if bc2_map is not None else None

    print("Creating BC mapper")
    # Create a transcript to cell barcode linking map
    linker = Link10xv31(bc1_pattern, bc2_len, bc2_seed, is_zipped=is_zipped, bc1_whitelist=bc1_whitelist,
                        bc2_whitelist=bc2_whitelist)

    print("Scanning FASTQs")
    bcs = linker.parse_fastq_mp(bc_fastq_1, bc_fastq_2)

    # Convert the output to a dataframe
    print("Identified {n} cell-specific barcodes".format(n=len(bcs)))
    bc_df = create_10x_map_df(bcs, bc2_map=bc2_map)

    if out_file_path is not None:
        bc_df.to_csv(out_file_path, sep="\t", index=False)
    return bc_df


if __name__ == '__main__':
    main()
