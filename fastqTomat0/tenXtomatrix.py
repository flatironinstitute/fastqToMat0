from fastqTomat0.lib import tenX
from fastqTomat0.lib.barcode_indexer import IDX, BARCODE, GENOTYPE

import pandas as pd
import argparse


def main():
    ap = argparse.ArgumentParser(description="Process 10x Files into a Matrix File")
    ap.add_argument("-p", "--path", dest="path", help="Input File PATH", metavar="PATH", required=True)
    ap.add_argument("-b", "--bc", dest="bc", help="Input Barcode-Genotype Map FILE", metavar="FILE", default=None)
    ap.add_argument("-i", "--bc_idx", dest="bc_idx", help="Barcode Library INDEX", metavar="INDEX", default=None)
    ap.add_argument("-o", "--output", dest="out", help="Output TSV FILE", metavar="FILE", required=True)
    ap.add_argument("--bulkup", dest="bulk", help="Unzip FASTQ", action='store_const', const=True, default=False)
    args = ap.parse_args()
    tenX_to_matrix(args.path, bc_file=args.bc, bc_file_lib_index=args.bc_idx, outfile_path=args.out,
                   bulk_up_genotypes=args.bulk)


def tenX_to_matrix(tenX_path, bc_file=None, bc_file_lib_index=None, outfile_path=None, remove_doublets=True,
                   bulk_up_genotypes=False):

    if bc_file is not None:
        bc = pd.read_table(bc_file, sep="\t")
        bc = bc.loc[bc[IDX] == bc_file_lib_index]
        bc.index = bc[BARCODE]
        if remove_doublets:
            bc.drop_duplicates(subset=BARCODE, keep=False)
        df = tenX.tenXProcessor(file_path=tenX_path, allowed_barcodes=bc.index.tolist()).process_files()
        df = filter_barcodes(df, bc)

        if bulk_up_genotypes:
            df = df.groupby(GENOTYPE).sum()

    else:
        df = tenX.tenXProcessor(file_path=tenX_path).process_files()

    if outfile_path is not None:
        df.to_csv(outfile_path, sep="\t")
    return df


def filter_barcodes(tenX_df, barcode_df):
    tenX_df = tenX_df.loc[barcode_df.index.intersection(tenX_df.index)]
    tenX_df = tenX_df.merge(barcode_df[[GENOTYPE]], left_index=True, right_index=True)
    return tenX_df


if __name__ == '__main__':
    main()
