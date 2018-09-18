from fastqTomat0.lib import tenX
from fastqTomat0.lib.barcode_indexer import IDX, BARCODE, GENOTYPE

import pandas as pd
import argparse

def main():
    ap = argparse.ArgumentParser(description="Process 10x Files into a Matrix File")
    ap.add_argument("-p", "--path", dest="path", help="Input File PATH", metavar="PATH", required=True)
    ap.add_argument("-b", "--bc", dest="bc", help="Input Barcode-Genotype Map FILE", metavar="FILE", default=None)
    ap.add_argument("-o", "--output", dest="out", help="Output TSV FILE", metavar="FILE", required=True)



def tenX_to_matrix(tenX_path, bc_file=None, bc_file_lib_index=None, outfile_path=None):
    df = tenX.tenXProcessor(file_path=tenX_path).process_files()

    if bc_file is not None:
        bc = pd.read_table(bc_file, sep="\t")
        bc = bc.loc[bc[IDX] == bc_file_lib_index]
        bc.index = bc[BARCODE]
        bc.drop(columns=BARCODE)
        df = filter_barcodes(df, bc)

    if outfile_path is not None:
        df.to_csv(outfile_path, sep="\t")
    return df


def filter_barcodes(tenX_df, barcode_df):
    tenX_df = tenX_df.loc[barcode_df.index.intersection(tenX_df.index)]
    tenX_df = tenX_df.merge(barcode_df[GENOTYPE])
    return tenX_df