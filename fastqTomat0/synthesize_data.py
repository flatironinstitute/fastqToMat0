import argparse
import numpy as np
import pandas as pd
from scipy import stats

EXPRESSION_MATRIX_METADATA = ['Genotype', 'Genotype_Group', 'Replicate', 'Condition', 'tenXBarcode']
RANDOM_SEED = 42

def main():
    ap = argparse.ArgumentParser(description="Create a synthetic UMI count table")
    ap.add_argument("-d", "--dist_file", dest="file", help="Expression data table", metavar="FILE", default=None)
    ap.add_argument("-s", "--ss_file", dest="ssfile", help="Single-Cell Expression data table", metavar="FILE", required=True)
    ap.add_argument("-o", "--out", dest="out", help="Output count table", metavar="FILE", required=True)
    ap.add_argument("--log", dest="log", help="Data is log-transformed", action='store_const', const=True, default=False)

    args = ap.parse_args()

    synthesize_data(args.file, args.ssfile, args.out, dist_is_log=args.log)


def synthesize_data(distribution_file_name, single_cell_file_name, output_file_name, dist_is_log=False):

    np.random.seed(RANDOM_SEED)

    print("Reading single-cell data")
    ss_df = pd.read_csv(single_cell_file_name, sep="\t", header=0, index_col=0)

    n, _ = ss_df.shape
    meta_data = ss_df.loc[:, EXPRESSION_MATRIX_METADATA].copy()
    umi = ss_df.drop(EXPRESSION_MATRIX_METADATA, axis=1).sum(axis=1)
    ss_df = None

    if distribution_file_name is not None:
        print("Reading distribution data")
        expr_df = pd.read_csv(distribution_file_name, sep="\t", header=0, index_col=0)
        if dist_is_log:
            expr_df = np.exp2(expr_df)
        print("Fixing Data")
        _, k = expr_df.shape
        assert k == 1
        expr_df = expr_df.divide(expr_df.sum()).tolist()
    else:
        expr_df = (np.ones((ss_df.shape[1], 1), dtype=np.dtype(float)) / ss_df.shape[1]).tolist()

    print("Building Model")

    model = stats.rv_discrete(values=(range(len(expr_df)), expr_df))
    synthetic_data = np.zeros((n, g), dtype=np.uint32)

    print("Simming Data")

    for i, u in enumerate(umi):
        if i % 1000 == 0:
            print("\t[{i}/{tots}]".format(i=i, tots=n))
        reads = model.rvs(size=u)
        count_line = np.bincount(reads, minlength=g)
        synthetic_data[i, :] = count_line

    print("Writing Output")

    synth_df = pd.DataFrame(synthetic_data, index=ss_df.index, columns=expr_df.index)
    synth_df = pd.concat([expr_df, meta_data], axis=1)
    synth_df.to_csv(output_file_name, sep="\t", compression="gzip")

if __name__ == '__main__':
    main()
