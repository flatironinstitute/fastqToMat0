import argparse
import numpy as np
import pandas as pd
from scipy import stats

EXPRESSION_MATRIX_METADATA = ['Genotype', 'Genotype_Group', 'Replicate', 'Condition', 'tenXBarcode']
RANDOM_SEED = 42

def main():
    ap = argparse.ArgumentParser(description="Create a synthetic UMI count table")
    ap.add_argument("-d", "--dist_file", dest="file", help="Expression data table", metavar="FILE", default=None)
    ap.add_argument("-s", "--ss_file", dest="ssfile", help="Single-Cell Expression data table", metavar="FILE",
                    required=True)
    ap.add_argument("-o", "--out", dest="out", help="Output count table", metavar="FILE", required=True)
    ap.add_argument("--log", dest="log", help="Data is log-transformed", action='store_const', const=True,
                    default=False)
    ap.add_argument("--shuffle", dest="shuffle", help="Don't simulate; just reshuffle", action='store_const',
                    const=True, default=False)

    args = ap.parse_args()

    synthesize_data(args.file, args.ssfile, args.out, dist_is_log=args.log, reshuffle_data=args.shuffle)


def synthesize_data(distribution_file_name, single_cell_file_name, output_file_name, dist_is_log=False,
                    reshuffle_data=False):

    np.random.seed(RANDOM_SEED)

    print("Reading single-cell data")
    ss_df = pd.read_csv(single_cell_file_name, sep="\t", header=0, index_col=0)

    meta_data = ss_df.loc[:, EXPRESSION_MATRIX_METADATA].copy()
    ss_df = ss_df.drop(EXPRESSION_MATRIX_METADATA, axis=1)
    nrows, ncols = ss_df.shape

    if reshuffle_data:
        pass
    else:
        umi = ss_df.drop(EXPRESSION_MATRIX_METADATA, axis=1).sum(axis=1)
        if distribution_file_name is not None:
            print("Reading distribution data")
            expr_df = pd.read_csv(distribution_file_name, sep="\t", header=0, index_col=0)
            if dist_is_log:
                expr_df = np.exp2(expr_df)
            print("Fixing Data")
            _, k = expr_df.shape
            assert k == 1
            cols = expr_df.columns
            expr_df = expr_df.divide(expr_df.sum()).tolist()
        else:
            cols = list(range(ss_df.shape[1]))
            expr_df = map(lambda x: x/ss_df.shape[1], [1.0] * ncols)

        synthetic_data = simulate_data(expr_df, nrows, umi)

    print("Writing Output")
    synth_df = pd.DataFrame(synthetic_data, index=ss_df.index, columns=cols)
    synth_df = pd.concat([synth_df, meta_data], axis=1)
    synth_df.to_csv(output_file_name, sep="\t", compression="gzip")



def simulate_data(prob_dist, nrows, n_per_row):
    if np.sum(prob_dist) != 1:
        raise ValueError("Probability distribution does not sum to 1")

    ncols = len(prob_dist)

    print("Building Model")

    model = stats.rv_discrete(values=(range(ncols), prob_dist))
    synthetic_data = np.zeros((nrows, ncols), dtype=np.uint32)

    print("Simming Data")

    for i, u in enumerate(n_per_row):
        if i % 1000 == 0:
            print("\t[{i}/{tots}]".format(i=i, tots=nrows))
        reads = model.rvs(size=u)
        count_line = np.bincount(reads, minlength=ncols)
        synthetic_data[i, :] = count_line

    return synthetic_data



if __name__ == '__main__':
    main()
