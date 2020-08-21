import argparse
import numpy as np
import pandas as pd
import pandas.api.types as pat
from scipy import stats

EXPRESSION_MATRIX_METADATA = ['Genotype', 'Genotype_Group', 'Replicate', 'Condition', 'tenXBarcode']
RANDOM_SEED = 42


def main():
    ap = argparse.ArgumentParser(description="Create a synthetic UMI count table")
    ap.add_argument("-e", "--expression_file", dest="expr_file", help="Expression data table", metavar="FILE",
                    default=None)
    ap.add_argument("-d", "--dist_file", dest="dist_file", help="Single-Cell Expression data table", metavar="FILE",
                    default=None)
    ap.add_argument("-o", "--out", dest="out", help="Output count table", metavar="FILE", required=True)
    ap.add_argument("--metadata_cols", dest="metadata", nargs="+", default=None)
    ap.add_argument("--samples_are_cols", dest="flip", help="Samples are columns; genes are rows", action='store_const',
                    const=True, default=False)
    ap.add_argument("--sample", dest="base_samp", help="Which column should be sampled from", default=0, type=int)
    ap.add_argument("--seed", dest="seed", help="Random Seed", metavar="NUM", default=RANDOM_SEED)
    ap.add_argument("--gzip", dest="gzip", help="gzip output", action='store_const', const=True, default=False)
    args = ap.parse_args()
    synthesize_data(args.expr_file, args.out, distribution_file_name=args.dist_file, meta_data_cols=args.metadata,
                    random_seed=args.seed, flip_matrix=args.flip, base_samp=args.base_samp, gzip=args.gzip)


def synthesize_data(expr_file_name, output_file_name, distribution_file_name=None, meta_data_cols=None,
                    random_seed=RANDOM_SEED, flip_matrix=False, base_samp=0, gzip=False):

    np.random.seed(random_seed)

    print("Reading expression data")

    expr_df = pd.read_csv(expr_file_name, sep="\t", header=0, index_col=0)

    object_cols = expr_df.dtypes == object

    if sum(object_cols) > 0:
        object_data = expr_df.loc[:, object_cols]
        expr_df.drop(expr_df.columns[object_cols], inplace=True, axis=1)
    else:
        object_data = None

    meta_data_cols = expr_df.columns.intersection(meta_data_cols) if meta_data_cols is not None else None

    if meta_data_cols is not None and len(meta_data_cols) > 0:
        object_data = pd.concat((object_data, expr_df.loc[:, meta_data_cols]))
        expr_df = expr_df.drop(meta_data_cols, axis=1)

    if flip_matrix:
        expr_df = expr_df.transpose()

    nrows, ncols = expr_df.shape

    if distribution_file_name is not None:
        dist_df = pd.read_csv(distribution_file_name, sep="\t", header=0, index_col=0)
        _, k = expr_df.shape
        assert k == 1

    else:
        assert base_samp < ncols
        dist_df = expr_df.iloc[base_samp, :].copy()

    if np.all([pat.is_integer_dtype(d) for d in expr_df.dtypes]):
        # Discrete sampling for count data
        umi = expr_df.sum(axis=1)
        dist_df = dist_df.divide(dist_df.sum()).tolist()
        synthetic_data = simulate_integer_data(dist_df, nrows, umi)

    else:
        gene_sd = expr_df.std(axis=0)
        gene_center = dist_df
        assert (gene_sd.index == gene_center.index).all()
        synthetic_data = simulate_floats_data(nrows, gene_center.tolist(), gene_sd.tolist())

    print("Writing Output")
    synthetic_data = pd.DataFrame(synthetic_data, index=expr_df.index, columns=expr_df.columns)

    if flip_matrix:
        synthetic_data = synthetic_data.transpose()

    if object_data is not None:
        synthetic_data = pd.concat((synthetic_data, object_data), axis=1)

    if gzip:
        synthetic_data.to_csv(output_file_name, sep="\t", compression="gzip")
    else:
        synthetic_data.to_csv(output_file_name, sep="\t")


def simulate_integer_data(prob_dist, nrows, n_per_row):
    if not np.isclose(np.sum(prob_dist), 1.):
        raise ValueError("Probability distribution does not sum to 1")
    ncols = len(prob_dist)
    print("Building Discrete Model")
    model = stats.rv_discrete(values=(range(ncols), prob_dist))
    synthetic_data = np.zeros((nrows, ncols), dtype=np.uint32)
    print("Simming Discrete Data")
    for i, u in enumerate(n_per_row):
        if i % 1000 == 0:
            print("\t[{i}/{tots}]".format(i=i, tots=nrows))
        reads = model.rvs(size=u)
        count_line = np.bincount(reads, minlength=ncols)
        synthetic_data[i, :] = count_line
    return synthetic_data


def simulate_floats_data(nrows, gene_centers, gene_sds):
    ncols = len(gene_centers)
    assert ncols == len(gene_sds)
    synthetic_data = np.zeros((nrows, ncols), dtype=np.float64)
    print("Simming Gaussian Data")
    for i, (gene_cen, gene_std) in enumerate(zip(gene_centers, gene_sds)):
        synthetic_data[:, i] = np.random.normal(loc=gene_cen, scale=gene_std, size=nrows)
    return synthetic_data


if __name__ == '__main__':
    main()
