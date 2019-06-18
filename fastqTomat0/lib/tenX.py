import numpy as np
import pandas as pd
from scipy.io import mmread

import os

GENE_FILE = "genes.tsv"
BARCODE_FILE = "barcodes.tsv"
MATRIX_FILE = "matrix.mtx"

class tenXProcessor:
    check_barcodes = False  # bool
    allowed_barcodes = None  # list(str)
    file_path = None  # str

    gene_list = None  # list (str)
    barcode_list = None  # list(str)

    def __init__(self, allowed_barcodes=None, file_path=None):
        if allowed_barcodes is not None:
            self.allowed_barcodes = allowed_barcodes
        self.file_path = file_path

    def process_files(self, gene_file=GENE_FILE, barcode_file=BARCODE_FILE, matrix_file=MATRIX_FILE):
        """

        :param gene_file: str
        :param barcode_file: str
        :param matrix_file: str
        :return: pd.DataFrame
            Dataframe indexed by barcode (str) with genes as columns. Values are UMI counts.
        """
        self.read_genes(gene_file)
        self.read_barcodes(barcode_file)
        return self.read_matrix(matrix_file)

    def read_genes(self, gene_file):
        with self.open_wrapper(gene_file, mode="r") as gene_fh:
            self.gene_list = pd.read_table(gene_fh, header=None).iloc[:, 0].tolist()

    def read_barcodes(self, barcode_file):
        with self.open_wrapper(barcode_file, mode="r") as bc_fh:
            self.barcode_list = pd.read_table(bc_fh, header=None).iloc[:, 0].str.replace("-1", "").tolist()

    def read_matrix(self, matrix_file):

        # Read the Market Matrix file
        with self.open_wrapper(matrix_file, mode="rb") as mat_fh:
            data = mmread(mat_fh).todense()

        if data.shape[1] == len(self.barcode_list) and data.shape[0] != len(self.barcode_list):
            data = data.transpose()
        data = pd.DataFrame(data, index=self.barcode_list, columns=self.gene_list)

        # Remove any barcodes with absolutely no reads
        has_reads = data.sum(axis=1) > 0
        data = data.loc[has_reads]

        # Remove any non-whitelisted barcodes if there's a whitelist
        if self.allowed_barcodes is not None:
            keepers = data.index.intersection(self.allowed_barcodes)
            data = data.loc[keepers]

        return data

    def open_wrapper(self, file_name, mode="r"):
        if self.file_path is not None:
            file_name = os.path.join(self.file_path, file_name)

        if file_name.endswith(".bz2"):
            import bz2
            return bz2.BZ2File(file_name, mode=mode)
        elif file_name.endswith(".gz"):
            import gzip
            return gzip.open(file_name, mode=mode)
        else:
            return open(file_name, mode=mode)
