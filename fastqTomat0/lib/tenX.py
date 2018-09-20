import numpy as np
import pandas as pd
from scipy.io import mmread


import os

class tenXProcessor:
    check_barcodes = False  # bool
    allowed_barcodes = None  # list(str)
    file_path = None # str

    gene_map = None # dict (str: str)
    barcode_list = None # list(str)

    def __init__(self, allowed_barcodes=None, file_path=None):
        if allowed_barcodes is not None:
            self.allowed_barcodes = allowed_barcodes
            self.check_barcodes = True
        self.file_path = file_path

    def process_files(self, gene_file="genes.tsv", barcode_file="barcodes.tsv", matrix_file="matrix.mtx"):
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
        self.gene_map = dict()
        with self.open_wrapper(gene_file, mode="r") as gene_fh:
            for line in gene_fh:
                try:
                    larr = line.strip().split()
                    self.gene_map[larr[0]] = larr[1]
                except IndexError:
                    continue

    def read_barcodes(self, barcode_file):
        self.barcode_list = list()
        with self.open_wrapper(barcode_file, mode="r") as bc_fh:
            for line in bc_fh:
                larr = line.strip().split("-")
                try:
                    self.barcode_list.append(larr[0])
                except IndexError:
                    continue

    def read_matrix(self, matrix_file):

        # Read the Market Matrix file
        with self.open_wrapper(matrix_file, mode="r") as mat_fh:
            data = mmread(mat_fh).todense().transpose()

        data = pd.DataFrame(data, index = self.barcode_list, columns = list(self.gene_map.keys()))
        has_reads = data.sum(axis=1) > 0
        data = data.loc[has_reads]

        # Remove any non-whitelisted barcodes if there's a whitelist
        if self.check_barcodes:
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
