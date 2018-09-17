from __future__ import print_function

from fastqTomat0.lib.fasta import fasta_to_dict
from fastqTomat0.lib.degenerate_tools import sequence_mismatch_degenerate, search_list_degenerate
from fastqTomat0.lib.barcode_indexer import Link10xBCCounts, nest_dict

import argparse
import multiprocessing

# These are the minimum quality scores required to consider a sequence as a valid barcode
BC_MIN_QUAL = 25
UMI_MIN_QUAL = 15

# This is the pattern for 10x indexes
PATTERN = "BBBBBBBBBBBBBBBBUUUUUUUUUU"

# This is the pattern for transcript level genomic barcodes
BC2_LENGTH = 19
BC2_SEED = 'ttttctaaggatcca'

UNKNOWN_STR = "UNK"


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


def link_barcodes(bc_fastq_1, bc_fastq_2, bc_fastq_3, out_file_path, is_zipped=False, allowed_indexes=None,
                  max_index_mismatch=1, bc1_pattern=PATTERN, bc2_seed=BC2_SEED, bc2_len=BC2_LENGTH, bc2_mapfile=None):
    assert len(bc_fastq_1) == len(bc_fastq_2)

    if bc2_mapfile is not None:
        with open(bc2_mapfile, mode="r") as mapfh:
            bc2_map = fasta_to_dict(mapfh, key_type='seq')
    else:
        bc2_map = None

    linker = Link10xBCCounts(bc1_pattern, bc2_len, bc2_seed, is_zipped=is_zipped)
    mp_pool = multiprocessing.Pool(processes=len(bc_fastq_1), maxtasksperchild=100)

    bcs = dict()
    for bc_dict in mp_pool.imap_unordered(linker.parse_fastq_mp, zip(bc_fastq_1, bc_fastq_2, bc_fastq_3)):
        bcs = merge_bcs(bcs, bc_dict)

    output_bcs(out_file_path, bcs, allowed_indexes=allowed_indexes, max_index_mismatch=max_index_mismatch,
               bc2_map=bc2_map)


def output_bcs(out_file_path, bc_dict, allowed_indexes=None, bc2_map=None, max_index_mismatch=1, max_bc_mismatch=1):

    with open(out_file_path, mode="w") as outfh:
        for idx, idx_dict in bc_dict.items():
            if allowed_indexes is not None:
                # See if the index is in the keep list, or if it's within max_mismatch of the keep list
                try:
                    idx = reindex_for_mismatches(idx, allowed_indexes, max_mismatch=max_index_mismatch)
                except IndexError:
                    continue
            else:
                pass

            # If the index is OK, print all of the BC1, BC2, UMI counts for that index
            for bc1, bc2_dict in idx_dict.items():

                # Merge together the UMI counts from barcodes that have an acceptable number of mismatches
                if bc2_map is not None and max_bc_mismatch > 0:
                    bc2_dict = merge_mismatches(bc2_dict, list(bc2_map.keys()), max_mismatch=max_bc_mismatch)

                # Print the Index, Barcode 1, Barcode 2, and UMI counts as a TSV
                for bc2, umi_set in bc2_dict.items():
                    if bc2_map is not None:
                        try:
                            bc2 = bc2_map[bc2]
                        except KeyError:
                            pass

                    print("\t".join([idx, bc1, bc2, str(len(umi_set))]), file=outfh)


def merge_bcs(bc_dict1, bc_dict2):
    for idx, level_1_dict in bc_dict2.items():
        for bc1, level_2_dict in level_1_dict.items():
            for bc2, umi_set in level_2_dict.items():
                try:
                    bc_dict1[idx][bc1][bc2].union(umi_set)
                except KeyError:
                    nest_dict(bc_dict1, idx, bc1, bc2)
                    bc_dict1[idx][bc1][bc2] = umi_set
    return bc_dict1


def merge_mismatches(sequence_dict, allowed_sequences, max_mismatch=1):
    for seq in list(sequence_dict.keys()):
        if seq in allowed_sequences:
            pass
        else:
            closest = search_list_degenerate(seq, allowed_sequences, max_mismatch=max_mismatch)
            if closest is not None:
                merge = sequence_dict.pop(seq)
                try:
                    sequence_dict[closest] = sequence_dict[closest].union(merge)
                except KeyError:
                    sequence_dict[closest] = merge
            else:
                pass

    return sequence_dict


def reindex_for_mismatches(idx, allowed_idx, max_mismatch=1):
    if idx in allowed_idx:
        return idx
    for kidx in allowed_idx:
        if sequence_mismatch_degenerate(idx, kidx) <= max_mismatch:
            return kidx
    raise IndexError


if __name__ == '__main__':
    main()
