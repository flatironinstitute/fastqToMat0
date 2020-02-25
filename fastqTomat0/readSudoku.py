from fastqTomat0.lib.barcode_indexer import LinkSudokuBC

import pandas as pd
import textdistance

I5_COL = "i5"
I7_COL = "i7"
GROUP_COL = "Group"
LOC_COL = "Location"
BC_COL = "Barcode"
COUNT_COL = "Count"


def process_sudoku_reads(output_file, fastqr1, fastqr2, fastqi1, fastqi2, bc_seed, bc_length, index_table_file,
                         is_zipped=True, threshold_count=5, merge_dist=1):

    assert len(fastqr1) == len(fastqr2)
    assert len(fastqr1) == len(fastqi1)
    assert len(fastqr1) == len(fastqi2)

    print("Loading index table")

    idx_table = pd.read_csv(index_table_file, index_col=None, sep="\t")
    allowed_i5 = idx_table[I5_COL].tolist()
    allowed_i7 = idx_table[I7_COL].tolist()

    print("Loaded {n} pools with {i5} i5 indices and {i7} i7 indices".format(n=idx_table.shape[0],
                                                                             i5=len(allowed_i5),
                                                                             i7=len(allowed_i7)))

    linker = LinkSudokuBC(bc_seed, bc_length, is_zipped=is_zipped)
    bcs = linker.parse_fastq_mp(fastqi1, fastqi2, fastqr1)

    print("Barcode extraction complete")

    # Only keep indices that are in the table
    for i5 in list(bcs.keys()):
        if i5 not in allowed_i5:
            del bcs[i5]
        else:
            for i7 in list(bcs[i5].keys()):
                if i7 not in allowed_i7:
                    del bcs[i5][i7]

    print("Index restriction to whitelist complete")

    hamming = textdistance.Hamming()

    outputs = []
    for i5 in list(bcs.keys()):

        i5_ref = bcs[i5]

        for i7 in list(i5_ref.keys()):
            i57_ref = bcs[i5][i7]

            try:
                group = idx_table.loc[(idx_table[I5_COL] == i5) & (idx_table[I7_COL] == i7), GROUP_COL].tolist()[0]
                pool = idx_table.loc[(idx_table[I5_COL] == i5) & (idx_table[I7_COL] == i7), LOC_COL].tolist()[0]
            except IndexError:
                continue

            print("Processing {g} {p}".format(g=group, p=pool))

            # Keep only BCs over a count threshold
            valid_bcs = []
            for bc in list(i57_ref.keys()):
                if i57_ref[bc] > threshold_count:
                    valid_bcs.append(bc)

            # Merge all BCs that are less than merge_dist hamming distance
            bc_merge_toward = set()
            bc_merge_away = set()
            merged_bcs = set()
            for i in range(len(valid_bcs)):
                bc = valid_bcs[i]
                end_check = False
                for opp in valid_bcs[i+1:]:
                    if hamming.distance(bc, opp) < merge_dist:
                        end_check = True
                        if i57_ref[bc] > i57_ref[opp]:
                            bc_merge_toward.add(bc)
                            bc_merge_away.add(opp)
                        else:
                            bc_merge_toward.add(opp)
                            bc_merge_away.add(bc)

                if not end_check:
                    merged_bcs.add(bc)

            merged_bcs = merged_bcs.union(bc_merge_toward.difference(bc_merge_away))

            merged_bcs = list(merged_bcs)
            merged_counts = [i57_ref[bc] for bc in merged_bcs]

            output_df = pd.DataFrame(list(zip(merged_bcs, merged_counts)), columns=[BC_COL, COUNT_COL])
            output_df[GROUP_COL] = group
            output_df[LOC_COL] = pool
            outputs.append(output_df)

    outputs = pd.concat(outputs)
    outputs = outputs.reindex(labels=[GROUP_COL, LOC_COL, BC_COL, COUNT_COL], axis=1)
    outputs.sort_values(by=[GROUP_COL, LOC_COL], inplace=True)
    outputs.reset_index(inplace=True, drop=True)
    outputs.to_csv("raw_" + output_file, sep="\t")

    pivoted = map_barcodes_to_wells(outputs)
    pivoted.to_csv(output_file, sep="\t")

    return outputs, pivoted


def map_barcodes_to_wells(data, group_col=GROUP_COL, bc_col=BC_COL):

    deduplicated =[]
    uniques = None
    for gname, level in data.groupby(group_col):
        level = level.drop_duplicates(subset=bc_col, keep=False)
        deduplicated.append(level)
        uniques = uniques.intersection(set(level[BC_COL])) if uniques is not None else set(level[BC_COL])

    print("Final count: {n} Unique, Mappable Barcodes".format(n=len(uniques)))

    uniques = pd.DataFrame(uniques, columns=[BC_COL])
    deduplicated = pd.concat(deduplicated)
    deduplicated = pd.merge(deduplicated, uniques, how="inner", on=[BC_COL])
    deduplicated = deduplicated.pivot(index=BC_COL, columns=GROUP_COL, values=LOC_COL)

    return deduplicated




