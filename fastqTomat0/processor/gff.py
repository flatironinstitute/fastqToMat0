from __future__ import print_function

import pandas as pd
import numpy as np
import csv

GFF_START = "start"
GFF_STOP = "end"
GFF_FEATURE = "feature"
GFF_ATTRIBUTE = "attribute"
GFF_COLUMNS = ["seqname", "source", GFF_FEATURE, GFF_START, GFF_STOP, "score", "strand", "frame", GFF_ATTRIBUTE]
ATTRIBUTE_COLUMN = 8
ATTRIBUTE_EXTRACTION_COLUMNS = ["gene_id"]
PARENT_TAG = 'Parent'
ID_TAG = 'ID'
GROUPNAME_TAG = 'locus_tag'

MODIFY_FEATURES = ["gene", "mRNA", "exon", "transcript"]
MODIFY_SEARCH_KEY = 'gene_id'
MODIFY_GROUP_KEY = 'gene_id'


def load_gff(fh, extract_cols=ATTRIBUTE_EXTRACTION_COLUMNS, gff_version=2):
    """
    Load a GFF file in as a data frame
    :param fh: file handle
    :return df: pd.DataFrame
    """
    if gff_version == 2:
        att_proc = gff2_attributes
    elif gff_version == 3:
        att_proc = gff3_attributes
    else:
        raise ValueError("Only GFF2 and GFF3 are supported")

    df = []
    headers = []

    for i, line in enumerate(fh):
        line = line.strip()
        if line.startswith("#"):
            headers.append(line)
            continue
        if len(line) == 0:
            continue
        line_arr = line.split("\t")
        try:
            atts = att_proc(line_arr[ATTRIBUTE_COLUMN])
        except IndexError:
            print("Malformed Line {i}".format(i=i))
            print(line)
            raise
        for i in range(len(extract_cols)):
            try:
                val = atts[extract_cols[i]]
                if len(val) == 1:
                    val = val[0]
                line_arr.append(val)
            except KeyError:
                line_arr.append(None)
        df.append(line_arr)

    return pd.DataFrame(df, columns=GFF_COLUMNS + extract_cols), headers


def write_gff(fh, df, headers=None):
    if headers is not None:
        for line in headers:
            print(line, file=fh)
    df.loc[:, GFF_COLUMNS].to_csv(fh, sep="\t", index=False, header=False, quoting=csv.QUOTE_NONE, doublequote=False)


def gff3_attributes(att):
    """
    Process the attributes column into a dict of key - list(values) pairs
    :param att: str
    :return attribute_dict: dict
    """
    att_arr = att.split(";")  # Break K=V;K=V;K=V into key-value pairs
    if len(att_arr) == 0:
        return None
    attribute_dict = dict()
    for a in att_arr:
        k, v = a.split("=")  # Break Key = Value pair
        attribute_dict[k] = list()
        for sub_v in v.split(","):
            if ":" in sub_v:
                sub_sub_v = sub_v.split(":")
                attribute_dict[k].append((sub_sub_v[0], sub_sub_v[1]))
            else:
                attribute_dict[k].append(sub_v)
    return attribute_dict


def gff2_attributes(att):
    att_arr = att.split(";")  # Break K=V;K=V;K=V into key-value pairs
    if len(att_arr) == 0:
        return None
    attribute_dict = dict()
    for a in att_arr:
        if len(a) < 2:
            continue
        k, v = a.strip().split(" ")  # Break Key = Value pair
        v = v.strip("\" \n")
        attribute_dict[k] = v
    return attribute_dict


def process_gff_into_groups(df):
    parent_dict = dict()
    df_copy = df[[PARENT_TAG, ID_TAG, GROUPNAME_TAG]].copy()
    group = pd.Series(np.repeat(None, df_copy.shape[0]), index=df.index)
    group_name = pd.Series(np.repeat(None, df_copy.shape[0]), index=df.index)
    for i in range(df_copy.shape[0]):
        record = df_copy.iloc[i]
        cid, cparent = record[ID_TAG], record[PARENT_TAG]
        parent_dict[cid] = cparent
    for i in range(df_copy.shape[0]):
        record = df_copy.iloc[i]
        cid, cparent, cgroup = record[ID_TAG], record[PARENT_TAG], record[GROUPNAME_TAG]
        if cparent is not None:
            while cparent is not None:
                parid = cparent
                cparent = parent_dict[cparent]
            group.iloc[i] = parid
        else:
            group.iloc[i] = cid
        if cgroup is not None:
            group_name.iloc[i] = cgroup
    return df.assign(group=group, group_name=group_name)


def modify_gff_locations(df, new_positions, features=MODIFY_FEATURES, strict=True):
    for k, (left, right) in new_positions.items():
        # Get the group ID associated with the gene
        try:
            group = df.loc[df[MODIFY_SEARCH_KEY] == k][MODIFY_GROUP_KEY].tolist()[0]
        except IndexError:
            print("Unable to locate gene {} in GFF".format(k))
            continue
        # Select the group matching the ID from the entire data set and restrict to selected features
        dfg = df.loc[df[MODIFY_GROUP_KEY] == group]
        dfg = dfg.loc[dfg[GFF_FEATURE].isin(features)]
        # Get the original transcript start and stop
        change_df_locations(dfg, left, right)
        df.update(dfg)


def change_df_locations(df, left, right, change_att_string=False, only_make_longer=True):
    oleft, oright = df[GFF_START].astype(int).min(), df[GFF_STOP].astype(int).max()
    if only_make_longer and left < oleft:
        df.loc[df.loc[:, GFF_START].astype(int) == oleft, GFF_START] = left
    if only_make_longer and right > oright:
        df.loc[df.loc[:, GFF_STOP].astype(int) == oright, GFF_STOP] = right
    if change_att_string:
        df.loc[:, GFF_ATTRIBUTE] = df.loc[:, GFF_ATTRIBUTE].apply(lambda x: x.replace(str(oleft), str(left)))
        df.loc[:, GFF_ATTRIBUTE] = df.loc[:, GFF_ATTRIBUTE].apply(lambda x: x.replace(str(oright), str(right)))
