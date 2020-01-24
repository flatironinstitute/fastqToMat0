import pandas as pd

STEINMETZ_COLUMNS = ["chr", "strand", "t5", "t3", "ypd", "gal", "type", "class", "name"]
STEINMETZ_ID = 'name'
STEINMETZ_STRAND = 'strand'
STEINMETZ_COUNT = ["ypd", "gal"]
STEINMETZ_MIN_COUNT = 20


def load_steinmetz(fh):
    """
    Load the supplemental transcript files from Steinmetz 2013
    :param fh: file handle
    :return: pd.DataFrame
    """
    df = []
    for i, line in enumerate(fh):
        line_arr = line.strip().split()
        if i == 0:
            continue
        temp_data = line_arr[0:6]
        temp_id = line_arr.pop()
        temp_data.extend([" ".join(line_arr[6:]), line_arr[-1], temp_id])
        df.append(temp_data)
    return pd.DataFrame(df, columns=STEINMETZ_COLUMNS)


def process_steinmetz_data(stein_df, minimum_count=STEINMETZ_MIN_COUNT):
    """
    Process the dataframe from Steinmetz 2013 into a dict of (Start, Stop) tuples, keyed by gene
    :param stein_df: pd.DataFrame
    :param minimum_count: int
    :return transcript: {Gene:(Start, Stop)}
    """
    transcript = dict()
    genes = stein_df.loc[stein_df['class'] == 'ORF'][STEINMETZ_ID].unique().tolist()
    for gene in genes:
        gdata = stein_df.loc[stein_df[STEINMETZ_ID] == gene]
        try:
            strand = gdata.loc[:, STEINMETZ_STRAND].tolist()[0]
            counts = gdata.loc[:, STEINMETZ_COUNT].astype(int).sum(axis=1) > minimum_count
        except IndexError:
            continue
        if strand == "+":
            left, right = 't5', 't3'
        elif strand == "-":
            left, right = 't3', 't5'
        else:
            raise ValueError
        if counts.sum() == 0:
            continue
        gdata = gdata.loc[counts]
        transcript[gene] = (gdata[left].astype(int).min(), gdata[right].astype(int).max())
    return transcript
