from __future__ import unicode_literals

import sys
import string
import itertools

DNA_EQUALITY = {"A": ["A", "N", "R", "M", "W", "V", "D", "H"],
                "T": ["T", "N", "Y", "K", "W", "B", "D", "H"],
                "G": ["G", "N", "R", "K", "K", "B", "D", "V"],
                "C": ["C", "N", "Y", "M", "K", "B", "V", "H"],
                "N": ["A", "T", "G", "C", "N", "R", "Y", "M", "K", "W", "V", "D", "H", "B"],
                "R": ["R", "A", "G"],
                "Y": ["Y", "C", "T"],
                "M": ["M", "A", "C"],
                "S": ["G", "C", "S"],
                "W": ["W", "A", "T"],
                "K": ["G", "T", "K"],
                "V": ["V", "A", "G", "C"],
                "D": ["D", "A", "G", "T"],
                "H": ["H", "A", "T", "C"],
                "B": ["B", "T", "G", "C"]}

SIMPLE_REGEX = {"A": ["A"],
                "T": ["T"],
                "G": ["G"],
                "C": ["C"],
                "N": ["A", "T", "G", "C"],
                "R": ["A", "G"],
                "Y": ["C", "T"],
                "M": ["A", "C"],
                "S": ["G", "C"],
                "W": ["A", "T"],
                "K": ["G", "T"],
                "V": ["A", "G", "C"],
                "D": ["A", "G", "T"],
                "H": ["A", "T", "C"],
                "B": ["T", "G", "C"]}

HAMMING_1 = {"A": ["T", "G", "C"],
             "C": ["A", "G", "T"],
             "G": ["A", "T", "C"],
             "T": ["A", "G", "C"]}

#TODO: Figure out a non hacky way to make this work. Or dont.
if sys.version_info[0] >= 3:
    RC_TABLE = str.maketrans(dict(A="T", G="C", T="A", C="G"))
else:
    RC_TABLE = string.maketrans("ATGC", "TACG")


def search_list_degenerate(str1, list_of_strings, alphabet=DNA_EQUALITY, max_mismatch=0):
    mm_list = []

    for str2 in list_of_strings:
        mm_list.append(sequence_mismatch_degenerate(str1, str2, alphabet=alphabet))

    if min(mm_list) <= max_mismatch:
        return list_of_strings[mm_list.index(min(mm_list))]
    else:
        return None


def sequence_identical_degenerate(str1, str2, alphabet=DNA_EQUALITY):
    str1, str2 = str1.strip(), str2.strip()
    if len(str1) != len(str2):
        return False

    if sequence_mismatch_degenerate(str1, str2, alphabet=alphabet) == 0:
        return True
    else:
        return False


def sequence_mismatch_degenerate(str1, str2, alphabet=DNA_EQUALITY):
    str1 = str1.strip().upper()
    str2 = str2.strip().upper()

    mm = 0

    for i, j in itertools.zip_longest(str1, str2):
        if i is None or j is None:
            mm += 1
        else:
            if i.upper() in alphabet[j.upper()]:
                continue
            else:
                mm += 1

    return mm


def make_regex_str(str1, alphabet=SIMPLE_REGEX, num_mismatch=0):
    regex = ""
    for chr1 in str1.upper():
        regex += "[" + "".join(alphabet[chr1]) + "]"

    if num_mismatch > 0:
        regex = "(" + regex + "){s<=" + str(num_mismatch) + "}"

    return regex


def convert_pattern(pattern):
    pattern = pattern.upper()
    pattern_idx = {}
    for c in string.ascii_uppercase:
        idx = find_all_char(pattern, c)
        if len(idx) > 0:
            pattern_idx[c] = (min(idx), max(idx) + 1)
    return pattern_idx


def find_all_char(sstr, character):
    idx = []
    for pos, chr in enumerate(sstr):
        if character == chr:
            idx.append(pos)
    return idx


def rc_str(str1):
    return str1.upper().translate(RC_TABLE)[::-1]


def make_merge_map(bc_list):

    merge_map = {}

    def hamming_dist_1_gen(seq_str):
        for i, c in enumerate(seq_str.upper()):
            for n in HAMMING_1[c]:
                yield seq_str[0:i] + n + seq_str[i + 1:]

    for bc in bc_list:
        bc = bc.upper()
        merge_map[bc] = bc
        for degen in hamming_dist_1_gen(bc):
            merge_map[degen] = bc

    for bc in bc_list:
        merge_map[bc] = bc

    return merge_map



