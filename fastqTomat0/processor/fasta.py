def fasta_to_dict(fh, key_type='id', require_unique=True):
    if key_type not in ['id', 'seq']:
        raise ValueError("Key type {} not allowed".format(key_type))

    fasta_dict = dict()
    current_key = None
    current_seq = ""

    for line in fh:
        line = line.strip()

        if len(line) == 0:
            continue
        elif line[0] == "#" or line[0] == "!" or line[0] == "@":
            continue
        elif line[0] == ">":
            if current_key is not None:
                add_to_dict(fasta_dict, current_key, current_seq, require_unique, key_type)

            current_key = line[1:]
            current_seq = ""
        else:
            current_seq += line

    if current_key is not None and current_seq != "":
        add_to_dict(fasta_dict, current_key, current_seq, require_unique, key_type)

    return fasta_dict


def add_to_dict(fasta_dict, key, seq, require_unique=True, key_type='id'):
    if key_type is 'id':
        if require_unique and key in fasta_dict:
            raise KeyError("Sequence ID is not unique")
        fasta_dict[key] = seq
    if key_type is 'seq':
        if require_unique and seq in fasta_dict:
            raise KeyError("Sequence ID is not unique")
        fasta_dict[seq] = key
