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

                if key_type is 'id':
                    if require_unique and current_key in fasta_dict:
                        raise KeyError("Sequence ID is not unique")
                    fasta_dict[current_key] = current_seq
                if key_type is 'seq':
                    if require_unique and current_seq in fasta_dict:
                        raise KeyError("Sequence ID is not unique")
                    fasta_dict[current_seq] = current_key

            current_key = line[1:]
            current_seq = ""
        else:
            current_seq += line

    return fasta_dict