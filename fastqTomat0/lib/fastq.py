# Converts a quality ASCII string to a list of qualities
# This is the 33-offset illumina quality scoring
def convert_qual_illumina(qstr):
    qual = []
    for ch in qstr:
        qual.append(ord(ch) - 33)
    return qual

# Dict to map strings to quality score functions
PHRED_SCORE = {"illumina": convert_qual_illumina}

# Wrapper for fastqProcessor that takes a single file handle
# Returns a tuple of (control, sequence, quality) instead of a list of tuples
def fastq_gen(fh, phred_type='illumina'):
    for rec in fastqProcessor(phred_type=phred_type).fastq_gen(fh):
        yield rec[0]

# The fastqProcessor class takes an arbitrary number of linked reads and yields a list of tuples for each sequence
# The list is in the same order as the files that were given as arguments to fastq_gen
class fastqProcessor:

    phred = None
    verify_ids = True

    def __init__(self, phred_type='illumina', verify_ids=True):
        """
        Set the phred scoring function and verify_ids flags
        """

        try:
            phred = PHRED_SCORE[phred_type]
        except KeyError:
            raise ValueError("Score type {} unknown".format(phred_type))
        self.phred = phred
        self.verify = verify_ids

    def fastq_gen(self, *fhs):
        """
        Generator that yields a list of N tuples (sequence ID, sequence, quality) where N is the number of file handles
        given as arguments
        :param fhs: file handle
            File handles to the fastQ files
        :yield record: list[(str, str, list(int)]
            A list of N tuples where each tuple is a linked read
        """

        i = 0
        while True:
            i += 1
            try:
                record = []
                cid = None
                for fh in fhs:
                    (c, s, q) = self.fastq_process_file(fh, self.phred)

                    # If verify_ids was set, assert that the sequence IDs all match
                    if self.verify and cid is not None:
                        try:
                            assert self.extract_control_id(c) == cid
                        except AssertionError:
                            print("ID mismatch (Record {i}): {cid} != {oid}".format(i=i,
                                                                                    cid=cid,
                                                                                    oid=self.extract_control_id(c)))
                            raise
                    cid = self.extract_control_id(c)

                    record.append((c, s, q))
                yield record
            except StopIteration:
                break

    @staticmethod
    def fastq_process_file(fh, phred):
        """
        Read the next record and return it as a tuple of (context, sequence, quality)
        """
        cont, seq, qual = None, None, None
        line_id = -1
        for line in fh:
            line = line.strip()
            try:
                if line[0] == "+":
                    continue
                elif line[0] == "@":
                    line_id = 0
                    cont = line
                elif line_id == 0:
                    line_id = 1
                    seq = line
                elif line_id == 1:
                    line_id = 2
                    qual = phred(line)
                    return cont, seq, qual
            except IndexError:
                continue
        raise StopIteration

    @staticmethod
    def extract_control_id(con):
        try:
            return con.strip().split()[0]
        except IndexError:
            return None
