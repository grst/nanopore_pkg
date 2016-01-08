"""
Since skbio and Biopython are overkill and slightly to complicated most of the time
I came up with this really simple fasta-io class.
"""
from itertools import groupby


class FastaReader(object):
    def __init__(self, file):
        if not hasattr(file, 'read'):
            self.file = open(file, 'r')
        else:
            self.file = file

    def get_entries(self):
        """
        Get the next Entry from the fasta file.

        Returns: Generator, which yields (header, sequence) tuples

        """
        for isheader, group in groupby(self.file, lambda line: line[0] == ">"):
            if isheader:
                header = next(group)[1:]
            else:
                seq = "".join(line.strip() for line in group)
                yield header, seq

    def close(self):
        self.file.close()


class FastaWriter(object):
    """
    Very simple fasta file format writer.
    """
    SPLIT = 80

    def __init__(self, file):
        if not hasattr(file, 'write'):
            self.file = open(file, 'w')
        else:
            self.file = file

    def write_entry(self, header, sequence):
        """
        Write Entry to File

        Args:
            header: >sequence_header
            sequence: ACTGATT...
        """
        sequence = [sequence[i:i+self.SPLIT] for i in range(0, len(sequence), self.SPLIT)]
        self.file.write(">{0}\n".format(header))
        for s in sequence:
            self.file.write(s + "\n")

    def close(self):
        self.file.close()

