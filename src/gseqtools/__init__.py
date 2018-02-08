from Bio.Seq import Seq
from pysam import FastaFile


class SequenceExtractor(object):
    """

    Extracting sequences from FASTA file by interval objects.

    """

    def __init__(self, path):
        self._fasta = FastaFile(path)

    def get_sequence(self, gi):
        seqs = []
        for x, y in gi.blocks:
            seqs.append(self._fasta.fetch(gi.chrom, x, y))
        seq = Seq("".join(seqs))
        if gi.reverse:
            seq = seq.reverse_complement()
        return seq

    def close(self):
        if not self._fasta.closed:
            self._fasta.close()

    def __del__(self):
        self.close()
