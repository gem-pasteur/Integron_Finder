import os
from Bio import Seq
from Bio import SeqIO

def read_fasta(path):
    """
    :param path:the path to the fasta file
    :return: the sequence parsed
    :rtype: :class:`Bio.SeqRecord.SeqRecord` object
    """
    seq = SeqIO.read(path, "fasta",  alphabet=Seq.IUPAC.unambiguous_dna)
    seq.name = os.path.splitext((os.path.split(path)[1]))[0]
    return seq
