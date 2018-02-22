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


def model_len(path):
    if not os.path.exists(path):
        raise RuntimeError("Path to model_attc '{}' doe snot exists".format(path))
    with open(path) as model_file:
        for line in model_file:
            if line.startswith('CLEN'):
                model_len = int(line.split()[1])
                return model_len
        raise RuntimeError("CLEN not found in '{}', maybe it's not infernal model file".format(path))


def get_name_from_path(path):
    return os.path.splitext(os.path.split(path)[1])[0]
