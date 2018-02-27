import os
from collections import namedtuple
from Bio import Seq
from Bio import SeqIO


def make_single_fasta_reader(alphabet):

    def read_fasta(path):
        """
        :param path:the path to the fasta file
        :return: the sequence parsed
        :rtype: :class:`Bio.SeqRecord.SeqRecord` object
        """
        seq = SeqIO.read(path, "fasta",  alphabet=alphabet)
        seq.name = get_name_from_path(path)
        return seq
    return read_fasta


read_single_dna_fasta = make_single_fasta_reader(Seq.IUPAC.unambiguous_dna)

read_single_prot_fasta = make_single_fasta_reader(Seq.IUPAC.protein)


def make_multi_fasta_reader(alphabet):

    def fasta_iterator(path):
        """
        :param path:the path to the fasta file
        :return: the sequence parsed
        :rtype: :class:`Bio.SeqRecord.SeqRecord` object
        """
        name = get_name_from_path(path)
        seq_it = SeqIO.parse(path, "fasta",  alphabet=alphabet)
        for seq in seq_it:
            seq.name = name
            yield seq

    return fasta_iterator


read_multi_dna_fasta = make_multi_fasta_reader(Seq.IUPAC.unambiguous_dna)

read_multi_prot_fasta = make_multi_fasta_reader(Seq.IUPAC.protein)


def model_len(path):
    if not os.path.exists(path):
        raise IOError("Path to model_attc '{}' does not exists".format(path))
    with open(path) as model_file:
        for line in model_file:
            if line.startswith('CLEN'):
                model_length = int(line.split()[1])
                return model_length
        raise RuntimeError("CLEN not found in '{}', maybe it's not infernal model file".format(path))


def get_name_from_path(path):
    return os.path.splitext(os.path.split(path)[1])[0]


ProtDesc = namedtuple('ProtDesc', ('id', 'strand', 'start', 'stop'))


def gembase_parser(description):
    desc = description.split(" ")
    id_, strand, start, stop = desc[:2] + desc[4:6]
    strand = 1 if desc[1] == "D" else -1
    start = int(start)
    stop = int(stop)
    return ProtDesc(id_, strand, start, stop)


def non_gembase_parser(description):
    id_, start, stop, strand, _ = description.split(" # ")
    start = int(start)
    stop = int(stop)
    strand = int(strand)
    return ProtDesc(id_, strand, start, stop)