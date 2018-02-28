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


read_multi_prot_fasta = make_multi_fasta_reader(Seq.IUPAC.protein)


class FastaIterator(object):
    """
    :param path:the path to the fasta file
    :return: the sequence parsed
    :rtype: :class:`Bio.SeqRecord.SeqRecord` object
    """
    def __init__(self, alphabet):
        self.alphabet = alphabet
        self.name = None
        self.seq_index = None

    def __iter__(self):
        for id_ in self.seq_index.keys():
            seq = self.seq_index[id_]
            seq.name = self.name
            yield seq

    def __len__(self):
        return len(self.seq_index)

    def __call__(self, path):
        self.name = get_name_from_path(path)
        self.seq_index = SeqIO.index(path, "fasta", alphabet=self.alphabet)
        return self


read_multi_dna_fasta = FastaIterator(Seq.IUPAC.unambiguous_dna)

read_multi_prot_fasta = FastaIterator(Seq.IUPAC.protein)


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


SeqDesc = namedtuple('ProtDesc', ('id', 'strand', 'start', 'stop'))


def gembase_parser(description):
    desc = description.split(" ")
    id_, strand, start, stop = desc[:2] + desc[4:6]
    strand = 1 if desc[1] == "D" else -1
    start = int(start)
    stop = int(stop)
    return SeqDesc(id_, strand, start, stop)


def non_gembase_parser(description):
    id_, start, stop, strand, _ = description.split(" # ")
    start = int(start)
    stop = int(stop)
    strand = int(strand)
    return SeqDesc(id_, strand, start, stop)