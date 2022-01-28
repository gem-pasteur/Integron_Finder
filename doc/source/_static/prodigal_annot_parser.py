from typing import Tuple


def description_parser(seq_header: str) -> Tuple[int, int, int, int]:
    """
    Extract description features from protein sequence fasta header

    :param str seq_header: the fasta header of a sequence (without '>' char)
    :return: features
                - sequence id (string)
                - start the beginning of the protein (positive integer)
                - stop the end position of the protein (positive integer)
                - strand 1 if prot is coded on direct strand or -1 if it's on reverse
    :rtypes: tuple (str, int, int, 1/-1)
    """
    id_, start, stop, strand, *_ = seq_header.split(" # ")
    start = int(start)
    stop = int(stop)
    strand = int(strand)
    return id_, start, stop, strand


