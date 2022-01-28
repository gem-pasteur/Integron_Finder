from typing import Tuple



def description_parser(seq_header:str) -> Tuple[int,int,int,int]:
    """
    """
    id_, start, stop, strand, *_ = seq_header.split(" # ")
    start = int(start)
    stop = int(stop)
    strand = int(strand)
    return (id_, start, stop, strand)


