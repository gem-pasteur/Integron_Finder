from . import register, SeqDesc


@register
def prodigal_base(fasta_header):
    id_, start, stop, strand, _ = fasta_header.split(" # ")
    start = int(start)
    stop = int(stop)
    strand = int(strand)
    return SeqDesc(id_, strand, start, stop)


@register
def gembase_old(fasta_header):
    desc = fasta_header.split(" ")
    id_, strand, start, stop = desc[:2] + desc[4:6]
    strand = 1 if desc[1] == "D" else -1
    start = int(start)
    stop = int(stop)
    return SeqDesc(id_, strand, start, stop)


@register
def gembase(fasta_header):
    desc = fasta_header.split('.')
    id_, strand, start, stop = desc[:2] + desc[4:6]
    strand = 1 if desc[1] == "D" else -1
    start = int(start)
    stop = int(stop)
    return SeqDesc(id_, strand, start, stop)

