import pandas as pd


class toto:

    def __init__(self, p, id):
        df = pd.read_table(p, header=None, names=['start', 'end', 'strand', 'type', 'seq_id'], usecols=(0, 1, 2, 3, 4), dtype={'start': 'int', 'end': 'int', 'strand': 'str', 'type': 'str', 'seq_id': 'str'})
        specie, date, strain, contig = id.split('.')
        pattern = '{}\.{}\.{}\.\w?{}'.format(specie, date, strain, contig)
        print("pattern =", pattern)
        genome_info = df[df['seq_id'].str.contains(pattern, regex=True)]
        print(genome_info.shape)
        self._info = genome_info[genome_info['type'] == 'CDS']
        print(self._info.shape)


    def __iter__(self):
        return (id for id in self._info['seq_id'])


db = toto('../LSTINFO/ACBA.0917.00019.lst', 'ACBA.0917.00019.0001' )
