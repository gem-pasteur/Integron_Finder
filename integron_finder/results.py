import pandas as pd
from pandas.io.common import EmptyDataError


def merge_results(*results_file):
    """

    :param results_file: The path of the files to merge.
                         The files can be parsed by pandas as DataFrame
                         and have the same columns.
                         It is used to merge the integrons files (.integrons)
                         or summary files (.summary) from different replicons.
    :type results_file: str
    :return: all results aggregated in one :class:`pandas.DataFrame` object
    :rtype: a :class:`pandas.DataFrame` object.
    """
    all_res = []
    for one_result in results_file:
        try:
            res = pd.read_table(one_result, sep="\t", comment='#')
        except EmptyDataError:
            continue
        all_res.append(res)
    agg_results = pd.concat(all_res)
    return agg_results


def summary(result):
    dropped = result.drop_duplicates(subset=['ID_integron', 'ID_replicon'])
    summary = pd.crosstab([dropped.ID_replicon, dropped.ID_integron], dropped.type)
    empty_cols = set(['CALIN', 'complete', 'In0']).difference(summary.columns)
    for empty_col in empty_cols:
        summary[empty_col] = 0
    summary = summary.reset_index()
    summary.columns.name = None
    return summary[['ID_replicon', 'ID_integron', 'complete', 'In0', 'CALIN']]


def filter_calin(result, threshold=2):
    #     In [87]: %%timeit
    #     ...: d = pd.read_table('Results_Integron_Finder_multi_fasta/multi_fasta.integrons', sep='\t')
    #     ...: d.set_index(['ID_integron', 'ID_replicon'], inplace=True)
    #     ...: idx = d[(d.type_elt=='attC') & (d.type=='CALIN')].groupby(level=['ID_replicon', 'ID_integron']).filter(lambda x: x['type'].size<=1).index
    #     ...: d.loc[d.index.difference(idx)].drop_duplicates()
    #     ...:
    # 10 loops, best of 3: 30.4 ms per loop

    # In [88]: %%timeit
    #     ...: d = pd.read_table('Results_Integron_Finder_multi_fasta/multi_fasta.integrons', sep='\t')
    #     ...: d['ID'] = d.ID_replicon + '_' + d.ID_integron
    #     ...: idx = d[(d.type_elt=='attC') & (d.type=='CALIN')].groupby('ID').filter(lambda x: x['type'].size<2).ID
    #     ...: d[~d.ID.isin(idx)]
    #     ...:
    # 100 loops, best of 3: 13.8 ms per loop

    idx = result[(result.type_elt == 'attC') & (result.type == 'CALIN')].\
        groupby('ID_integron').\
        filter(lambda x: x['type'].size < threshold).ID_integron
    filtered = result[~result.ID_integron.isin(idx)]
    return filtered
