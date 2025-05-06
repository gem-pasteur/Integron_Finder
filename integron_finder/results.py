####################################################################################
# Integron_Finder - Integron Finder aims at detecting integrons in DNA sequences   #
# by finding particular features of the integron:                                  #
#   - the attC sites                                                               #
#   - the integrase                                                                #
#   - and when possible attI site and promoters.                                   #
#                                                                                  #
# Authors: Jean Cury, Bertrand Neron, Eduardo PC Rocha                             #
# Copyright (c) 2015 - 2025  Institut Pasteur, Paris and CNRS.                     #
# See the COPYRIGHT file for details                                               #
#                                                                                  #
# integron_finder is free software: you can redistribute it and/or modify          #
# it under the terms of the GNU General Public License as published by             #
# the Free Software Foundation, either version 3 of the License, or                #
# (at your option) any later version.                                              #
#                                                                                  #
# integron_finder is distributed in the hope that it will be useful,               #
# but WITHOUT ANY WARRANTY; without even the implied warranty of                   #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                    #
# GNU General Public License for more details.                                     #
#                                                                                  #
# You should have received a copy of the GNU General Public License                #
# along with this program (COPYING file).                                          #
# If not, see <http://www.gnu.org/licenses/>.                                      #
####################################################################################

import colorlog
import pandas as pd
from pandas.errors import EmptyDataError

_log = colorlog.getLogger(__name__)

"""
utilities to manage results
"""


def integrons_report(integrons):
    """

    :param integrons: list of integrons used to generate a report
    :type integrons: list of :class:`integron_finder.integron.Integron` object.
    :return: a report off all integrons from a replicon
    :rtype: :class:`pandas.DataFrame` object.
            this datafame have following columns:

            "ID_integron", "ID_replicon", "element",
            "pos_beg", "pos_end", "strand", "evalue",
            "type_elt", "annotation", "model",
            "type", "default", "distance_2attC", "considered_topology"
    """
    integrons_describe = pd.concat([i.describe() for i in integrons])
    dic_id = {id_: "{:02}".format(j) for j, id_ in
              enumerate(integrons_describe.sort_values("pos_beg").ID_integron.unique(), 1)}
    integrons_describe.ID_integron = ["integron_" + dic_id[id_] for id_ in integrons_describe.ID_integron]
    integrons_describe = integrons_describe[["ID_integron", "ID_replicon", "element",
                                             "pos_beg", "pos_end", "strand", "evalue",
                                             "type_elt", "annotation", "model",
                                             "type", "default", "distance_2attC", "considered_topology"]]
    integrons_describe['evalue'] = integrons_describe.evalue.astype(float)
    integrons_describe.sort_values(["ID_integron", "pos_beg", "evalue"], inplace=True)
    integrons_describe.index = list(range(len(integrons_describe)))
    return integrons_describe


def merge_results(*results_file):
    """

    :param results_file: The path of the files to merge.
                         The files can be parsed by pandas as DataFrame
                         and have the same columns.
                         It is used to merge the integrons files (.integrons)
                         or summary files (.summary) from different replicons.
    :type results_file: str
    :return: all results aggregated in one :class:`pandas.DataFrame` object.
             if there is no results to merge, return an empty DataFrame.
    :rtype: a :class:`pandas.DataFrame` object.
    """
    _log.debug("merging_results " + ' '.join(results_file))
    dtype = {"ID_integron": 'str',
             "ID_replicon": 'str',
             "element": 'str',
             "pos_beg": 'int',
             "pos_end": 'int',
             "strand": 'int',
             "evalue": 'float',
             "type_elt": 'str',
             "annotation": 'str',
             "model": 'str',
             "type": 'str',
             "default": 'str',
             "distance_2attC": 'float',
             "considered_topology": 'str'
             }
    all_res = []
    for one_result in results_file:
        try:
            res = pd.read_csv(one_result, sep="\t", comment='#', dtype=dtype)
        except EmptyDataError:
            continue
        all_res.append(res)
    if all_res:
        agg_results = pd.concat(all_res, sort=False)
        if agg_results.shape[1] == 6 or agg_results.shape[1] == 4:
            # it's a summary file
            # 4 fields old Integron finder
            # 6 fields new Integron finder (2 new cols 'topology' and 'size')
            agg_results = agg_results.set_index('ID_replicon')
    else:
        agg_results = pd.DataFrame(columns=['ID_integron', 'ID_replicon', 'element',
                                            'pos_beg', 'pos_end', 'strand', 'evalue',
                                            'type_elt', 'annotation', 'model', 'type', 'default',
                                            'distance_2attC', 'considered_topology'])
    return agg_results


def summary(result):
    """
    Create a summary of an integron report.
    Count the number of 'CALIN', 'In0' or 'complete' for each replicon.

    :param result: the integron to summarize
    :return: a :class:`pandas.DataFrame` object.
             with columns 'ID_replicon', 'ID_integron', 'complete', 'In0', 'CALIN'
    """
    dropped = result.drop_duplicates(subset=['ID_integron', 'ID_replicon'])
    summary = pd.crosstab([dropped.ID_replicon, dropped.ID_integron], dropped.type)
    empty_cols = set(['CALIN', 'complete', 'In0']).difference(summary.columns)
    for empty_col in empty_cols:
        summary[empty_col] = 0
    summary = summary.reset_index()
    summary.columns.name = None
    aggregate = summary.groupby('ID_replicon').sum()

    aggregate = aggregate.reindex(['CALIN', 'complete', 'In0'], axis=1)
    return aggregate


def filter_calin(result, threshold=2):
    """
    filter integron report, remove 'CALIN' integron where number of attC sites is lower than threshold.

    :param result: the output of :func:`integrons_report`
    :type result: :class:`pandas.dataFrame` object
    :param int threshold: the integron CALIN with less attc site than *threshold* are removed
    :return: filtered integron report
    :rtype: :class:`pandas.dataFrame` object
    """
    #     In [87]: %%timeit
    #     ...: d = pd.read_csv('Results_Integron_Finder_multi_fasta/multi_fasta.integrons', sep='\t')
    #     ...: d.set_index(['ID_integron', 'ID_replicon'], inplace=True)
    #     ...: idx = d[(d.type_elt=='attC') & (d.type=='CALIN')].groupby(level=['ID_replicon', 'ID_integron']).filter(lambda x: x['type'].size<=1).index
    #     ...: d.loc[d.index.difference(idx)].drop_duplicates()
    #     ...:
    # 10 loops, best of 3: 30.4 ms per loop

    # In [88]: %%timeit
    #     ...: d = pd.read_csv('Results_Integron_Finder_multi_fasta/multi_fasta.integrons', sep='\t')
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
