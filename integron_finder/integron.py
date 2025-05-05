# -*- coding: utf-8 -*-

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

import os
import colorlog
import numpy as np
import pandas as pd

from matplotlib import use as m_use
m_use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors

from Bio import Seq
from Bio import SeqIO
from Bio import motifs

from .hmm import read_hmm
from .infernal import read_infernal
from .attc import search_attc

_log = colorlog.getLogger(__name__)


def find_integron(replicon, prot_db, intI_file, phageI_file, cfg, attc_file=None, attc=None):
    """
    Function that looks for integrons given rules :
        * presence of intI
        * presence of attC
        * d(intI-attC) <= 4 kb
        * d(attC-attC) <= 4 kb

    It returns the list of all integrons, be they complete or not.
    found in attC files + integrases file which are formatted as follow :
    intI_file : Accession_number    ID_prot    strand    pos_beg    pos_end    evalue
    attc_file : Accession_number    attC    cm_debut    cm_fin    pos_beg    pos_end    sens    evalue

    :param replicon: the name of the replicon
    :type replicon: :class:`Bio.Seq.SeqRecord` object
    :param prot_db: the protein database corresponding to the replicon translation
    :type prot_db: a :class:`integron_finder.prot_db.ProteinDB` object.
    :param str intI_file: the output of hmmsearch with the integrase model
    :param str phageI_file: the output of hmmsearch with the phage model
    :param cfg: configuration
    :param attc_file: the output of cmsearch or the result of parsing of this file by read_infernal
    :type attc_file: path to cmsearch output or :class:`pd.Dataframe`
    :param attc: the result of parsing of cmsearch output file by read_infernal.
                 attc is used when search with local_max
                 attc_file when "normal search
                 attc and attc_file are mutually exclusive parameter
    :type attc: :class:`pd.Dataframe` object
    :type cfg: a :class:`integron_finder.config.Config` object
    :returns: list of all integrons, be they complete or not
    :retype: list of :class:`Integron` object
    """
    assert not(attc is not None and attc_file is not None), "attc and attc_file are mutually exclusive parameter"
    assert not((attc is None) and (not attc_file)), "Either attc or attc_file parameter must be provided"
    if not cfg.no_proteins:
        intI = read_hmm(replicon.id, prot_db, intI_file, cfg)
        intI.sort_values(["Accession_number", "pos_beg", "evalue"], inplace=True)

        phageI = read_hmm(replicon.id, prot_db, phageI_file, cfg)
        phageI.sort_values(["Accession_number", "pos_beg", "evalue"], inplace=True)

        tmp = intI[intI.ID_prot.isin(phageI.ID_prot)].copy()

        if not tmp.empty:
            tmp.loc[:, "query_name"] = "intersection_tyr_intI"

        if cfg.union_integrases:
            intI_ac = intI[intI.ID_prot.isin(tmp.ID_prot) == 0].merge(phageI[phageI.ID_prot.isin(tmp.ID_prot) == 0],
                                                                      how="outer").merge(tmp, how="outer")
            intI_ac.sort_values(['pos_beg', 'pos_end'], ascending=True, inplace=True)
        else:
            intI_ac = tmp
    else:
        intI_ac = pd.DataFrame(columns=["Accession_number", "query_name", "ID_query",
                                        "ID_prot", "strand", "pos_beg", "pos_end",
                                        "evalue", "hmmfrom", "hmmto", "alifrom",
                                        "alito", "len_profile"])
    if attc is not None:
        # attc is a pandas.DataFrame result of integron_max
        local_max_done = True
        attc.sort_values(["Accession_number", "pos_beg", "evalue"], inplace=True)
    elif attc_file:
        # it call after default search
        local_max_done = False
        attc = read_infernal(attc_file,
                             replicon.id,
                             len(replicon),
                             cfg.model_len,
                             evalue=cfg.evalue_attc,
                             size_max_attc=cfg.max_attc_size,
                             size_min_attc=cfg.min_attc_size)
        attc.sort_values(["Accession_number", "pos_beg", "evalue"], inplace=True)

    # attc_cluster_list = list of Dataframe, each have an array of attC
    attc_cluster_list = search_attc(attc, cfg.keep_palindromes, cfg.distance_threshold,
                                    len(replicon), replicon.topology)
    integrons = []
    if not intI_ac.empty and attc_cluster_list:
        attc_cluster_nb = len(attc_cluster_list)
        # If an array hasn't been clustered with an Integrase
        # or if an integrase lacks an array
        # redundant info, we could check for len(attc_cluster_list)==0
        # -> to remove
        for i, id_int in enumerate(intI_ac.ID_prot.values):  # For each Integrase
            if attc_cluster_nb == 0:  # No more array to attribute to an integrase

                integrons.append(Integron(replicon, cfg))
                integrons[-1].add_integrase(intI_ac.pos_beg.values[i],
                                            intI_ac.pos_end.values[i],
                                            id_int,
                                            int(intI_ac.strand.values[i]),
                                            intI_ac.evalue.values[i],
                                            intI_ac.query_name.values[i])

            else:  # we still have several attC cluster and int :
                # Look for array of attC (cluster) where intI would fall inside it

                # array_2_split is a boolean with True when intI is within an array
                array_2_split = [(attc.pos_beg.values[0] < intI_ac.pos_beg.values[i] and
                                  intI_ac.pos_beg.values[i] < attc.pos_end.values[-1])
                                  for attc in attc_cluster_list]
                # get the index of array to split
                split_index = np.where(array_2_split)[0]
                # for each of the attc cluster to split
                # pop it, split it and add the 2 new arrays back.
                for index in split_index:
                    poped_attc = attc_cluster_list.pop(index)
                    whr_split = np.searchsorted(poped_attc.pos_beg.values, intI_ac.pos_beg.values[i])
                    for split_item in poped_attc.iloc[:whr_split], poped_attc.iloc[whr_split:]:
                        # when there is only one attC in the cluster
                        # the split generate an emtpy dataframe
                        if not split_item.empty:
                            attc_cluster_list.append(split_item)
                    attc_cluster_nb += 1  # new attC array

                attc_left = np.array([i_attc.pos_beg.values[0] for i_attc in attc_cluster_list])
                attc_right = np.array([i_attc.pos_end.values[-1] for i_attc in attc_cluster_list])

                if replicon.topology == 'circ':
                    distances = np.array([(attc_left - intI_ac.pos_end.values[i]),
                                          (intI_ac.pos_beg.values[i] - attc_right)]) % len(replicon)
                else:
                    distances = np.array([abs(attc_left - intI_ac.pos_end.values[i]),
                                          abs(intI_ac.pos_beg.values[i] - attc_right)])

                if attc_cluster_list:
                    # tmp = (distances /
                    #       np.array([[len(aac) for attc in attc_cluster_list]]))
                    side, idx_attc = np.where(distances == distances.min())
                    # side : 0 <=> left; 1 <=> right
                    # index of the closest and biggest attC array to the integrase
                    # exactly tmp = dist(cluster to integrase) / size cluster
                    # to make a decision between 2 equally distant arrays
                    # Usually they are on the same side but on 2 different strands

                    # If they are exactly similar (same distance, same number of attC, take the first one arbitrarily
                    # Or just flatten from idx_attc=[i] to idx_attc=i
                    idx_attc = idx_attc[0]
                    side = side[0]

                else:
                    idx_attc = 0
                    side = np.argmin(distances)

                if distances[side, idx_attc] < cfg.distance_threshold:
                    integrons.append(Integron(replicon, cfg))
                    integrons[-1].add_integrase(intI_ac.pos_beg.values[i],
                                                intI_ac.pos_end.values[i],
                                                id_int,
                                                int(intI_ac.strand.values[i]),
                                                intI_ac.evalue.values[i],
                                                intI_ac.query_name.values[i])

                    attc_tmp = attc_cluster_list.pop(idx_attc)
                    for a_tmp in attc_tmp.values:
                        integrons[-1].add_attC(a_tmp[4],  # pos_beg
                                               a_tmp[5],  # pos_end
                                               1 if a_tmp[6] == "+" else -1,  # sens
                                               a_tmp[7],  # evalue
                                               cfg.model_attc_name
                                               )
                    attc_cluster_nb -= 1

                else:  # no array close to the integrase on both side
                    integrons.append(Integron(replicon, cfg))
                    integrons[-1].add_integrase(intI_ac.pos_beg.values[i],
                                                intI_ac.pos_end.values[i],
                                                id_int,
                                                int(intI_ac.strand.values[i]),
                                                intI_ac.evalue.values[i], intI_ac.query_name.values[i])

        if attc_cluster_nb > 0:  # after the integrase loop (<=> no more integrases)
            for attc_array in attc_cluster_list:
                integrons.append(Integron(replicon, cfg))

                for a_tmp in attc_array.values:
                    integrons[-1].add_attC(a_tmp[4],
                                           a_tmp[5],
                                           1 if a_tmp[6] == "+" else -1,
                                           a_tmp[7], cfg.model_attc_name)

    elif intI_ac.pos_end.values.size == 0 and attc_cluster_list:  # If attC only
        for attc_array in attc_cluster_list:
            integrons.append(Integron(replicon, cfg))
            for a_tmp in attc_array.values:
                integrons[-1].add_attC(a_tmp[4],
                                       a_tmp[5],
                                       1 if a_tmp[6] == "+" else -1,
                                       a_tmp[7], cfg.model_attc_name)

    elif intI_ac.pos_end.values.size >= 1 and not attc_cluster_list:  # If intI only
        for i, id_int in enumerate(intI_ac.ID_prot.values):
            integrons.append(Integron(replicon, cfg))
            integrons[-1].add_integrase(intI_ac.pos_beg.values[i],
                                        intI_ac.pos_end.values[i],
                                        id_int,
                                        int(intI_ac.strand.values[i]),
                                        intI_ac.evalue.values[i],
                                        intI_ac.query_name.values[i])

    #########################################
    # filter CALIN integron on attc number  #
    #########################################
    # Only after local_max if it will be called
    if (cfg.local_max and local_max_done) or not cfg.local_max:
        _log.debug("filter out 'CALIN' with less attC sites than {}".format(cfg.calin_threshold))
        integrons = [i for i in integrons if i.type() != 'CALIN' or len(i.attC) >= cfg.calin_threshold]

    ###############
    # log summary #
    ###############
    _log.info("In replicon {}, there are:".format(replicon.id))
    _log.info("- {} complete integron(s) found with a total {} attC site(s)".format(sum(
        [1 if i.type() == "complete" else 0 for i in integrons]),
        sum([len(i.attC) if i.type() == "complete" else 0 for i in integrons])))
    _log.info("- {} CALIN element(s) found with a total of {} attC site(s)".format(sum(
        [1 if i.type() == "CALIN" else 0 for i in integrons]),
        sum([len(i.attC) if i.type() == "CALIN" else 0 for i in integrons])))
    _log.info("- {} In0 element(s) found with a total of {} attC site".format(sum(
        [1 if i.type() == "In0" else 0 for i in integrons]),
        sum([len(i.attC) if i.type() == "In0" else 0 for i in integrons])))

    return integrons


class Integron(object):
    """Integron object represents an object composed of an integrase, attC sites and gene cassettes.
    Each element is characterized by their coordinates in the replicon, the strand (+ or -),
    the ID of the gene (except attC).
    The object Integron is also characterized by the ID of the replicon."""

    def __init__(self, replicon, cfg):
        """
        :param replicon: The replicon where integrons has been found
        :type replicon: a :class:`Bio.Seq.SeqRecord` object
        :param cfg: the configuration
        :type cfg: a :class:`integron_finder.config.Config` object
        """
        self.cfg = cfg
        self.replicon = replicon
        self.replicon_size = len(self.replicon)
        self._columns = ["pos_beg", "pos_end", "strand", "evalue", "type_elt", "model", "distance_2attC", "annotation"]
        self._dtype = {"pos_beg": "int",
                       "pos_end": "int",
                       "strand": "int",
                       "evalue": "float",
                       "type_elt": "str",
                       "model": "str",
                       "distance_2attC": "float",
                       "annotation": "str"}
        self.integrase = pd.DataFrame(columns=self._columns)
        self.integrase = self.integrase.astype(dtype=self._dtype)

        self.attC = pd.DataFrame(columns=self._columns)
        self.attC = self.attC.astype(dtype=self._dtype)

        self.promoter = pd.DataFrame(columns=self._columns)
        self.promoter = self.promoter.astype(dtype=self._dtype)

        self.attI = pd.DataFrame(columns=self._columns)
        self.attI = self.attI.astype(dtype=self._dtype)

        self.proteins = pd.DataFrame(columns=self._columns)
        self.proteins = self.proteins.astype(dtype=self._dtype)

        self.sizes_cassettes = None

    @property
    def dtype(self):
        return {k: v for k, v in self._dtype.items()}

    def add_integrase(self, pos_beg_int, pos_end_int, id_int, strand_int, evalue, model):
        """Adds integrases to the integron. Should be called once.

        :param int pos_beg_int: the position on the replicon of the beginning integrase site
        :param int pos_end_int: the position on replicon of the end of the integrase site
        :param str id_int: The protein id corresponding to the integrase
        :param int strand_int: the strand where is found the attc 1 for forward, -1 for reverse
        :param float evalue: the evalue associated to this attc site
        :param str model: the name of integrase model (for instance intersection_tyr_intI)
        """

        if not self.integrase.empty:
            raise RuntimeError("add_integrase should be called once.")
        tmp_df = pd.DataFrame(columns=self._columns)
        tmp_df = tmp_df.astype(dtype=self._dtype)
        tmp_df["pos_beg"] = [pos_beg_int]
        tmp_df["pos_end"] = [pos_end_int]
        tmp_df["strand"] = [strand_int]
        tmp_df["evalue"] = [evalue]
        tmp_df["type_elt"] = "protein"
        tmp_df["annotation"] = "intI"
        tmp_df["model"] = [model]
        tmp_df.index = [id_int]
        tmp_df["distance_2attC"] = [np.nan]
        self.integrase = pd.concat([self.integrase, tmp_df], ignore_index=False)


    def add_attC(self, pos_beg_attC, pos_end_attC, strand, evalue, model):
        """Adds attC site to the Integron object.

        :param int pos_beg_attC: the position on the replicon of the beginning attc site
        :param int pos_end_attC: the position on replicon of the end of the attc site
        :param int strand: the strand where is found the attc 1 for forward, -1 for reverse
        :param float evalue: the evalue associated to this attc site
        :param str model: the name of attc model (for instance attc4)
        """
        tmp_df = pd.DataFrame(columns=self._columns)
        tmp_df = tmp_df.astype(dtype=self._dtype)
        tmp_df["pos_beg"] = [pos_beg_attC]
        tmp_df["pos_end"] = [pos_end_attC]
        tmp_df["strand"] = [strand]
        tmp_df["evalue"] = [evalue]
        tmp_df["type_elt"] = "attC"
        tmp_df["annotation"] = "attC"
        tmp_df["model"] = [model]
        self.attC = pd.concat([self.attC, tmp_df], ignore_index=True)
        attC_len = len(self.attC)
        if attC_len < 2:
            self.sizes_cassettes = [np.nan]
        else:
            self.sizes_cassettes.append((self.attC.iloc[attC_len - 1].pos_beg -
                                         self.attC.iloc[attC_len - 2].pos_end) % len(self.replicon))
        self.attC["distance_2attC"] = self.sizes_cassettes

        # self.attC.sort_values(["pos_beg"], inplace = True)G669
        self.attC.index = ["attc_%03i" % int(j + 1) for j in self.attC.index]


    def type(self):
        """
        :returns: The type of the integrons:

                    - 'complete' : Have one integrase and at least one attC
                    - 'CALIN' : Have at least one attC
                    - 'In0' : Just an integrase intI
        :rtype: str
        """
        if not self.attC.empty and not self.integrase.empty:
            return "complete"
        elif self.attC.empty and not self.integrase.empty:
            return "In0"
        elif not self.attC.empty and self.integrase.empty:
            return "CALIN"


    def add_promoter(self):
        """
        Looks for known promoters if they exists within your integrons element.
        It takes 1s for about 13kb.
        """
        dist_prom = 500  # pb distance from edge of the element for which we seek promoter

        ######## Promoter of integrase #########

        if self.has_integrase():
            # PintI1
            p_intI1_seq = Seq.Seq("TTGCTGCTTGGATGCCCGAGGCATAGACTGTACA")
            p_intI1_mot = motifs.create([p_intI1_seq])
            p_intI1_mot.name = "P_intI1"

            # PintI2
            # Not known

            # PintI3
            # Not known

            motifs_Pint = [p_intI1_mot]

            seq_p_int = self.replicon.seq[int(self.integrase.pos_beg.min()) - dist_prom:
                                          int(self.integrase.pos_end.max()) + dist_prom]

            for m in motifs_Pint:
                if self.integrase.strand.values[0] == 1:
                    generator_motifs = seq_p_int[:dist_prom].search(p_intI1_mot.alignment.sequences)
                    for pos, s in generator_motifs:
                        tmp_df = pd.DataFrame(columns=self._columns)
                        tmp_df = tmp_df.astype(dtype=self._dtype)
                        tmp_df["pos_beg"] = [self.integrase.pos_beg.values[0] - dist_prom + pos]
                        tmp_df["pos_end"] = [self.integrase.pos_beg.values[0] - dist_prom + pos + len(s)]
                        tmp_df["strand"] = [self.integrase.strand.values[0]]
                        tmp_df["evalue"] = [np.nan]
                        tmp_df["type_elt"] = "Promoter"
                        tmp_df["annotation"] = "Pint_%s" %(m.name[-1])
                        tmp_df["model"] = "NA"
                        tmp_df.index = [m.name]
                        tmp_df["distance_2attC"] = [np.nan]
                        self.promoter = pd.concat([self.promoter, tmp_df])
                else:
                    # generator_motifs = m.instances.reverse_complement().search(seq_p_int[-dist_prom:])
                    generator_motifs = seq_p_int[-dist_prom:].search(m.reverse_complement().alignment.sequences)
                    for pos, s in generator_motifs:
                        tmp_df = pd.DataFrame(columns=self._columns)
                        tmp_df = tmp_df.astype(dtype=self._dtype)
                        tmp_df["pos_beg"] = [self.integrase.pos_end.max() + pos]
                        tmp_df["pos_end"] = [self.integrase.pos_end.max() + pos + len(s)]
                        tmp_df["strand"] = [self.integrase.strand.values[0]]
                        tmp_df["evalue"] = [np.nan]
                        tmp_df["type_elt"] = "Promoter"
                        tmp_df["annotation"] = "Pint_%s" % (m.name[-1])
                        tmp_df["model"] = "NA"
                        tmp_df.index = [m.name]
                        tmp_df["distance_2attC"] = [np.nan]
                        self.promoter = pd.concat([self.promoter, tmp_df])
            integrase_start = int(self.integrase.pos_beg.values[0])
            integrase_end = int(self.integrase.pos_end.values[-1])
        ######## Promoter of K7 #########

        # Pc-int1
        motifs_Pc = []

        pc = SeqIO.parse(os.path.join(self.cfg.model_dir, "variants_Pc_intI1.fst"), "fasta")
        d = {}
        for seq_rec in pc:
            seq_len = len(seq_rec)
            d[seq_len] = d.get(seq_len, []) + [seq_rec.seq.upper()]
        if not pc.stream.closed:
            # Biopython open file in SeqIO.parse but not close it
            pc.stream.close()
        
        for k, i in d.items():
            motifs_Pc.append(motifs.create(i))
            motifs_Pc[-1].name = "Pc_int1"

        # Pc-int2
        # Not known

        # Pc-int3
        pc_intI3 = motifs.create([Seq.Seq("TAGACATAAGCTTTCTCGGTCTGTAGGCTGTAATG"),
                                  Seq.Seq("TAGACATAAGCTTTCTCGGTCTGTAGGATGTAATG")])
        pc_intI3.name = "Pc_int3"
        motifs_Pc.append(pc_intI3)

        if not self.attC.empty:
            attc_start = int(self.attC.pos_beg.values[0])
            attc_end = int(self.attC.pos_end.values[-1])

        if self.type() == "complete":
            if self.replicon.topology == 'circ':
                if ((attc_start - integrase_end) % self.replicon_size >
                        (integrase_start - attc_end) % self.replicon_size):
                    # if integrase after attcs (on the right)
                    left = attc_end
                    right = integrase_start
                else:
                    left = integrase_end
                    right = attc_start
            else: # replicon is linear
                if attc_end < integrase_start:
                    # integrase on the right of attC cluster.
                    left = attc_end
                    right = integrase_start
                else:
                    left = integrase_end
                    right = attc_start
            strand_array = self.attC.strand.unique()[0]

        elif self.type() == "In0":
            left = integrase_start
            right = integrase_end
            strand_array = "both"

        elif self.type() == "CALIN":
            left = attc_start
            right = attc_end
            strand_array = self.attC.strand.unique()[0]

        if left < right:
            seq_Pc = self.replicon.seq[left - dist_prom:right + dist_prom]
        else:
            seq_Pc1 = self.replicon.seq[left - dist_prom:self.replicon_size]
            seq_Pc2 = self.replicon.seq[:right + dist_prom]
            seq_Pc = seq_Pc1 + seq_Pc2

        for m in motifs_Pc:
            if strand_array == 1:
                mot = [m]
            elif strand_array == "both":
                mot = [m.reverse_complement(), m]
            else:
                mot = [m.reverse_complement()]

            for sa, mo in enumerate(mot):
                # for pos, s in mo.instances.search(seq_Pc):
                for pos, s in seq_Pc.search(mo.alignment.sequences):
                    tmp_df = pd.DataFrame(columns=self._columns)
                    tmp_df = tmp_df.astype(dtype=self._dtype)
                    tmp_df["pos_beg"] = [(left - dist_prom + pos) % self.replicon_size]
                    tmp_df["pos_end"] = [(left - dist_prom + pos + len(s)) % self.replicon_size]
                    tmp_df["strand"] = [strand_array] if strand_array != "both" else [sa * 2 - 1]
                    tmp_df["evalue"] = [np.nan]
                    tmp_df["type_elt"] = "Promoter"
                    tmp_df["annotation"] = "Pc_%s" % (m.name[-1])
                    tmp_df["model"] = "NA"
                    tmp_df.index = [m.name]
                    tmp_df["distance_2attC"] = [np.nan]
                    self.promoter = pd.concat([self.promoter, tmp_df])


    def add_attI(self):
        """
        Looking for Att1 sites and add them to this integron.
        """
        dist_atti = 500

        # attI1
        instances_attI1 = [Seq.Seq('TGATGTTATGGAGCAGCAACGATGTTACGCAGCAGGGCAGTCGCCCTAAAACAAAGTT')]
        attI1 = motifs.create(instances_attI1)
        attI1.name = "attI1"

        # attI2
        instances_attI2 = [Seq.Seq('TTAATTAACGGTAAGCATCAGCGGGTGACAAAACGAGCATGCTTACTAATAAAATGTT')]
        attI2 = motifs.create(instances_attI2)
        attI2.name = "attI2"

        # attI3
        instances_attI3 = [Seq.Seq('CTTTGTTTAACGACCACGGTTGTGGGTATCCGGTGTTTGGTCAGATAAACCACAAGTT')]
        attI3 = motifs.create(instances_attI3)
        attI3.name = "attI3"

        motif_attI = [attI1, attI2, attI3]

        if self.type() in ("CALIN", "complete"):
            attc_start = self.attC.pos_beg.values[0]
            attc_end = self.attC.pos_end.values[-1]

        if self.type() in ("complete", "In0"):
            integrase_start = self.integrase.pos_beg.values[0]
            integrase_end = self.integrase.pos_end.values[-1]

        if self.type() == "complete":
            if self.replicon.topology == 'circ':
                if ((attc_start - integrase_end) % self.replicon_size >
                        (integrase_start - attc_end) % self.replicon_size):
                    # if integrase after attcs (on the right)
                    left = attc_end
                    right = integrase_start
                else:
                    left = integrase_end
                    right = attc_start
            else: # replicon is linear
                if attc_end < integrase_start:
                    # integrase on the right of attC cluster.
                    left = attc_end
                    right = integrase_start
                else:
                    left = integrase_end
                    right = attc_start

            strand_array = self.attC.strand.unique()[0]

        elif self.type() == "In0":
            left = int(self.integrase.pos_beg.iloc[0])
            right = int(self.integrase.pos_end.iloc[0])
            strand_array = "both"
        elif self.type() == "CALIN":
            left = attc_start
            right = attc_end
            strand_array = self.attC.strand.unique()[0]

        if left < right:
            seq_attI = self.replicon.seq[left - dist_atti:right + dist_atti]
        else:
            seq_attI1 = self.replicon.seq[left - dist_atti:self.replicon_size]
            seq_attI2 = self.replicon.seq[:right + dist_atti]
            seq_attI = seq_attI1 + seq_attI2

        for m in motif_attI:

            if strand_array == 1:
                mot = [m]
            elif strand_array == "both":
                mot = [m.reverse_complement(), m]
            else:
                mot = [m.reverse_complement()]

            for sa, mo in enumerate(mot):
                for pos, s in seq_attI.search(mo.alignment.sequences):  #  mo.instances.search(seq_attI) is depprecated
                    tmp_df = pd.DataFrame(columns=self._columns)
                    tmp_df = tmp_df.astype(dtype=self._dtype)
                    tmp_df["pos_beg"] = [(left - dist_atti + pos) % self.replicon_size]
                    tmp_df["pos_end"] = [(left - dist_atti + pos + len(s)) % self.replicon_size]
                    tmp_df["strand"] = [strand_array] if strand_array != "both" else [sa * 2 - 1]
                    tmp_df["evalue"] = [np.nan]
                    tmp_df["type_elt"] = "attI"
                    tmp_df["annotation"] = f"attI_{m.name[-1]}"
                    tmp_df["model"] = "NA"
                    tmp_df.index = [m.name]
                    tmp_df["distance_2attC"] = [np.nan]
                    self.attI = pd.concat([self.attI, tmp_df])


    def add_proteins(self, prot_db):
        """
        :param prot_db: The protein db corresponding to the translation of the replicon
        :type prot_db: :class:`integron.prot_db.ProteinDB` object.
        """


        def to_add(window_start, window_end, prot_attr):
            """
            decide if we keep the protein or not

            We keep proteins (<--->) if start (<) and end (>) follows that scheme:

                 ok:            <--->         <--->
                 ok:  <--->                                    <--->
                          ^ 200pb v                    v 200pb ^
                                  |------integron------|
                                window_start                 fin
            """
            if self.replicon.topology == 'circ':
                s_int = (window_end - window_start) % self.replicon_size
                return ((window_end - prot_attr.stop) % self.replicon_size < s_int) or \
                       ((prot_attr.start - window_start) % self.replicon_size < s_int)
            else:
                return window_start < prot_attr.stop < window_end or window_start < prot_attr.start < window_end

        attc_start = self.attC.pos_beg.values[0]
        attc_end = self.attC.pos_end.values[-1]

        if self.has_integrase():
        
            integrase_start = self.integrase.pos_beg.values[0]
            integrase_end = self.integrase.pos_end.values[-1]

            if self.replicon.topology == 'circ':
                if ((attc_start - integrase_end) % self.replicon_size >
                        (integrase_start - attc_end) % self.replicon_size):
                    # integrase on the right of attC cluster.
                    window_start = attc_start - 200
                    window_end = self.integrase.pos_beg.min()
                else:
                    window_start = self.integrase.pos_end.max()
                    window_end = attc_end + 200
            else:  # replicon is linear
                if attc_end < integrase_start:
                    # integrase on the right of attC cluster.
                    window_start = max(attc_start - 200, 0)
                    window_end = self.integrase.pos_beg.min()

                else:
                    # integrase on the left of attC cluster.
                    window_start = self.integrase.pos_end.max()
                    window_end = min(attc_end + 200, self.replicon_size)

        else:
            # To allow the first protein after last attC to aggregate.
            window_start = attc_start - 200
            window_end = attc_end + 200

        for prot_id in prot_db:
            prot_attr = prot_db.get_description(prot_id)
            if to_add(window_start, window_end, prot_attr):
                prot_annot = "protein"
                prot_evalue = np.nan
                prot_model = "NA"
                self.proteins.loc[prot_attr.id] = [prot_attr.start, prot_attr.stop, prot_attr.strand, prot_evalue,
                                                   "protein", prot_model, np.nan, prot_annot]

            intcols = ["pos_beg", "pos_end", "strand"]
            floatcols = ["evalue", "distance_2attC"]
            self.proteins[intcols] = self.proteins[intcols].astype(int)
            self.proteins[floatcols] = self.proteins[floatcols].astype(float)

    def describe(self):
        """
        :returns: DataFrame describing the integron object
                  The columns are:

                  "pos_beg", "pos_end", "strand", "evalue", "type_elt", "model",
                  "distance_2attC", "annotation", "considered_topology"

        """
        full = pd.concat([self.integrase, self.attC, self.promoter, self.attI, self.proteins])
        full["pos_beg"] = full["pos_beg"].astype(int)
        full["pos_end"] = full["pos_end"].astype(int)
        full["strand"] = full["strand"].astype(int)
        full["distance_2attC"] = full["distance_2attC"].astype(float)
        full = full.reset_index()
        full.columns = ["element"] + list(full.columns[1:])
        full["type"] = self.type()
        full["ID_replicon"] = self.replicon.id
        full["ID_integron"] = id(self)  # uniq identifier of a given Integron
        full["default"] = "Yes" if not self.cfg.local_max else "No"
        try:
            # when replicon has been got using utils.FastaIterator
            full["considered_topology"] = self.replicon.topology
        except AttributeError:
            # if replicon is a bare Bio.SeqRecord
            full["considered_topology"] = self.cfg.default_topology

        full.drop_duplicates(subset=["element"], inplace=True)
        return full


    def draw_integron(self, file=None):
        """
        Represent the different element of the integrons if file is provide
        save the drawing on the file otherwise display it on screen.

        :param str file: the path to save the integron schema (in pdf format)
        """
        full = self.describe()
        full["evalue"] = full["evalue"].astype("float")
        h = [i + (0.5*i) if j == "Promoter" else i for i, j in zip(full.strand, full.type_elt)]
        fig, ax = plt.subplots(1, 1, figsize=(16, 9))
        alpha = [i if i < 1 else 1 for i in (
                (np.log10(full.evalue) - np.ones(len(full)) * -1) /
                (np.ones(len(full)) * -10 - np.ones(len(full)) * -1)
                * (1 - 0.2) + 0.2).fillna(1).tolist()]
        # normalize alpha value with 0.2 as min value

        colors = ["#749FCD" if i == "attC" else
                  "#DD654B" if i == "intI" else
                  "#6BC865" if (i[-2:] == "_1" and j == "Promoter") else
                  "#D06CC0" if (i[-2:] == "_2" and j == "Promoter") else
                  "#C3B639" if (i[-2:] == "_3" and j == "Promoter") else
                  "#e8950e" if i != "protein" else
                  "#d3d3d3" for (i, j) in zip(full.annotation,
                                              full.type_elt)]

        # colors_alpha = [j+[i] for j, i in zip([[ord(c) / 255. for c in i[1:].decode("hex")] for i in colors],
        #                                      alpha)]

        colors_alpha = [matplotlib.colors.to_rgba_array(c, a)[0].tolist() for c, a in zip(colors, alpha)]

        # ec = ["red" if i =="attC" else
        #      "white" for i in full.type_elt]
        # z_order = [100 if i == "attC" else 1 for i in full.type_elt]
        z_order = 10
        ax.barh(np.zeros(len(full)), full.pos_end-full.pos_beg,
                height=h, left=full.pos_beg,
                color=colors_alpha, zorder=z_order, ec=None,
                align="edge")  # edgecolor=ec,
        xlims = ax.get_xlim()
        for color, label in zip(["#749FCD", "#DD654B", "#6BC865", "#D06CC0", "#C3B639", "#e8950e", "#d3d3d3"],
                        ["attC", "integrase", "Promoter/attI class 1",
                         "Promoter/attI class 2", "Promoter/attI class 3",
                         "Functional Annotation", "Hypothetical Protein"]):
            ax.bar(0, 0, color=color, label=label)
        plt.legend(loc=[1.01, 0.4])
        ax.set_xlim(xlims)
        fig.subplots_adjust(left=0.05, right=0.80)
        ax.hlines(0, ax.get_xlim()[0], ax.get_xlim()[1], "lightgrey", "--")
        ax.grid(True, "major", axis="x")
        ax.set_ylim(-4, 4)
        ax.get_yaxis().set_visible(False)
        if file:
            fig.savefig(file, format="pdf")
            plt.close(fig)
        else:
            fig.show()


    def has_integrase(self):
        """
        :return: True if integron has integrase False otherwise.
        """
        return not self.integrase.empty


    def has_attC(self):
        """
        :return: True if integron has attc sites False otherwise.
        """
        return not self.attC.empty
