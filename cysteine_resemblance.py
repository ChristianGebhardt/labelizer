# -*- coding: utf-8 -*-
# Copyright (C) 2020, Christian Gebhardt (chrirstian.gebhardt@bio.lmu.de)
#
# This file is part of the labelizer distribution and governed by your
# choice of the "MIT license".
# Please see the LICENSE file that should have been included as part of this
# package.

"""
Calculation of cysteine ressemblance labelling parameter.
"""

from .labeling_parameter import LabelingParameter


class CysteineResemblance(LabelingParameter):
    """
    Class implements the cystein resemblance labelling parameter as subclass of
    abstract class LabellingParameter. It overrides function calc_parameter_score.

    """

    # constants / parameter
    file_tag = "cr"
    # RESEMBLANCES = {'ALA': 1., 'ARG': 0.5, 'ASN': 0.5, 'ASP': 0.5, 'CYS': 1.,
    #                 'GLU': 0.5, 'GLN': 0.5, 'GLY': 0.75, 'HIS': 0.5, 'ILE': 0.75,
    #                 'LEU': 0.75, 'LYS': 0.5, 'MET': 0.75, 'PHE': 0.75, 'PRO': 0.75,
    #                 'SER': 1., 'THR': 0.5, 'TRP': 0.75, 'TYR': 0.5, 'VAL': 0.75}
    RES_SHORT_DICT = {
        "ALA": "A",
        "ARG": "R",
        "ASN": "N",
        "ASP": "D",
        "CYS": "C",
        "GLU": "E",
        "GLN": "Q",
        "GLY": "G",
        "HIS": "H",
        "ILE": "I",
        "LEU": "L",
        "LYS": "K",
        "MET": "M",
        "PHE": "F",
        "PRO": "P",
        "SER": "S",
        "THR": "T",
        "TRP": "W",
        "TYR": "Y",
        "VAL": "V",
    }

    # RESEMBLANCES_SHORT = {'A': 0., 'R': 0., 'N': 0., 'D': 0., 'C': 1.,
    #                       'E': 0., 'Q': 0., 'G': 0., 'H': 0., 'I': 0.,
    #                       'L': 0., 'K': 0., 'M': 0., 'F': 0., 'P': 0.,
    #                       'S': 0., 'T': 0., 'W': 0., 'Y': 0., 'V': 0.}
    # class variables

    def calc_parameter_score(self, chain_id, res_id_int):
        """
        Funcion to calculate parameter score for cystein resemblance for one residue.
        Overwritten function from LabelingParameter.

        :param chain_id: Chain ID in model.
        :type chain_id: string
        :param res_id_int: Residue number within chain.
        :type res_id_int: int
        :return: Cystein resemblance score for residue
                 (value betwenn 0 and 1, or -1. on error).
        :rtype: float
        """
        res_name = self.model[chain_id][res_id_int].get_resname()
        cr_score = self._resemblance_score(res_name)
        return cr_score

    def _resemblance_score(self, residue_name):
        if self.parameter_model == "C_CR1_Name":
            value = self.RES_SHORT_DICT[residue_name]
        else:
            raise NotImplementedError("Only 'C_CR1_Name' is implemented so far")
            # if residue_name in self.RESEMBLANCES:
            #     cr_score = self.RESEMBLANCES[residue_name]
            # else:
            #     cr_score = -1.
        cr_score = self.calc_score_frequency(value)
        return cr_score
