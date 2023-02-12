# -*- coding: utf-8 -*-
# Copyright (C) 2020, Christian Gebhardt (chrirstian.gebhardt@bio.lmu.de)
#
# This file is part of the labelizer distribution and governed by your
# choice of the "MIT license".
# Please see the LICENSE file that should have been included as part of this
# package.

"""
Calculation of solvent exposure labelling parameter.
"""
import logging
import os
import platform
import warnings

from Bio import PDB
from Bio.PDB import DSSP, PDBParser
from Bio.PDB import HSExposure
from Bio.PDB.ResidueDepth import get_surface, residue_depth

from .auxiliary_functions import find_msms
from .labeling_parameter import LabelingParameter


class SolventExposure(LabelingParameter):
    """
    Class implements the solventexposure labelling parameter as subclass of
    abstract class LabellingParameter. It overrides function calc_parameter_score.

    """

    # constants / parameter
    file_tag = "se"
    HSE_RADIUS = 13.0
    SUM_HSE_THRESH = 30
    SUM_HSE_MAX = 45
    SLOPE_OFFSET = 5
    HSE1_THRESH = 10
    HSE1_MAX = 20
    half_sphere_exp = None
    dssp = None
    msms: str = "msms"

    def __init__(self) -> None:
        super().__init__()
        self.msms = find_msms(self.msms)
        self.surface = None

    def set_up(self, labelizer, protein, sensitivity="m", para_model=None):
        """
        Load and prepare stucture to be analyzed and set analysis settings
        and calculate half-sqhere-exposure.
        Function overwrites function from LabelingParameter.

        :param labelizer: Labelling analysis class with analysis settings.
        :type labelizer: Labelizer
        :param protein: name of protein pdb-file to be analysized.
        :type protein: string
        :param sensitivity: Sensitivity of parameter score.
                            Must be high 'h', medium 'm', or low 'l',
                            defaults to 'm'.
        :type sensitivity: string, optional
        """
        super(SolventExposure, self).set_up(labelizer, protein, sensitivity, para_model)
        if "HSE" in para_model:
            self._calc_hse()
        filepath = os.path.join(
            os.path.dirname(self.file_path), self.identifier + ".pdb"
        )

        if para_model == "N_SE1_RSA_WILKE":
            warnings.simplefilter("ignore")  # suppress dssp user warnings
            # warning
            # Residue RES N A has X instead of expected Y sidechain atoms.
            # Calculation solvent accesibility refers to incomplete sidechain.
            if platform.system() == "Windows":
                logging.info("-- windows --")
                warnings.simplefilter("ignore")
                parser = PDB.PDBParser()
                # io = PDB.PDBIO()
                structure = parser.get_structure("X", filepath)
                model = structure[0]
                self.dssp = DSSP(model, filepath, dssp="dsspcmbi")
            elif platform.system() == "Linux":
                logging.info("-- linux --")
                warnings.simplefilter("ignore")
                parser = PDB.PDBParser()
                # io = PDB.PDBIO()
                structure = parser.get_structure("X", filepath)
                model = structure[0]
                self.dssp = DSSP(model, filepath)
            else:
                raise Exception("Unknown platform")
            # self.dssp = DSSP(self.model, self.file_path, dssp="dssp\dsspcmbi", acc_array='Wilke')
        if "SURFACE_DIST" in para_model:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")  # suppress dssp user warnings
                # warning
                # Residue RES N A has X instead of expected Y sidechain atoms.
                # Calculation solvent accessibility refers to incomplete sidechain.
                parser = PDBParser()
                structure = parser.get_structure("X", filepath)
                model = structure[0]
                self.surface = get_surface(model, MSMS=self.msms)
                # TODO check if surface is not complete (several not connected protein parts)
                # print(surface)
                # with open(r"C:\Users\gebha\Downloads\surface_all.xyz","w") as f:
                #    for p in surface:
                #        f.write("S "+str(p[0])+" "+str(p[1])+" "+str(p[2])+'\n')

    def calc_parameter_score(self, chain_id: str, res_id_int: int) -> float:
        """
        Function to calculate parameter score for solvent exposure for one residue.
        Overwritten function from LabelingParameter.

        :param chain_id: Chain ID in model.
        :param res_id_int: Residue number within chain.
        :return: Solvent exposure score for residue (value between 0 and 1).
        """
        if "HSE" in self.parameter_model:
            se_value = self.half_sphere_exp[chain_id + str(res_id_int)]
        elif self.parameter_model == "N_SE1_RSA_WILKE":
            se_value = self.dssp[(chain_id, res_id_int)][3]
        elif self.parameter_model == "N_SE11_MEAN_SURFACE_DIST":
            res = self.model[chain_id][res_id_int]
            se_value = residue_depth(res, self.surface)
        else:
            raise NotImplementedError(
                "Only 'N_CS2_Score' and 'N_SE1_RSA_WILKE' are implemented so far"
            )

        se_score = self._solvent_exposure(se_value)
        return se_score

    def _calc_hse(self):
        hse_radius = int(self.parameter_model[-3:-1])
        # hse_cb = HSExposure.HSExposureCB(self.model, self.HSE_RADIUS)
        hse_cb = HSExposure.HSExposureCB(self.model, hse_radius)
        self.half_sphere_exp = {}
        for key in hse_cb.keys():
            chain_id = key[0]
            amino_nbr = key[1][1]
            self.half_sphere_exp[chain_id + str(amino_nbr)] = hse_cb[key][0:2]

    def _solvent_exposure(self, se_value):
        if self.parameter_model == "I_SE4_HSE1_10A":
            value = se_value[0]
        elif self.parameter_model == "N_SE1_RSA_WILKE":
            value = se_value
        elif self.parameter_model == "N_SE11_MEAN_SURFACE_DIST":
            value = min(se_value, 4.0)
        else:
            raise NotImplementedError(
                "Only 'N_CS2_Score' and 'N_SE1_RSA_WILKE' are implemented so far"
            )
        con_score = self.calc_score_frequency(value)
        return con_score
