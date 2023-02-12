# -*- coding: utf-8 -*-
# Copyright (C) 2020, Christian Gebhardt (chrirstian.gebhardt@bio.lmu.de)
#
# This file is part of the labelizer distribution and governed by your
# choice of the "MIT license".
# Please see the LICENSE file that should have been included as part of this
# package.

"""
Calculation of secondary structure labelling parameter.
"""
import logging
import os
import platform
import warnings

# if not sys.warnoptions:
#     print("SUPPRESS")
#     import warnings
#     warnings.simplefilter("ignore")
from Bio import PDB
from Bio.PDB import DSSP

from .labeling_parameter import LabelingParameter

logger = logging.getLogger(__name__)


class SecondaryStructure(LabelingParameter):
    """
    Class implements the secondary structure labelling parameter as subclass of
    abstract class LabellingParameter. It overrides function calc_parameter_score.

    """

    # constants / parameter
    file_tag = "ss"
    # STRUCTURE_VALUE = {'H': alpha-helix, 'B': beta-bridge,
    #                    'E': strand (beta-sheet), 'G': 3-10-helix,
    #                    'I': pi-helix, 'T': turn, 'S': bend, '-': none}
    STRUCTURE_VALUE = {
        "H": 0.3,
        "B": 0.1,
        "E": 0.2,
        "G": 0.05,
        "I": 0.15,
        "T": 1.0,
        "S": 1.0,
        "-": 1.0,
    }
    # relative accesible surface area
    # https://en.wikipedia.org/wiki/Relative_accessible_surface_area
    MAX_ASA = {
        "ALA": 106.0,
        "ARG": 248.0,
        "ASN": 157.0,
        "ASP": 163.0,
        "CYS": 135.0,
        "GLN": 198.0,
        "GLU": 194.0,
        "GLY": 84.0,
        "HIS": 184.0,
        "ILE": 169.0,
        "LEU": 164.0,
        "LYS": 205.0,
        "MET": 188.0,
        "PHE": 197.0,
        "PRO": 136.0,
        "SER": 130.0,
        "THR": 142.0,
        "TRP": 227.0,
        "TYR": 222.0,
        "VAL": 142.0,
    }
    # class variables
    dssp = None
    structurePosition = None

    def set_up(self, labelizer, protein, sensitivity="m", para_model=None):
        """
        Load and prepare stucture to be analyzed and set analysis settings
        and secondary structure with dssp.
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
        super(SecondaryStructure, self).set_up(
            labelizer, protein, sensitivity, para_model
        )
        warnings.simplefilter("ignore")  # suppress dssp user warnings
        # warning
        # Residue RES N A has X instead of expected Y sidechain atoms.
        # Calculation solvent accessibility refers to incomplete sidechain.
        filepath = os.path.join(
            os.path.dirname(self.file_path), self.identifier + ".pdb"
        )
        logger.debug(f"Using pdb {filepath} for SS")
        if platform.system() == "Windows":
            logger.info("-- windows --")
            # script_path = os.path.dirname(os.path.realpath(__file__))
            # dssp_path = os.path.join(script_path,"dssp","dsspcmbi")
            # self.dssp = DSSP(self.model, filepath, dssp=dssp_path)
            # sys.path.insert(0,os.path.join(script_path,"dssp"))
            warnings.simplefilter("ignore")
            parser = PDB.PDBParser()
            structure = parser.get_structure("X", filepath)
            model = structure[0]
            self.dssp = DSSP(model, filepath, dssp="dsspcmbi")
        elif platform.system() == "Linux":
            warnings.simplefilter("ignore")
            parser = PDB.PDBParser()
            structure = parser.get_structure("X", filepath)
            model = structure[0]
            self.dssp = DSSP(model, filepath)
            logger.info(f"KEYS {self.dssp.keys()}")
        else:
            raise RuntimeError("Unknown platform")

    def calc_parameter_score(self, chain_id: str, res_id_int: int):
        """
        Funcion to calculate parameter score for secondary structure for one residue.
        Overwritten function from LabelingParameter.

        :param chain_id: Chain ID in model.
        :type chain_id: string
        :param res_id_int: Residue number within chain.
        :type res_id_int: int
        :return: Secondary structure score for residue (value betwenn 0 and 1).
        :rtype: float
        """
        # lastPos = 0
        entry = self.dssp[(chain_id, res_id_int)]
        # while(self.listDSSP[lastPos][0]<=res_id_int):
        #     if(self.listDSSP[lastPos][0]==res_id_int):
        #         strTag = self.listDSSP[lastPos][2]
        #     lastPos += 1
        str_tag = entry[2]
        ss_score = self._structure_score_simple(str_tag)

        return ss_score

    def _structure_score_simple(self, structure_code):
        # print("PARA CS ",self.parameter_model)
        if self.parameter_model == "C_SS1_SS":
            value = structure_code
        else:
            raise NotImplementedError("Only 'C_SS1_SS' is implemented so far")
        cr_score = self.calc_score_frequency(value)
        return cr_score
