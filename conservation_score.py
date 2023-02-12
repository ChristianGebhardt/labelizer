# -*- coding: utf-8 -*-
# Copyright (C) 2020, Christian Gebhardt (chrirstian.gebhardt@bio.lmu.de)
#
# This file is part of the labelizer distribution and governed by your
# choice of the "MIT license".
# Please see the LICENSE file that should have been included as part of this
# package.

"""
Calculation of conservation score labelling parameter.
"""

# import logging
import os
import warnings

from Bio import PDB

from .labeling_parameter import LabelingParameter


class ConservationScore(LabelingParameter):
    """
    Class implements the conservation score labelling parameter as subclass of
    abstract class LabellingParameter. It overrides function calc_parameter_score.

    """

    # constants / parameter
    file_tag = "cs"
    # different scores for sensitivity ['l', 'm', 'h']
    MAX_SCORE = [0.78, 0.84, 1.0]  # threshold for conservation score 1
    # (lower limit when parameter is 1)
    MIN_SCORE = [-0.04, 0.12, 0.28]  # threshold for conservation score 4
    # (theoretical limit when parameter is 0)
    THRESHOLD = [0.2, 0.36, 0.52]  # threshold for conservation score 3
    # (0 below threshold)

    # class variables
    cs_files = {}

    def __init__(self):
        self.cs_files = {}

    def set_up(self, labelizer, protein, sensitivity="m", para_model=None):
        """
        Load and prepare stucture to be analyzed and set analyse settings.
        Load conservation scores.
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
        super(ConservationScore, self).set_up(
            labelizer, protein, sensitivity, para_model
        )
        protein_position = labelizer.proteins.index(protein)
        self.cs_files = {}  # needed because __init__ is only called once
        for chain, prot in zip(labelizer.chains1, labelizer.prot_cs[protein_position]):
            if chain not in self.cs_files:
                # file_path = os.path.join(labelizer.folder, labelizer.file_extension(p,'pdb'))
                file_path = os.path.join(labelizer.folder, prot)
                self._load_cs_file(labelizer, file_path, protein, chain)

    def calc_parameter_score(self, chain_id, res_id_int):
        """
        Funcion to calculate parameter score for conservation score for one residue.
        Overwritten function from LabelingParameter.

        :param chain_id: Chain ID in model.
        :type chain_id: string
        :param res_id_int: Residue number within chain.
        :type res_id_int: int
        :return: Conservation score (value betwenn 0 and 1, on error -1).
        :rtype: float
        """
        b_factor = self.cs_files[chain_id][chain_id + str(res_id_int)]
        cs_score = self._conservation_score(b_factor)
        return cs_score

    def _load_cs_file(self, labelizer, file_path, identifier, tag):
        cs_dict = {}
        try:
            ext = file_path.split(".")[-1]
        except Exception:
            raise ValueError("unknown file extension: {}".format("None"))
        if ext == "pdb":
            cs_dict = self._load_cs_pdb_file(file_path, identifier, tag)
        elif ext == "grades":
            cs_dict = self._load_cs_grades_file(file_path)
        else:
            raise ValueError("unknown file extension: {}".format(ext))
        self.cs_files[tag] = cs_dict

    @staticmethod
    def _load_cs_pdb_file(file_path, identifier, tag):
        cs_dict = {}
        parser = PDB.PDBParser(PERMISSIVE=1)
        structure = None
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            structure = parser.get_structure(identifier, file_path)
        for chain in structure[0]:
            chain_id = chain.get_id()
            if chain_id == tag:
                for residue in chain:
                    try:
                        res_id_int = str(residue.get_id()[1])
                        b_factor = next(residue.get_atoms()).get_bfactor()
                        cs_dict[chain_id + res_id_int] = b_factor
                    except Exception:
                        continue
        return cs_dict

    @staticmethod
    def _load_cs_grades_file(file_path):
        cs_dict = {}
        with open(file_path, "r") as grades:
            lines = grades.readlines()
            for line in lines[15:]:
                try:
                    line_split = line.split("\t")
                    pos = line_split[0].strip()
                    if pos.isdigit():
                        pos = int(pos)
                        res_id_int = line_split[2].strip().split(":")[0][3:]
                        chain_id = line_split[2].strip().split(":")[1]
                        score = float(line_split[3].strip())
                        cs_dict[chain_id + res_id_int] = score
                    else:
                        continue
                except Exception:
                    continue
        return cs_dict

    def _conservation_score(self, conservation):
        if self.parameter_model == "N_CS2_Score":
            value = conservation
        else:
            raise NotImplementedError("Only 'N_CS2_Score' is implemented so far")
        con_score = self.calc_score_frequency(value)
        return con_score
