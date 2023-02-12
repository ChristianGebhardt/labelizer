# -*- coding: utf-8 -*-
# Copyright (C) 2020, Christian Gebhardt (chrirstian.gebhardt@bio.lmu.de)
#
# This file is part of the labelizer distribution and governed by your
# choice of the "MIT license".
# Please see the LICENSE file that should have been included as part of this
# package.

"""
The module provides an abstract class for all labelling parameters, which inherit from it.
"""

import csv
import json
import logging
import os
from abc import ABC
from . import pdbhelper

logger = logging.getLogger(__name__)


class LabelingParameter(ABC):
    """
    Abstract class, which provides generic functionality to analyse each
    residue according to a parameter and assign parameter score to it.
    """

    file_tag: str
    sensitivity = "m"
    sens_idx = 1
    parameter_score = {}
    identifier = None
    file_path = None

    structure = None
    model = None
    chains = None

    SENSITIVITY_RANGE = ["l", "m", "h"]

    parameter_model = None

    save_pdb = False

    # class variable
    skript_folder = os.path.dirname(os.path.realpath(__file__))
    logic_folder = os.path.join(skript_folder, "resources", "probabilities")

    @staticmethod
    def factory(tag: str):
        """
        Construct specific parameter class as subclass from LabelingParameter
        specified by its parameter tag.

        :param tag: Two letter tag for labelling parameter
                  (e.g 'se' for solvent exposure).
        :type tag: string
        :return: New instance of LabelingParameter subclass with tag t.
                 If tag t unknown, instance of LabelingParameter.
        :rtype: LabellingParameter, subclass
        """
        # Avoid cyclic imports
        from . import available_parameter

        mapping = {klass.file_tag: klass for klass in available_parameter}
        return mapping[tag]()

    # @abstractmethod
    def set_up(self, labelizer, protein, sensitivity="m", para_model=None):
        """
        Load and prepare structure to be analyzed and set analyze settings.
        Load parameter logic.
        If further preparation is needed for some parameter, this class
        can be overwritten.

        :param labelizer: Labelling analysis class with analysis settings.
        :type labelizer: Labelizer
        :param protein: name of protein pdb-file to be analysized.
        :type protein: string
        :param sensitivity: Sensitivity of parameter score.
                            Must be high 'h', medium 'm', or low 'l',
                            defaults to 'm'.
        :type sensitivity: string, optional
        """
        path1 = os.path.join(labelizer.folder, labelizer.file_ext(protein, "pdb"))
        self._load_pdb(path1, protein)
        pos = labelizer.proteins.index(protein)
        if pos == 0:
            self.chains = labelizer.chains1
        elif pos == 1:
            self.chains = labelizer.chains2
        else:
            NotImplementedError("Maximum 2 protein structures are supported!")
        self.save_pdb = labelizer.SAVE_PDB
        assert sensitivity in self.SENSITIVITY_RANGE, "invalid sensitivity <{}>".format(
            sensitivity
        )
        self.sensitivity = sensitivity
        self.sens_idx = self.SENSITIVITY_RANGE.index(self.sensitivity)

        self.parameter_model = para_model
        self._load_logic(para_model)

    # @abstractmethod
    def clean_up(self, save_pdb=True):
        """
        Save parameter scores and structure with scores in pdb-function (if specified).
        If further clean-up is needed for some parameter, this class can be overwritten.

        """
        if save_pdb:
            self._save_pdb()
        self._save_csv()

    def calc_parameter_scores(self):
        """
        Iterate over all chains for analysis and over all residues therein.
        Function calls abstract function calc_parameter_score for each residue,
        which needs to be overwritten in subclass.

        """
        for chain in self.model:
            chain_id = chain.get_id()
            if chain_id not in self.chains:
                logger.debug("Skip chain %s", chain_id)
                for residue in chain:
                    for atom in residue:
                        atom.set_bfactor(-1)
            else:
                logger.info("%s %s %s", self.__class__, self.identifier, chain_id)
                for residue in chain:
                    res_id = residue.get_id()
                    parameter_key = chain_id + str(res_id[1])
                    try:
                        b_factor = self.calc_parameter_score(chain_id, res_id[1])
                        for atom in residue:
                            atom.set_bfactor(b_factor)
                        self.parameter_score[parameter_key] = b_factor
                    except Exception as e:
                        logger.exception(
                            f"Parameter calculation failed {self.__class__.__name__} {chain_id} {res_id}: {e}"
                        )
                        # TBD DEBUG ONLY WITH DEBUG FLAG
                        # logger.debug(e)
                        self.parameter_score[parameter_key] = -1
                        for atom in residue:
                            atom.set_bfactor(-1)

    # @abstractmethod
    def calc_parameter_score(self, chain_id, res_id_int):
        """
        Abstract function to calculate specific labelling parameter for one residue.
        This needs to be overwritten in subclass.

        :param chain_id: Chain ID in model.
        :type chain_id: string
        :param res_id_int: Residue number within chain.
        :type res_id_int: int
        :return: Labelling parameter score for residue
                 (value betwenn 0 and 1, or -1. on error).
        :rtype: float
        """
        return -1.0

    def calc_score_frequency(self, value):
        if self.para_type == "C":
            return self.para_dict[value]
        else:
            return self.para_dict.get(
                value,
                self.para_dict[
                    min(
                        self.para_dict.keys(),
                        key=lambda k: abs(float(k) - float(value)),
                    )
                ],
            )

    def _load_logic(self, ls_para):
        """
        Load parameter logic from file.

        :param ls_para: Paremeter identifier as in publication.
        :type ls_para: string
        """
        dict_file = os.path.join(
            LabelingParameter.logic_folder, ls_para + "_P_l_after_s.json"
        )
        with open(dict_file, "r") as f:
            p_dict = json.loads(f.read())
            # p_dict = json.loads(open(dict_file))
            if ls_para[:2] == "C_":
                self.para_dict = {k: float(v) for k, v in p_dict.items()}
            else:
                self.para_dict = {float(k): float(v) for k, v in p_dict.items()}
            self.para_type = ls_para[:1]

    def _load_pdb(self, file_path, identifier):
        """
        Load pdb-structure from file and initiate variables in class.

        :param file_path: File path to pdb structure.
        :type file_path: string
        :param identifier: Four letter pdb identifier.
        :type identifier: string
        """
        self.file_path = file_path
        self.identifier = identifier

        self.structure, self.model = pdbhelper.load_pdb(
            identifier, file_path, remove_hetatm=True
        )

        self.parameter_score = {}

    def _save_pdb(self):
        """
        Save pdb-structure with parameter scores in b-factor.
        """
        if self.save_pdb:
            filepath = os.path.join(
                os.path.dirname(self.file_path),
                self.identifier + "_" + self.file_tag + ".pdb",
            )
            pdbhelper.save_pdb(self.structure, filepath)

    def _save_csv(self):
        """
        Save labelling parameter scores into csv-file.
        """
        header = ["ID", "Parameter Score (" + self.file_tag.upper() + ")"]
        filepath = os.path.join(
            os.path.dirname(self.file_path),
            self.identifier + "_" + self.file_tag + ".csv",
        )
        with open(filepath, "w", newline="") as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(header)
            for key, value in self.parameter_score.items():
                writer.writerow([key, value])

    @classmethod
    def get_parameter_probabilities(cls):
        names = []
        for entry in os.scandir(cls.logic_folder):
            names.append(entry.name.split(".")[0])
        names = list(set(names))
        names.sort()
        return names
