# -*- coding: utf-8 -*-
# Copyright (C) 2020, Christian Gebhardt (chrirstian.gebhardt@bio.lmu.de)
#
# This file is part of the labelizer distribution and governed by your
# choice of the "MIT license".
# Please see the LICENSE file that should have been included as part of this
# package.

"""
Calculation of a combined labeling score based on all labelling parameters.
"""

import csv
import logging
import os
from dataclasses import dataclass
from typing import Dict

import numpy as np
from Bio import PDB
from Bio.PDB import PDBIO

logger = logging.getLogger(__name__)


@dataclass
class LabelParameterConfig:
    """Every parameter (conservation, solvent exposure, etc.) has one of those"""

    # 0/1/2
    weight: int
    # high 'h', medium 'm', or low 'l'
    sensitivity: str


@dataclass
class AAExclusionParameterConfig(LabelParameterConfig):
    """Methionine Exclusion has some extra parameters"""

    distance: float = 7
    solvent_exposure: float = 0.35
    amino_acid: str = "MET"


class LabelingScore:
    """
    Calculation of a combined labeling score based on all labelling parameters.
    """

    folder = ""
    protein = ""
    chain_combi = []

    label_score_model: Dict[str, LabelParameterConfig] = {}
    label_params = {}
    label_parameter_score = {}
    label_score = {}

    def __init__(
        self,
        folder,
        protein,
        chain_combi,
        label_score_model: Dict[str, LabelParameterConfig],
    ):
        self.folder = folder
        self.protein = protein
        self.chain_combi = chain_combi
        self.label_score_model = label_score_model
        self.weights = self._extract_weights()
        self.para_nbr = self._count_parameter()

    def calc_labeling_scores(self):
        """
        Load calculated parameter scores, loop over all residue keys and
        calculate joined labeling score.
        """
        #            logging.info(protein)
        self.label_params = {}
        aa_keys = []
        for key in self.label_score_model:
            if self.label_score_model[key].weight == 0:
                continue
            path = os.path.join(self.folder, self.protein + r"_" + key + ".csv")

            self.label_params[key] = self._load_parameter(path)
            aa_keys.extend(self.label_params[key].keys())
        # take all occuring aa_keys in one list and sort list
        aa_keys = list(set(aa_keys))
        aa_keys.sort(key=lambda aa_key: aa_key[0] + "{:5d}".format(int(aa_key[1:])))

        #        aa_keys = self._get_aa_keys()
        label_score = {}
        label_parameter_score = {}
        for aa_key in aa_keys:
            parameters = []
            try:
                score, parameters = self._calc_labeling_score(aa_key)
            #                score = 0
            #                parameters = []
            #                para_err = False
            #                score_zero = False  # set labelling score zero if one parameter is zero
            #                for key in self.label_score_model:
            #                    if self.label_score_model[key][0] == 0:
            #                        continue
            #                    try:
            #                        if not para_err and float(label_params[key][aa_key]) >= 0:
            #                            if float(label_params[key][aa_key]) == 0:
            #                                score_zero = True
            #                            parameters.append(label_params[key][aa_key])
            #                            score += float(self.label_score_model[key][0]/weight_sum)*float(label_params[key][aa_key])
            #                        else:  # on error for one parameter: keep score -1
            #                            score = -1.0
            #                            parameters.append(label_params[key][aa_key])
            #                            para_err = True
            #                    except Exception:  # on error: set score and paramter in parameters to -1
            #                        score = -1.0
            #                        parameters.append(-1.0)
            #                        para_err = True
            #                if score_zero:
            #                    score = 0.0
            except NotImplementedError as e:
                raise e
            except Exception:
                score = -1.0  # set score to -1
                logger.debug(len(parameters))
                for _ in range(len(parameters), self.para_nbr):
                    parameters.append(
                        -1.0
                    )  # fill parameters list to normal length (with -1)
            parameters = [score] + parameters
            label_parameter_score[aa_key] = parameters  # keep all entries in long list
            if score > 0:  # TODO check if all values should be kept
                # label_parameter_score[aa_key] = parameters
                label_score[aa_key] = [score]

            self.label_parameter_score = label_parameter_score
            self.label_score = label_score

    def _load_parameter(self, path):
        with open(path, "r") as csv_file:
            reader = csv.reader(csv_file)

            parameters = dict(reader)
            if "ID" in parameters:
                parameters.pop("ID")
        return parameters

    def _extract_weights(self):
        weights = []
        for key, parameter_config in self.label_score_model.items():
            if parameter_config.weight > 0:
                weights.append(parameter_config.weight)
        return weights

    def _count_parameter(self):
        para_nbr = sum(x.weight > 0 for x in list(self.label_score_model.values()))
        return para_nbr

    def _calc_labeling_score(self, aa_key, avg="geometric"):
        score = 0
        parameters = []
        para_err = False
        score_zero = False  # set labelling score zero if one parameter is zero
        for key in self.label_score_model:
            if self.label_score_model[key].weight == 0:
                continue
            try:
                if not para_err and float(self.label_params[key][aa_key]) >= 0:
                    if float(self.label_params[key][aa_key]) == 0:
                        score_zero = True
                    parameters.append(float(self.label_params[key][aa_key]))
                #                    score += float(self.label_score_model[key][0]/sum(self.weights))*float(self.label_params[key][aa_key])
                else:  # on error for one parameter: keep score -1
                    #                    score = -1.0
                    parameters.append(float(self.label_params[key][aa_key]))
                    para_err = True
            except Exception:  # on error: set score and paramter in parameters to -1
                score = -1.0
                parameters.append(-1.0)
                para_err = True
        if score_zero:
            score = 0.0
        elif para_err:
            score = -1.0
        else:
            score = self._labeling_score(parameters, avg_type=avg)
        return score, parameters

    def _labeling_score(self, parameters, avg_type="geometric"):
        if avg_type == "arithmetic":
            ws = np.array(self.weights)
            ps = np.array(parameters)
            score = np.dot(ws, ps) / sum(ws)
        elif avg_type == "geometric":
            params = np.array(parameters)
            score = params.prod() ** (1.0 / len(params))
        elif avg_type == "harmonic":
            raise NotImplementedError(
                "Only arithmetic and geometric mean values supported "
                'avg_type needs to be "geometric" or "arithmetic"'
            )
        else:
            raise NotImplementedError(
                "Only arithmetic and geometric mean values supported "
                'avg_type needs to be "geometric" or "arithmetic"'
            )
        return score

    def save_csv(self, long=False):
        """
        Save labeling scores to file.

        :param long: Save parameter scores in addition to labeling scores.
        :type long: bool, defaults to `False`.
        """
        header = ["ID", "Labeling Score"]
        if not long:
            path = os.path.join(self.folder, self.protein + r"_LS.csv")
            score_dict = self.label_score
        else:
            path = os.path.join(self.folder, self.protein + r"_LSlong.csv")
            score_dict = self.label_parameter_score
            for key in self.label_score_model:
                if self.label_score_model[key].weight == 0:
                    continue
                header.append(key.upper())
        with open(path, "w", newline="") as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(header)

            for key, value in score_dict.items():
                line = []
                line.append(key)
                for val in value:
                    line.append(val)
                writer.writerow(line)

    def save_pdb(self):
        """
        Save labeling scores to b-factor entry in pdb file.
        """
        path = os.path.join(self.folder, self._file_ext(self.protein, "pdb"))
        parser = PDB.PDBParser(PERMISSIVE=1)
        structure = parser.get_structure(self.protein, path)
        model = structure[0]
        for chain in model:
            chain_id = chain.get_id()
            # chainCombi = [self.chainCombination, self.chainCombination2][prot_idx]
            if chain_id not in self.chain_combi:
                logger.debug("Skip chain %s", chain_id)
                for residue in chain:
                    for atom in residue:
                        atom.set_bfactor(-1)
                continue
            for residue in chain:
                # TODO ACCELERATE CHAINS WHICH ARE NOT IN LIST
                res_id = residue.get_id()
                try:
                    parameter_key = chain.get_id() + str(res_id[1])
                    b_factor = self.label_parameter_score[parameter_key][0]
                    for atom in residue:
                        atom.set_bfactor(b_factor)
                except Exception:
                    for atom in residue:
                        atom.set_bfactor(-1)

        filepath = os.path.join(self.folder, self.protein + r"_LS.pdb")
        with open(filepath, "w") as f:
            # warnings.simplefilter("ignore")
            io = PDBIO()
            io.set_structure(structure)
            io.save(f)

    @staticmethod
    def _file_ext(filename, ext):
        if ext[0] != ".":
            ext = "." + ext
        assert "." not in ext[1:], "invalid file extension"
        if len(filename) >= len(ext) and filename[-len(ext) :] == ext:
            filename_all = filename
        else:
            filename_all = filename + ext
        return filename_all
