# -*- coding: utf-8 -*-
# Copyright (C) 2020, Christian Gebhardt (chrirstian.gebhardt@bio.lmu.de)
#
# This file is part of the labelizer distribution and governed by your
# choice of the "MIT license".
# Please see the LICENSE file that should have been included as part of this
# package.

"""
The labelizer class, representing a labeling analysis with all analysis steps.
"""

import csv
import json
import logging

import os
from typing import List, Dict, Optional, Union
from zipfile import ZipFile

from . import Fluorophore
from .fret_score import FRETScore
from .labeling_parameter import LabelingParameter
from .labeling_score import (
    LabelingScore,
    LabelParameterConfig,
    AAExclusionParameterConfig,
)


# LABELLING_PARAMETER = ["SolventExposure", "ConservationScore",
#                        "TryptophanProximity", "CysteineResemblance",
#                        "SecondaryStructure", "ChargeEnvironment"]
# for lp in LABELLING_PARAMETER:  # load all subclasses of LabelingParameter
#     exec("from ."+lp+" import "+lp)


# I'm not sure if the default Labelizer.label_score_model is used anywhere so i'm rolling my own with the actual
# published defaults
label_score_model_paper: Dict[str, LabelParameterConfig] = {
    "cs": LabelParameterConfig(1, "m"),
    "se": LabelParameterConfig(1, "m"),
    "cr": LabelParameterConfig(1, "m"),
    "ss": LabelParameterConfig(1, "m"),
}


class Labelizer:
    """
    Perform labeling analysis with all analysis steps.
    """

    mode = "f"
    """Can be fret (f), single FRET (s) or single label (l)"""
    protein1 = "1OMP"  # apo
    protein2 = "1ANF"  # holo
    proteins = [protein1, protein2]
    chains1 = ["A", "A"]
    chains2 = ["A", "A"]
    chains = [chains1, chains2]
    label_score_model: Dict[str, LabelParameterConfig] = {
        "cs": LabelParameterConfig(3, "m"),
        "se": LabelParameterConfig(3, "m"),
        "cr": LabelParameterConfig(0, "m"),
        "ss": LabelParameterConfig(0, "m"),
    }
    """parameter -> (weight, sensitivity). Note that those the wrong defaults (see label_score_model_paper)"""
    parameter_score_model = {
        "cs": "N_CS2_Score",
        "se": "N_SE11_MEAN_SURFACE_DIST",
        # "I_SE4_HSE1_10A", #"N_SE1_RSA_WILKE" #N_SE11_MEAN_SURFACE_DIST
        "cr": "C_CR1_Name",
        "ss": "C_SS1_SS",
        "me": "N_ME11_Methionin_Exclusion_Dummy",
    }
    parameter_scores = {}

    fluorophore1: Union[Fluorophore, Dict[str, float]] = {
        "ll": 20.5,
        "lw": 4.5,
        "R1": 8.1,
        "R2": 4.2,
        "R3": 1.5,
    }
    fluorophore2: Union[Fluorophore, Dict[str, float]] = {
        "ll": 20.5,
        "lw": 4.5,
        "R1": 8.1,
        "R2": 4.2,
        "R3": 1.5,
    }

    # noinspection PyUnresolvedReferences
    folder = os.path.join("..", "media", "test")
    cons_serv_add = "_With_Conservation_Scores.pdb"
    prot1_cs = ""
    prot2_cs = ""
    prot_cs = [prot1_cs, prot2_cs]

    LONG_OUTPUT = False
    SAVE_PDB = False

    def __init__(
        self,
        mode: str,
        protein1: str,
        protein2: Optional[str],
        chain1: List[str],
        chain2: List[str],
        label_score_dict: Dict[str, LabelParameterConfig],
        workdir: str = os.getcwd(),
        cons_serv_add: str = "_With_Conservation_Scores.pdb",
        prot1_cs: Optional[List[str]] = None,
        prot2_cs: Optional[List[str]] = None,
        long_output: bool = False,
        save_pdb: bool = False,
    ):
        """
        Initialize the class.
        """

        self.mode = mode
        self.protein1 = self.file_ext(protein1, "pdb")[:-4]  # protein1
        if protein2 is not None:  # protein2
            self.protein2 = self.file_ext(protein2, "pdb")[:-4]
        else:
            self.protein2 = None
        self.chains1 = [c.upper() for c in chain1]
        if mode == "f":
            self.proteins = [self.protein1, self.protein2]
            self.chains2 = [c.upper() for c in chain2]
            self.chains = [self.chains1, self.chains2]
        elif mode in ("s", "l"):
            self.proteins = [self.protein1]
            self.chains = [self.chains1]
        else:
            raise NotImplementedError(
                (
                    "Only fret (f), single FRET (s),"
                    "and single label (l) modes are valid"
                )
            )
        if prot1_cs is None:
            self.prot1_cs = [
                self.protein1 + self.file_ext(cons_serv_add, "pdb")
                for _chains in chain1
            ]
        else:
            self.prot1_cs = prot1_cs
        if mode == "f":
            if prot2_cs is None:
                self.prot2_cs = [
                    self.protein2 + self.file_ext(cons_serv_add, "pdb")
                    for _chains in chain2
                ]
            else:
                self.prot2_cs = prot2_cs
            self.prot_cs = [self.prot1_cs, self.prot2_cs]
        else:
            self.prot_cs = [self.prot1_cs]
        #        self.chains1 = [c.upper() for c in chain1]  # [chain1,chain2]
        #        if mode == 'f':
        #            self.chains2 = [c.upper() for c in chain2]
        #            self.chains = [self.chains1, self.chains2]
        #        elif mode in ('s', 'l'):
        #            self.chains = [self.chains1]
        #        else:
        #            raise NotImplementedError(('Only fret (f), single FRET (s),'
        #                                       'and single label (l) modes are valid'))
        assert isinstance(label_score_dict, dict), (
            "label_score_dict is a dict of the form {'tag':(w,s)} "
            + "w as weight (float,>=0) and s as sensitivity 'h','m', or 'l'"
        )
        for key, value in label_score_dict.items():
            assert isinstance(key, str)
            assert isinstance(value, LabelParameterConfig)
        self.label_score_model = label_score_dict
        self.parameter_scores = {}
        self.folder = workdir

        self.LONG_OUTPUT = long_output
        self.SAVE_PDB = save_pdb

    # TODO change fluorophores to fluorophore class
    def set_fluorophore(
        self, fluo1: Union[Fluorophore, dict], fluo2=Optional[Union[Fluorophore, dict]]
    ):
        """
        Set fluorophores for analysis.

        :param fluo1: Fluorophore parameter of first fluorophore as dict. Must be of
                      the form {'name':'Alexa546','ll':20.5,'lw':4.5,'R1':8.1,'R2':3.2,'R3':1.5}.
        :param fluo2: Fluorophore parameter of second fluorophore for FRET. Must be of
                      the form {'name':'Alexa647','ll':20.5,'lw':4.5,'R1':8.1,'R2':3.2,'R3':1.5}.
        """
        self.fluorophore1 = fluo1
        self.fluorophore2 = fluo2

    @staticmethod
    def file_ext(filename, ext):
        if ext[0] != ".":
            ext = "." + ext
        assert "." not in ext[1:], "invalid file extension"
        if len(filename) >= len(ext) and filename[-len(ext) :] == ext:
            filename_all = filename
        else:
            filename_all = filename + ext
        return filename_all

    def calc_parameter_score(self):
        """
        Calculate labelling parameter scores for selected parameters.
        """
        # Avoid circular import
        from labelizer import AAExclusion

        for tag, parameter_config in self.label_score_model.items():
            if parameter_config.weight > 0:
                for protein in self.proteins:
                    labelingParameter = LabelingParameter.factory(tag)
                    if tag == AAExclusion.file_tag:
                        assert isinstance(parameter_config, AAExclusionParameterConfig)
                        assert isinstance(labelingParameter, AAExclusion)
                        labelingParameter.max_bad_aa_distance = (
                            parameter_config.distance
                        )
                        labelingParameter.min_bad_aa_solvent_exposure = (
                            parameter_config.solvent_exposure
                        )
                        labelingParameter.bad_aa = parameter_config.amino_acid
                    else:
                        assert not isinstance(
                            parameter_config, AAExclusionParameterConfig
                        )
                    labelingParameter.set_up(
                        self,
                        protein,
                        parameter_config.sensitivity,
                        self.parameter_score_model[tag],
                    )
                    labelingParameter.calc_parameter_scores()
                    labelingParameter.clean_up(save_pdb=True)

    def calc_labeling_score(self):
        """
        Calculate labeling scores based on parameter scores for all structures.
        """
        for prot_idx, protein in enumerate(self.proteins):
            chain_combi = [self.chains1, self.chains2][prot_idx]
            ll_score = LabelingScore(
                self.folder, protein, chain_combi, self.label_score_model
            )
            ll_score.calc_labeling_scores()
            if self.SAVE_PDB:
                ll_score.save_pdb()  # save pdb with scores in b-factor field
            ll_score.save_csv()
            ll_score.save_csv(long=True)

    def calc_fret_score(
        self, foerster_radius=60, method="SIMPLE", off=8.0, save_AV=False, weights=True
    ):
        """
        Calculate FRET scores based on labeling scores, fluorophore distance, and Foerster radiues.

        :param foerster_radius: Foerster radius in Angstrom, defaults to 60.
        :type foerster_radius: float, optional
        :param method: Method to evaluate distances.
                       Must be "SIMPLE", "GEBHARDT", or "KALININ",
                       defaults to "SIMPLE"
        :type method: string, optional
        :param off: Offset for mean dye position from protein surface
                    (only for method "SIMPLE"), defaults to 8.
        :type off: float, optional
        :param save_AV: `True` if accessible volume clouds should be saved
                        (only for method "KALININ"), defaults to `False`
        :type save_AV: boolean, optional
        :param weights: I DON'T KNOW
        :type weights: I DON'T KNOW, optional
        """

        if self.mode == "f":
            self._calc_fret_score_double(
                foerster_radius=foerster_radius,
                method=method,
                off=off,
                save_AV=save_AV,
                weights=weights,
            )
        elif self.mode == "s":  # single FRET (only one pdb file)
            self._calc_fret_score_single(
                foerster_radius=foerster_radius,
                method=method,
                off=off,
                save_AV=save_AV,
                weights=weights,
            )
        else:
            pass

    def _calc_fret_score_double(
        self,
        foerster_radius="auto",
        method="SIMPLE",
        off=8.0,
        save_AV=False,
        weights=True,
    ):
        """
        Calculate FRET score (two conformations)
        """
        fret = FRETScore()
        fret.set_method(method)
        fret.set_offset_AV(off)
        fret.set_save_AV(save_AV)
        fret.set_fluorophore(self.fluorophore1, self.fluorophore2)
        fret.set_foerster_radius(foerster_radius)
        fret.set_folder(self.folder)
        fret.set_weights(weights)

        # load calculated labeling scores
        apo_path = os.path.join(self.folder, self.protein1 + "_LS.csv")
        holo_path = os.path.join(self.folder, self.protein2 + "_LS.csv")
        fret.load_labeling_score(apo_path, self.protein1, "apo")
        fret.load_labeling_score(holo_path, self.protein2, "holo")
        # load pdbs
        apo_pdb_path = os.path.join(self.folder, self.file_ext(self.protein1, "pdb"))
        holo_pdb_path = os.path.join(self.folder, self.file_ext(self.protein2, "pdb"))
        fret.load_pdb(apo_pdb_path, self.protein1, "apo")
        fret.load_pdb(holo_pdb_path, self.protein2, "holo")
        # calculate pairing score
        fret.calc_measurement_scores(self.chains1, self.chains2)
        # export results
        fret.save_csv()
        fret.save_csv("long")

        try:
            fret.save_csv("refine")
        except Exception as ex:
            logging.debug(type(ex).__name__, str(ex))
            pass

        try:
            filepath = os.path.join(
                self.folder, self.protein1 + "_" + self.protein2 + "_heatmap.json"
            )
            with open(filepath, "w") as file:
                json.dump(fret.heatmap_data, file, indent=2)
        except Exception:
            logging.info("Failed to save heatmap.")
            pass

    def _calc_fret_score_single(
        self,
        foerster_radius="auto",
        method="SIMPLE",
        off=8.0,
        save_AV=False,
        weights=True,
    ):
        """
        Calculate FRET score (single conformation)
        """
        fret = FRETScore()
        fret.set_method(method)
        fret.set_offset_AV(off)
        fret.set_save_AV(save_AV)
        fret.set_fluorophore(self.fluorophore1, self.fluorophore2)
        fret.set_foerster_radius(foerster_radius)
        fret.set_folder(self.folder)
        fret.set_weights(weights)

        # load calculated labeling scores
        path = os.path.join(self.folder, self.protein1 + "_LS.csv")
        fret.load_labeling_score(path, self.protein1, "apo")
        # load pdbs
        pdb_path = os.path.join(self.folder, self.file_ext(self.protein1, "pdb"))
        fret.load_pdb(pdb_path, self.protein1, "apo")
        # calculate pairing score
        fret.calc_measurement_scores_single(self.chains1)
        # export results
        fret.save_csv()
        fret.save_csv("long")

        try:
            fret.save_csv("refine")
        except:
            pass

    # FUTURE def calc_pife_score(self)
    #    """
    #    Calculate PIFE score
    #    """

    # FUTURE def calc_quenching_score(self)
    #    """
    #    Calculate quenching score
    #    """

    def get_all_file_paths(self, directory):
        """
        Search all files (not subdirectories) in folder and return paths of files.

        :param directory: directory to search files in.
        :type directory: string
        :return: list of file paths
        :rtype: List[string]
        """
        # initializing empty file paths list
        file_paths = []

        # crawling through directory and subdirectories
        for root, directories, files in os.walk(directory):
            for filename in files:
                # join the two strings in order to form the full filepath.
                filepath = os.path.join(root, filename)
                file_paths.append(filepath)

        # returning all file paths
        return file_paths

    def zip_files(self):
        """
        Combine all evaluation files (imput and output) and zip files to archive.

        :return: path or zip-file
        :rtype: string
        """
        # path to folder which needs to be zipped
        # calling function to get all file paths in the directory
        file_paths = self.get_all_file_paths(self.folder)
        logging.debug("%s files will be zipped:", len(file_paths))

        # writing files to a zipfile
        zip_path = os.path.join(self.folder, "result.zip")
        with ZipFile(zip_path, "w") as zip_:
            # writing each file one by one
            for file in file_paths:
                if file != zip_path:
                    head_tail = os.path.split(file)
                    zip_.write(file, head_tail[1])
                else:
                    logging.debug("Skip file")
            zip_.close()
        logging.info("All files zipped successfully!")
        return zip_path

    # helper function for testing
    # noinspection PyUnresolvedReferences
    def load_parameter(self, path):
        """
        Helper function to load parameter for testing class.

        :param path: file path to parameter score csv-file.
        :type path: string
        """
        with open(path, "r") as csv_file:
            reader = csv.reader(csv_file)
            label_params = dict(reader)
            if "ID" in label_params:
                label_params.pop("ID")
            # convert values to float
            for key in label_params.keys():
                label_params[key] = float(label_params[key])
        return label_params
