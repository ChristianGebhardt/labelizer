# -*- coding: utf-8 -*-
# Copyright (C) 2020, Christian Gebhardt (chrirstian.gebhardt@bio.lmu.de)
#
# This file is part of the labelizer distribution and governed by your
# choice of the "MIT license".
# Please see the LICENSE file that should have been included as part of this
# package.

"""
The module provides an abstract class for all assay scores, which inherit from it.
"""
import csv
import logging
import os
import time
from math import pi
from os.path import join as pthjoin

import numpy as np
import pandas as pd
from Bio.PDB import NeighborSearch
from Bio.PDB.vectors import rotaxis

from .config import cfg as ll_cfg
from .fluorophore import Fluorophore
from .measurement_score import MeasurementScore

try:
    # try to import LabelLib
    from . import label_lib_functions as llf
except ImportError:
    llf = None


# TBD
# * flag for results (e.g. extendended coord removal, failed, ...)


class FRETScore(MeasurementScore):
    """
    Class implements the FRET assay analysis as subclass of
    abstract class MeasurementScore.

    FRETScore includes method to simulate residue / accessible volume distances
    based on different methods and calculate a FRET assay score for a single
    conformation or two conformations.
    """

    # constants / parameter
    file_tag = "FRET"
    method_options = ["SIMPLE", "GEBHARDT", "KALININ"]
    method = "SIMPLE"  # SIMPLE, GEBHARDT or KALININ

    SIMPLE_AV = {}
    SIMPLE_AV_F2 = {}
    #    AV_CALC_RADIUS = 20.0
    off_alpha = "auto"
    save_av = False
    folder = None
    use_weights = True
    json_heatmap = None

    #    COLLISION_RADIUS = 1.7 # Used for carbon: https://de.wikipedia.org/wiki/Van-der-Waals-Radius

    DEBUG_FRET_SCORE = ll_cfg.DEBUG
    #    SAVE_DISTANCE_DATA = True

    #    KALININ_SPEED = ["fast", "precise"][1]
    #    N_KALININ = {"fast": 10000, "precise": 100000}
    #    GS = {"fast": 1.2, "precise": 0.8}
    # dye parameter: linker length, linker width radius 1, radius 2, radius 3
    fluo_parameter = [14.0, 4.5, 5.5, 5.5, 2.0]  # 23.5, 4.5, 8.1, 4.2, 1.5
    fluo_parameter2 = [14.0, 4.5, 5.5, 5.5, 2.0]

    # off = 5.2
    #    LS_THRESHOLD = 0.5

    def __init__(self):
        super(FRETScore, self).__init__()
        #        self.KALININ_SPEED = ll_cfg.SPEED
        if ll_cfg.SPEED == "fast":
            self.N_KALININ = ll_cfg.N_KALININ_FAST
            self.GS = ll_cfg.GS_FAST
        else:
            self.N_KALININ = ll_cfg.N_KALININ_PRECISE
            self.GS = ll_cfg.GS_PRECISE
        self.foerster_radius = "auto"
        self.donor = None
        self.acceptor = None
        self.measurement_score_header = ["ID", "FRET Score"]
        self.measurement_score_refine = None

    def set_foerster_radius(self, radius="auto"):
        """
        Set Förster radius for analysis.

        :param radius: Förster radius (in Angström).
        :type radius: Union[float, str]
        """
        if radius != "auto":
            self.foerster_radius = radius
        elif type(self.donor) == Fluorophore and type(self.acceptor) == Fluorophore:
            self.foerster_radius = self.donor.calc_foerster_radius(self.acceptor)
        else:
            raise TypeError(
                "Either specify a radius (float) or set attributes"
                "donor and acceptor to type Fluorophore"
            )

    def set_method(self, method):
        """
        Set calculation method for fluorophore distance estimation.

        :param method: Distance calculation method.
                       Must be "SIMPLE", "GEBHARDT", or "KALININ",
                       defaults to "SIMPLE"
        """
        if method in self.method_options:
            self.method = method
            # test import DELETE
            if method == "KALININ" and llf is not None:
                raise ImportError("AV module LabelLib not installed!")
        else:
            raise ValueError("method must be one of " + str(self.method_options))

    def set_offset_AV(self, offset):
        """Set offset of alpha-sphere for distance simulation."""
        self.off_alpha = offset

    def set_save_AV(self, save_av):
        """Set save accessible volumes."""
        self.save_av = save_av

    def set_folder(self, folder):
        """Set analysis folder."""
        self.folder = folder

    def set_weights(self, weights):
        """Set usage of weights in accessible volume simulation."""
        self.use_weights = weights

    def set_fluorophore(self, fluo1, fluo2):
        """
        Set fluorophores for analysis.

        :param fluo1: Fluorophore parameter of first fluorophore as dict.
                      Must be of the form
                      {'name':'Alexa546','ll':20.5,'lw':4.5,'R1':8.1,'R2':3.2,'R3':1.5}.
        :param fluo2: Fluorophore parameter of first fluorophore as dict.
                      Must be of the form
                      {'name':'Alexa647','ll':20.5,'lw':4.5,'R1':11.0,'R2':3.2,'R3':1.5}.
        """
        # backward compatibility (for fluorophores, which are only a dict)
        if type(fluo1) == Fluorophore:
            self.donor = fluo1
            fluo1 = fluo1.fluo_para
            self.fluo_parameter = [
                fluo1["ll"],
                fluo1["lw"],
                fluo1["r1"],
                fluo1["r2"],
                fluo1["r3"],
            ]
        else:
            self.fluo_parameter = [
                fluo1["ll"],
                fluo1["lw"],
                fluo1["R1"],
                fluo1["R2"],
                fluo1["R3"],
            ]
        if type(fluo2) == Fluorophore:
            self.acceptor = fluo2
            fluo2 = fluo2.fluo_para
            self.fluo_parameter2 = [
                fluo2["ll"],
                fluo2["lw"],
                fluo2["r1"],
                fluo2["r2"],
                fluo2["r3"],
            ]
        else:
            self.fluo_parameter2 = [
                fluo2["ll"],
                fluo2["lw"],
                fluo2["R1"],
                fluo2["R2"],
                fluo2["R3"],
            ]

    @staticmethod
    def joined_label_score(*scores):
        """
        Calculation of a joined label score jls. E.g. two residues for single FRET pair,
        or 4 residues for FRET pair in two conformational states.

        :math:`jls = N \\frac{ s_1 * s_2 * ... * s_N }{ s_2*s_3*...*s_N + s_1*s_3*...*s_N + ... + s_1*s_2*...*s_{N-1}}`

        :param scores: Labeling scores for all residues contributing to
                       the FRET assay.
                       Variable length score list.
        :type scores: List[float]
        :return: joined label score jls.
        :rtype: float
        """
        # normalization = 0.
        # product = 1.
        # for idx, score in enumerate(scores):
        #     product *= score
        #     reduced_product = 1.
        #     for red_score in scores[:idx]+scores[idx+1:]:
        #         reduced_product *= red_score
        #     normalization += reduced_product
        #     #python 3.8
        #     # normalization += math.prod([scores[:idx]+scores[idx+1:])

        params = np.array(scores)
        # OLD
        # score = params.prod()**(1.0/len(params))
        score = params.prod() ** (0.5)

        # return len(scores)*product/normalization
        return score

    def calc_measurement_score(self, label_scores, apo_dist, holo_dist, negative=False):
        """
        Calculation of measurement score ms for FRET assays with two conformations
        based on joined labeling score jls and FRET efficiencies :math:`E_i`
        for the two fluorophore distances.

        :math:`ms = jls * | E_1 -E_2 |`

        :param label_scores: Labeling scores for all residues contributing to
                       the FRET assay. Length of list should be 4.
        :type label_scores: List[float]
        :param apo_dist: fluorophore distance in first conformational state.
        :type apo_dist: float
        :param holo_dist: fluorophore distance in second conformational state.
        :type holo_dist: float
        :param negative: Calculate best score for negative mutants without distance change.
        :type negative: boolean
        :return: calculated measurement score.
        :rtype: float
        """
        apo_efficiency = self._fret_efficiency(apo_dist)
        holo_efficiency = self._fret_efficiency(holo_dist)
        #        scoring = self.collective_score(ls_apo, ls_apo2, ls_holo, ls_holo2)
        scoring = self.joined_label_score(*label_scores)
        if not negative:
            val = scoring * abs(apo_efficiency - holo_efficiency)
        else:
            val = (
                scoring
                * (1 - abs(1 - (apo_efficiency + holo_efficiency)))
                * max(0, (1 - 20 * abs(apo_efficiency - holo_efficiency)))
            )
        # assert val <= 1, "measurement score larger 1"
        return val

    def _fret_efficiency(self, radius):
        return 1 / (1 + (radius / self.foerster_radius) ** 6)

    def calc_measurement_score_single(self, label_scores, dist):
        """
        Calculation of measurement score ms for FRET assays with one conformation
        based on joined labeling score jls and FRET efficiency E for the
        fluorophore distance.

        :math:`ms = jls * ( 1 -| 2 * E - 1 | )`

        :param label_scores: Labeling scores for all residues contributing to
                       the FRET assay. Length of list should be 2.
        :type label_scores: List[float]
        :param dist: fluorophore distance based on simulation.
        :type dist: float
        :return: calculated measurement score.
        :rtype: float
        """
        efficiency = self._fret_efficiency(dist)
        #        scoring = self.collective_score_single(ls1, ls2)
        try:
            scoring = self.joined_label_score(*label_scores)
        except TypeError as exc:
            logging.debug(
                label_scores, label_scores[0], type(label_scores), type(label_scores[0])
            )
            raise exc
        val = scoring * (1 - 2 * abs(efficiency - 0.5))
        # assert val <= 1, "measurement score larger 1"
        return val

    @staticmethod
    def _get_gly_cb_vector(residue):
        """Return a pseudo CB vector for a Gly residue (private).
        The pseudoCB vector is centered at the origin.
        CB coord=N coord rotated over -120 degrees
        along the CA-C axis.

        Function taken from: Bio/PDB/HSExposure.py
        """
        try:
            n_v = residue["N"].get_vector()
            c_v = residue["C"].get_vector()
            ca_v = residue["CA"].get_vector()
        except KeyError:
            return None
        # center at origin
        n_v = n_v - ca_v
        c_v = c_v - ca_v
        # rotation around c-ca over -120 deg
        rot = rotaxis(-pi * 120.0 / 180.0, c_v)
        cb_at_origin_v = n_v.left_multiply(rot)
        # move back to ca position
        # cb_v = cb_at_origin_v + ca_v
        return cb_at_origin_v

    def _get_cb_position(self, residue):
        if residue.get_resname() != "GLY":
            cb_position = residue["CB"].get_coord()
        else:
            cb_position = self._get_gly_cb_vector(residue)
            cb_position = residue["CA"].get_coord() + np.array(list(cb_position))
        cb_position = np.array(cb_position).astype(float)
        return cb_position

    def _calc_simple_mean_position(self, coords, cb_pos, fluo_para):
        mean_atoms = np.mean(coords, axis=0)
        d = np.linalg.norm(cb_pos - mean_atoms)
        R = ll_cfg.AV_CALC_RADIUS
        ll = fluo_para[0]
        radii = fluo_para[2:]
        if self.off_alpha == "auto":
            b = max(ll_cfg.COLLISION_RADIUS, 2 * min(radii) - ll_cfg.COLLISION_RADIUS)
            off = b + 0.54 * ll - 0.0225 * ll**2
        else:
            off = self.off_alpha
        alpha = np.arccos(1 - 8 / 3 * d / R)
        norm_vector = (cb_pos - mean_atoms) / np.linalg.norm(cb_pos - mean_atoms)
        mean_av = cb_pos + (3.0 / 4.0 - d / R) * (ll + off) * norm_vector

        return mean_av, alpha, off

    def _save_simple_av(self, tag, fluo_parameter, off):
        if "holo" in self.identifier:
            folder_name = (
                self.identifier["apo"]
                + "_"
                + self.identifier["holo"]
                + "_radii_{}_{}_{}".format(*fluo_parameter[2:5])
            )
        else:
            folder_name = self.identifier["apo"] + "_radii_{}_{}_{}".format(
                *fluo_parameter[2:5]
            )
        directory = pthjoin(os.path.dirname(self.file_path), folder_name)
        if not os.path.exists(directory):
            os.makedirs(directory)
        if self.off_alpha == "auto":
            filepath = pthjoin(
                directory,
                self.identifier[tag]
                + r"_AV_GEBHARDT_ll{0:.1f}_offAUTO_.csv".format(fluo_parameter[0]),
            )
        else:
            filepath = pthjoin(
                directory,
                self.identifier[tag]
                + r"_AV_GEBHARDT_ll{0:.1f}_off{1:.1f}_.csv".format(
                    fluo_parameter[0], off
                ),
            )
        with open(filepath, "w", newline="") as csv_file:
            writer = csv.writer(csv_file)
            for key, value in self.SIMPLE_AV[tag].items():
                line = []
                line.append(key)
                for val in value:
                    line.append(val)
                writer.writerow(line)

    def _generate_coords(self, model):
        coords = []
        for chain in model:
            for residue in chain:

                # TBD handle other non-amino-acid molecules
                if "CA" not in residue:
                    continue

                for atom in residue:
                    vec = atom.get_vector()
                    coords.append([*vec, ll_cfg.COLLISION_RADIUS])
        logging.debug("Length coords: %d", len(coords))
        return coords

    def _reduce_coordinates(self, coords, cb_pos, lw):
        vdw = ll_cfg.COLLISION_RADIUS
        gs = self.GS  # self.GS[self.KALININ_SPEED]
        coords_res = []
        for c in coords:
            if np.linalg.norm(np.array(c[:3]) - cb_pos) >= lw / 2.0 + vdw + gs:
                coords_res.append(c)
        return coords_res

    def _calc_single_av(
        self, coords, tag, key, fluo_parameter, save=False, filename=""
    ):
        chain = key[0]
        res_nbr = int(key[1:])
        res = self.model[tag][chain][res_nbr]
        cb_position = self._get_cb_position(res)
        coords_res = self._reduce_coordinates(coords, cb_position, fluo_parameter[1])
        if save:
            av = llf.calcAV(
                coords_res,
                cb_position,
                linker_length=fluo_parameter[0],
                linker_width=fluo_parameter[1],
                dye_radii=fluo_parameter[2:],
                grid_size=self.GS,  # self.GS[self.KALININ_SPEED],
                save=True,
                filename=filename,
            )
        else:
            av = llf.calcAV(
                coords_res,
                cb_position,
                linker_length=fluo_parameter[0],
                linker_width=fluo_parameter[1],
                dye_radii=fluo_parameter[2:],
                grid_size=self.GS,
            )  # self.GS[self.KALININ_SPEED])
        if not self.use_weights:
            av = llf.removeWeights(av)
        return av

    #    def calc_and_save_avs(self, tag, fluo_parameter):
    #        av_folder = pthjoin(self.folder, 'AVs')
    #        if not (os.path.exists(av_folder) and os.path.isdir(av_folder)):
    #            os.mkdir(av_folder)
    #        folder_name = self.identifier[tag] + '_AV' + \
    #                      '_ll_{}_radii_{}_{}_{}'.format(fluo_parameter[0],
    #                                                     *fluo_parameter[2:])
    #        temp_path = pthjoin(self.folder, 'AVs', folder_name)
    #        if not (os.path.exists(temp_path) and os.path.isdir(temp_path)):
    #            os.mkdir(temp_path)
    #
    #        coords = self._generate_coords(self.model[tag])
    #
    #        for key, ls in self.labeling_scores['apo'].items():
    #            try:
    #                filename = pthjoin(temp_path, key+'.grid')
    #                av = self._calc_single_av(coords, tag, key, fluo_parameter,
    #                                          save=True, filename=filename)
    #                filename_xyz = filename[:-4] + "xyz"
    #                llf.saveXYZ(filename_xyz, av)
    #            except:
    #                if self.DEBUG_FRET_SCORE:
    #                    raise
    #                logging.debug("Failed to generate AV %s", key)

    def _calc_or_load_av(self, tag, key, fluorophore, load=True, coords=None):
        if load:  # TBD load with fluorophore
            search_path = pthjoin(self.folder, tag + "_AV", key + ".grid")
            if os.path.exists(search_path):
                av = llf.loadAV(search_path)
            else:
                logging.debug("File not found: %s", search_path)
                av = None
        else:
            try:
                av = self._calc_single_av(coords, tag, key, fluorophore, save=False)
            except:
                if self.DEBUG_FRET_SCORE:
                    raise
                logging.debug("Failed to generate AV %s", key)
                av = None
        return av

    def calc_measurement_scores(self, chains_apo, chains_holo, N_refine=300):
        """
        Calculation of measurement scores for FRET assays with two conformations.
        Labeling scores for residues in the specified chains are loaded and
        FRET scores are calculated for all pairs of residues therein.

        :param chains_apo: Chain identifiers for first conformation. One
                           fluorophore is labelled in first chain (chains_apo[0]) and
                           the other fluorphore is labelled in second chain.
        :type chains_apo:  List[string]
        :param chains_holo: Chain identifiers for second conformation similar to chains_apo.
        :type chains_holo:  List[string]
        :param N_refine: Number of pairs to be refined with Kalinin-method
                         after fast screening.
        :type N_refine:  integer
        """
        self.measurement_score_long_header = [
            "ID",
            "FRET Score",
            "Distance 1",
            "Distance 2",
            "Mean Labeling Score",
            "Labeling Score 1 (Res. 1)",
            "Labeling Score 1 (Res. 2)",
            "Labeling Score 2 (Res. 1)",
            "Labeling Score 2 (Res. 2)",
        ]

        ls_dist_apo = self._calc_distances_single("apo", chains_apo)
        ls_dist_holo = self._calc_distances_single("holo", chains_holo)
        self.measurement_score = {}
        self.measurement_score_long = {}
        for ms_key, ls_dist in ls_dist_apo.items():
            ms_key_holo = (
                chains_holo[0]
                + ms_key.split("_")[0][1:]
                + "_"
                + chains_holo[1]
                + ms_key.split("_")[1][1:]
            )
            #            print(ms_key_holo)
            if ms_key_holo in ls_dist_holo:
                ls_dist2 = ls_dist_holo[ms_key_holo]
                ls_all = ls_dist[:2] + ls_dist2[:2]
                ms = self.calc_measurement_score(ls_all, ls_dist[2], ls_dist2[2])
                self.calc_measurement_score(
                    ls_all, ls_dist[2], ls_dist2[2], negative=True
                )
                self.measurement_score[ms_key] = [ms]
                # self.measurement_score_long[ms_key] = [ms, ls_dist[2], ls_dist2[2], self.joined_label_score(*ls_all), *ls_all, ms_neg]
                self.measurement_score_long[ms_key] = [
                    ms,
                    ls_dist[2],
                    ls_dist2[2],
                    self.joined_label_score(*ls_all),
                    *ls_all,
                ]

        if self.method != "KALININ" and N_refine > 0:
            try:
                # self.measurement_score_long_header.extend(["Refined FRET Score",
                #                                           "Refined Distance 1",
                #                                           "Refined Distance 2"])
                df = pd.DataFrame.from_dict(self.measurement_score_long, orient="index")
                df = df.sort_values(by=df.columns[0], ascending=False)
                refinement_list = list(df.index.values)[:N_refine]
                ms_refinement_apo = self._distances_refinement(
                    "apo", chains_apo, refinement_list
                )
                ms_refinement_holo = self._distances_refinement(
                    "holo", chains_holo, refinement_list
                )
                self.measurement_score_refine = {}
                for ms_key in refinement_list:
                    try:
                        dist_apo = ms_refinement_apo[ms_key]
                        dist_holo = ms_refinement_holo[ms_key]
                        ls_all = self.measurement_score_long[ms_key][4:8]
                        ms = self.calc_measurement_score(ls_all, dist_apo, dist_holo)
                        self.measurement_score_refine[ms_key] = [
                            ms,
                            dist_apo,
                            dist_holo,
                            self.joined_label_score(*ls_all),
                            *ls_all,
                        ]  # .extend([ms, dist_apo, dist_holo])
                    except Exception as ex:
                        # if self.DEBUG_FRET_SCORE:
                        #    raise
                        logging.debug(type(ex).__name__, ms_key)
            except Exception as ex:
                logging.info(type(ex).__name__, "refinement calculation failed")
        try:
            self.json_heatmap = self._calc_heat_map(chains_apo, chains_holo)
        except:
            logging.info("Failed to calculate heatmap.")

    def calc_measurement_scores_single(self, chains, N_refine=300):
        """
        Calculation of measurement scores for FRET assays with one conformation.
        Labeling scores for residues in the specified chains are loaded and
        FRET scores are calculated for all pairs of residues therein.

        :param chains: Chain identifiers for first conformation. One
                       fluorophore is labelled in first chain (chains[0]) and
                       the other fluorphore is labelled in second chain (chains[1]).
        :type chains:  List[string]
        :param N_refine: Number of pairs to be refined with Kalinin-method
                         after fast screening.
        :type N_refine:  integer
        """
        self.measurement_score_long_header = [
            "ID",
            "FRET Score",
            "Distance",
            "Mean Labeling Score",
            "Labeling Score (Res. 1)",
            "Labeling Score (Res. 2)",
        ]
        ls_distances = self._calc_distances_single("apo", chains)
        self.measurement_score = {}
        self.measurement_score_long = {}
        for ms_key, ls_dist in ls_distances.items():
            ms = self.calc_measurement_score_single(ls_dist[:2], ls_dist[2])
            self.measurement_score[ms_key] = [ms]
            self.measurement_score_long[ms_key] = [
                ms,
                ls_dist[2],
                self.joined_label_score(*(ls_dist[:2])),
                *(ls_dist[:2]),
            ]

        if self.method != "KALININ" and N_refine > 0:
            try:
                # self.measurement_score_long_header.extend(["Refined FRET Score",
                #                                           "Refined Distance"])
                df = pd.DataFrame.from_dict(self.measurement_score_long, orient="index")
                df = df.sort_values(by=df.columns[0], ascending=False)
                refinement_list = list(df.index.values)[:N_refine]
                ms_refinement = self._distances_refinement(
                    "apo", chains, refinement_list
                )
                self.measurement_score_refine = {}
                for ms_key in refinement_list:
                    try:
                        dist = ms_refinement[ms_key]
                        ms = self.calc_measurement_score_single(
                            self.measurement_score_long[ms_key][3:], dist
                        )
                        self.measurement_score_refine[ms_key] = [
                            ms,
                            dist,
                            *self.measurement_score_long[ms_key][2:],
                        ]
                        # self.measurement_score_long[ms_key] = [ms,ls_dist[2], self.joined_label_score(*(ls_dist[:2])), *(ls_dist[:2])]
                    except Exception as ex:
                        if self.DEBUG_FRET_SCORE:
                            raise
                        logging.debug(type(ex).__name__, ms_key)
            except Exception as ex:
                logging.info(type(ex).__name__, "refinement calculation failed")

    def _calc_distances_single(self, tag, chain1):
        logging.info("FRET score %s", chain1)
        ls_distance_data = {}

        # set up
        if self.method == "GEBHARDT":
            atom_list = list(self.model[tag].get_atoms())
            nb_search = NeighborSearch(atom_list)
        elif self.method == "KALININ":
            coords = self._generate_coords(self.model[tag])

        for key, ls in self.labeling_scores[tag].items():
            if self.method == "KALININ":  # DELETE
                logging.debug("KEY %s", key)
            # skip small labelling scores
            if float(ls) < ll_cfg.LS_THRESHOLD:
                continue
            start = time.time()
            len_pos = 0  # for timing
            res = self.model[tag][key[0]][int(key[1:])]
            try:
                if self.method == "SIMPLE":
                    av = self._get_cb_position(res)
                    av_f2 = av
                elif self.method == "GEBHARDT":
                    cb_position = self._get_cb_position(res)
                    # defined for `self.method == "GEBHARDT"`
                    # noinspection PyUnboundLocalVariable
                    search_result = nb_search.search(
                        cb_position, ll_cfg.AV_CALC_RADIUS, level="A"
                    )
                    coords_search = []
                    for atom in search_result:
                        coords_search.append(atom.get_coord())
                    # av, alpha, offset
                    av, _, _ = self._calc_simple_mean_position(
                        coords_search, cb_position, self.fluo_parameter
                    )
                    av_f2, _, _ = self._calc_simple_mean_position(
                        coords_search, cb_position, self.fluo_parameter2
                    )
                elif self.method == "KALININ":
                    if self.save_av and False:
                        av = self._calc_or_load_av(
                            tag, key, self.fluo_parameter, load=True
                        )
                        if self.fluo_parameter != self.fluo_parameter2:
                            av_f2 = self._calc_or_load_av(
                                tag, key, self.fluo_parameter2, load=True
                            )
                    else:
                        # defined for `self.method == "KALININ"`
                        # noinspection PyUnboundLocalVariable
                        av = self._calc_or_load_av(
                            tag, key, self.fluo_parameter, load=False, coords=coords
                        )
                        if self.fluo_parameter != self.fluo_parameter2:
                            av_f2 = self._calc_or_load_av(
                                tag,
                                key,
                                self.fluo_parameter2,
                                load=False,
                                coords=coords,
                            )
                    if self.fluo_parameter == self.fluo_parameter2:
                        av_f2 = av

                    if len(av.points()[0]) < 1000:
                        logging.debug("%s AV (D) size %d", key, len(av.points()[0]))
                        continue
                    if len(av_f2.points()[0]) < 1000:
                        logging.debug("%s AV (A) size %d", key, len(av_f2.points()[0]))
                        continue
            except:
                if self.DEBUG_FRET_SCORE:
                    raise
                av = None
                av_f2 = None

            for key2, ls2 in self.labeling_scores[tag].items():
                if key[0] == key2[0] and int(key[1:]) >= int(key2[1:]):
                    continue
                if float(ls2) < ll_cfg.LS_THRESHOLD:
                    continue
                ms_key = str(key) + "_" + str(key2)
                len_pos += 1
                res2 = self.model[tag][key2[0]][int(key2[1:])]
                if key != key2 and (chain1[0] == key[0] and chain1[1] == key2[0]):
                    if self.method == "SIMPLE":
                        av2 = self._get_cb_position(res2)
                        av2_f2 = av2
                    elif self.method == "GEBHARDT":
                        try:
                            cb_position = self._get_cb_position(res2)
                            search_result = nb_search.search(
                                cb_position, ll_cfg.AV_CALC_RADIUS, level="A"
                            )
                            coords_search = []
                            for atom in search_result:
                                coords_search.append(atom.get_coord())
                            av2, _, _ = self._calc_simple_mean_position(
                                coords_search, cb_position, self.fluo_parameter
                            )
                            av2_f2, _, _ = self._calc_simple_mean_position(
                                coords_search, cb_position, self.fluo_parameter2
                            )
                        except:
                            if self.DEBUG_FRET_SCORE:
                                raise
                            av2 = None
                            av2_f2 = None
                    elif self.method == "KALININ":
                        try:
                            if self.save_av and False:
                                av2 = self._calc_or_load_av(
                                    tag, key2, self.fluo_parameter, load=True
                                )
                                if self.fluo_parameter != self.fluo_parameter2:
                                    av2_f2 = self._calc_or_load_av(
                                        tag, key2, self.fluo_parameter2, load=True
                                    )
                            else:
                                av2 = self._calc_or_load_av(
                                    tag,
                                    key2,
                                    self.fluo_parameter,
                                    load=False,
                                    coords=coords,
                                )
                                if self.fluo_parameter != self.fluo_parameter2:
                                    av2_f2 = self._calc_or_load_av(
                                        tag,
                                        key2,
                                        self.fluo_parameter2,
                                        load=False,
                                        coords=coords,
                                    )
                            if self.fluo_parameter == self.fluo_parameter2:
                                av2_f2 = av2

                            if len(av2.points()[0]) < 1000:
                                logging.debug(
                                    "%s AV2 (D) size %d", key2, len(av2.points()[0])
                                )
                                continue
                            if len(av2_f2.points()[0]) < 1000:
                                logging.debug(
                                    "%s AV2 (A) size %d", key2, len(av2_f2.points()[0])
                                )
                                continue
                        except:
                            if self.DEBUG_FRET_SCORE:
                                raise
                            av2 = None
                            av2_f2 = None
                    else:
                        raise NotImplementedError(
                            "Method " + str(self.method) + " not implemented!"
                        )
                    try:
                        if self.method in ("SIMPLE", "GEBHARDT"):
                            dist = np.linalg.norm(np.array(av) - np.array(av2_f2))
                            dist2 = np.linalg.norm(np.array(av_f2) - np.array(av2))
                        #                            dist = np.linalg.norm(np.array(av)-np.array(av2_f2))
                        #                            dist2 = dist
                        else:
                            n_avg = self.N_KALININ  # self.N_KALININ[self.KALININ_SPEED]
                            # dist = llf.rdaDistance(av, av2_f2, n=n_avg)
                            dist = llf.effDistance(
                                av, av2_f2, self.foerster_radius, n=n_avg
                            )
                            #                            dist = llf.rmpDistance(av, av2_f2)
                            #                            dist2 = dist
                            if self.fluo_parameter != self.fluo_parameter2:
                                # dist2 = llf.rdaDistance(av_f2, av2, n=n_avg)
                                dist2 = llf.effDistance(
                                    av_f2, av2, self.foerster_radius, n=n_avg
                                )
                            else:
                                dist2 = dist
                        ls_distance_data[ms_key] = [
                            float(ls),
                            float(ls2),
                            round((dist + dist2) / 2.0, 1),
                        ]
                    except Exception as ex:
                        if self.DEBUG_FRET_SCORE:
                            raise
                        logging.debug(type(ex).__name__, key, key2)

            end = time.time()
            if self.method == "KALININ":
                logging.info("Time for %d steps: %.1f s", len_pos, end - start)
                logging.info(
                    "Estimated remaining time: %.0f s", len_pos / 2 * (end - start)
                )
        return ls_distance_data

    def _distances_refinement(self, tag, chain1, ms_keys):
        logging.info("FRET score refinement %s", chain1)
        ls_distance_data = {}

        coords = self._generate_coords(self.model[tag])

        for ms_key in ms_keys:
            logging.debug("Refine %s", ms_key)
            # print("Refine ", ms_key)
            # start = time.time()
            # len_pos = 0 # for timing
            key, key2 = ms_key.split("_")
            try:
                av = self._calc_or_load_av(
                    tag, key, self.fluo_parameter, load=False, coords=coords
                )
                if self.fluo_parameter != self.fluo_parameter2:
                    av_f2 = self._calc_or_load_av(
                        tag, key, self.fluo_parameter2, load=False, coords=coords
                    )
                else:
                    av_f2 = av

                if len(av.points()[0]) < 1000:
                    logging.debug("%s AV (D) size %d", key, len(av.points()[0]))
                    continue
                if len(av_f2.points()[0]) < 1000:
                    logging.debug("%s AV (A) size %d", key, len(av_f2.points()[0]))
                    continue
            except:
                if self.DEBUG_FRET_SCORE:
                    raise
                av = None
                av_f2 = None
            try:
                av2 = self._calc_or_load_av(
                    tag, key2, self.fluo_parameter, load=False, coords=coords
                )
                if self.fluo_parameter != self.fluo_parameter2:
                    av2_f2 = self._calc_or_load_av(
                        tag, key2, self.fluo_parameter2, load=False, coords=coords
                    )
                else:
                    av2_f2 = av2

                if len(av2.points()[0]) < 1000:
                    logging.debug("%s AV2 (D) size %d", key2, len(av2.points()[0]))
                    continue
                if len(av2_f2.points()[0]) < 1000:
                    logging.debug("%s AV2 (A) size %d", key2, len(av2_f2.points()[0]))
                    continue
            except:
                if self.DEBUG_FRET_SCORE:
                    raise
                av2 = None
                av2_f2 = None

            try:
                n_avg = self.N_KALININ
                # dist = llf.rdaDistance(av, av2_f2, n=n_avg)
                dist = llf.effDistance(av, av2_f2, self.foerster_radius, n=n_avg)
                # dist = llf.rmpDistance(av, av2_f2)

                if self.fluo_parameter != self.fluo_parameter2:
                    # dist2 = llf.rdaDistance(av_f2, av2, n=n_avg)
                    dist2 = llf.effDistance(av_f2, av2, self.foerster_radius, n=n_avg)
                    # dist2 = llf.rmpDistance(av_f2, av2)
                else:
                    dist2 = dist
                ls_distance_data[ms_key] = round((dist + dist2) / 2.0, 1)
            except Exception as ex:
                if self.DEBUG_FRET_SCORE:
                    raise
                logging.debug(type(ex).__name__, key, key2)

            # end = time.time()
            # if self.method == "KALININ":
            # logging.info("Time for %d steps: %.1f s", len_pos, end-start)
            # logging.info("Estimated remaining time: %.0f s", len_pos/2*(end-start))
        return ls_distance_data

    def _calc_heat_map(self, chains_apo, chains_holo):

        pos_apo_1 = {}
        for res in self.model["apo"][chains_apo[0]]:
            res_id = res.get_id()[1]
            c_beta = self._get_cb_position(res)
            pos_apo_1[res_id] = c_beta
        pos_apo_2 = {}
        for res in self.model["apo"][chains_apo[1]]:
            res_id = res.get_id()[1]
            c_beta = self._get_cb_position(res)
            pos_apo_2[res_id] = c_beta

        pos_holo_1 = {}
        for res in self.model["holo"][chains_apo[0]]:
            res_id = res.get_id()[1]
            c_beta = self._get_cb_position(res)
            pos_holo_1[res_id] = c_beta
        pos_holo_2 = {}
        for res in self.model["holo"][chains_apo[1]]:
            res_id = res.get_id()[1]
            c_beta = self._get_cb_position(res)
            pos_holo_2[res_id] = c_beta

        min_int = min(
            min(pos_apo_1.keys()),
            min(pos_apo_2.keys()),
            min(pos_holo_1.keys()),
            min(pos_holo_2.keys()),
        )
        max_int = max(
            max(pos_apo_1.keys()),
            max(pos_apo_2.keys()),
            max(pos_holo_1.keys()),
            max(pos_holo_2.keys()),
        )

        heatmap_array = np.zeros((max_int - min_int, max_int - min_int))
        for i in range(max_int - min_int):
            for j in range(max_int - min_int):
                try:
                    c_beta1 = pos_apo_1[i + min_int]
                    c_beta2 = pos_apo_2[j + min_int]
                    dist_apo = np.linalg.norm(np.array(c_beta1) - np.array(c_beta2))
                    c_beta1 = pos_holo_1[i + min_int]
                    c_beta2 = pos_holo_2[j + min_int]
                    dist_holo = np.linalg.norm(np.array(c_beta1) - np.array(c_beta2))
                    heatmap_array[i][j] = round(dist_holo - dist_apo, 1)
                except Exception:
                    # raise
                    heatmap_array[i][j] = 0  # defined for `self.method == "KALININ"`

        self._json_heatmap(min_int, min_int, heatmap_array)

    def _json_heatmap(self, x0, y0, matrix):
        json_data = {"x0": x0, "y0": y0, "matrix": matrix.tolist()}

        # with open(filepath, 'w') as file:
        #     json.dump(json_data, file, indent=2)
        self.heatmap_data = json_data
