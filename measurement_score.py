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
import os

from . import pdbhelper


class MeasurementScore:
    """
    Abstract class, which provides generic functionality to analyse each
    residue (pair) for the different measurements (e.g. FRET, PIFE, quenching).
    """

    file_tag = None
    measurement_score = {}
    measurement_score_long = {}
    identifier = {}
    labeling_scores = {}
    file_path = None

    structure = {}
    model = {}

    def __init__(self):
        self.measurement_score = {}
        self.measurement_score_long = {}
        self.measurement_score_header = ["ID", "Measurement Score"]
        self.measurement_score_long_header = ["ID", "Measurement Score"]
        self.identifier = {}
        self.labeling_scores = {}
        self.structure = {}
        self.model = {}

    # tag has to be apo or holo
    def load_labeling_score(self, file_path, identifier, tag):
        if tag in ("apo", "holo"):
            self.file_path = file_path
            self.identifier[tag] = identifier
            with open(file_path, "r") as csv_file:
                reader = csv.reader(csv_file)
                self.labeling_scores[tag] = dict(reader)
                if "ID" in self.labeling_scores[tag]:
                    self.labeling_scores[tag].pop("ID")
        else:
            raise ValueError("tag must be 'apo' or 'holo'")

    #    def set_pdb(self,file_path,identifier,structure):
    #        self.file_path = file_path
    #        self.identifier = identifier
    #        self.structure = structure
    #        self.model = self.structure[0]

    # form of file_path
    # file_path = r'G:\Programming\Labelizer\test2.pdb'
    # file_path = 'G:\\Programming\\Labelizer\\test2.pdb'
    def load_pdb(self, file_path, identifier, tag):
        # warnings.simplefilter("ignore")
        if tag in ("apo", "holo"):
            # parser = PDB.PDBParser(PERMISSIVE=1)

            #            self.file_path = file_path
            #            self.identifier = identifier
            # self.structure[tag] = parser.get_structure(identifier,file_path)
            # self.model[tag]=self.structure[tag][0]
            self.structure[tag], self.model[tag] = pdbhelper.load_pdb(
                identifier, file_path
            )
        else:
            raise ValueError("tag must be 'apo' or 'holo'")

    def save_csv(self, tag="short"):
        try:
            identifiers = self.identifier["apo"] + r"_" + self.identifier["holo"]
        except:
            identifiers = self.identifier["apo"]
        if tag == "short":
            data = self.measurement_score
            filepath = os.path.join(
                os.path.dirname(self.file_path),
                identifiers + r"_" + self.file_tag + ".csv",
            )
            header = self.measurement_score_header
        elif tag == "long":
            data = self.measurement_score_long
            filepath = os.path.join(
                os.path.dirname(self.file_path),
                identifiers + r"_" + self.file_tag + tag + ".csv",
            )
            header = self.measurement_score_long_header
        elif tag == "refine":
            data = self.measurement_score_refine
            filepath = os.path.join(
                os.path.dirname(self.file_path),
                identifiers + r"_" + self.file_tag + tag + ".csv",
            )
            header = self.measurement_score_long_header
        else:
            data = None
            filepath = None
            header = None
        if data is not None:
            with open(filepath, "w", newline="") as csv_file:
                writer = csv.writer(csv_file)
                writer.writerow(header)
                for key, value in data.items():
                    line = []
                    line.append(key)
                    for val in value:
                        line.append(val)
                    writer.writerow(line)

    #     def calc_measurement_scores(self,chain1,chain2):
    # #        apo_residues = list(self.model['apo'].get_residues())
    # #        holo_residues = list(self.model['holo'].get_residues())
    #         for apo_key, ls_apo  in self.labeling_scores['apo'].items():
    # #            apo_res = apo_residues[int(apo_key)]
    # #            print("START")
    #             apo_res = self.model['apo'][apo_key[0]][int(apo_key[1:])] #'A'
    # #            print("START LOOP")
    #             for apo_key2, ls_apo2  in self.labeling_scores['apo'].items():
    # #                apo_res2 = apo_residues[int(apo_key2)]
    # #                print("KEY")
    #                 apo_res2 = self.model['apo'][apo_key2[0]][int(apo_key2[1:])]
    # #                if apo_key=='A7' and apo_key2=='B111':
    # #                    print("A7_B111 START")
    # #                print("IF")
    #                 if apo_key!=apo_key2 and (chain1==apo_key[0] and chain2==apo_key2[0]):
    # #                    print("IN IF")
    #                     apo_distance = apo_res['CA'] - apo_res2['CA']
    # #                    if apo_key=='A7' and apo_key2=='B111':
    # #                        print(apo_distance)

    #                     try:
    #                         holo_res = self.model['holo'][apo_key[0]][int(apo_key[1:])] #holo_residues[apo_key]
    #                         ls_holo = self.labeling_scores['holo'][apo_key]
    #                         holo_res2 = self.model['holo'][apo_key2[0]][int(apo_key2[1:])] #holo_residues[apo_key2]
    #                         ls_holo2 = self.labeling_scores['holo'][apo_key2]
    #                         holo_distance = holo_res['CA'] - holo_res2['CA']
    #                         #
    #                         ms = self.calc_measurement_score(float(ls_apo),float(ls_apo2),float(ls_holo),float(ls_holo2),apo_distance,holo_distance)
    #                         if int(apo_key[1:])<int(apo_key2[1:]):
    #                             ms_key = str(apo_key) + '_' + str(apo_key2)
    #                         else:
    #                             ms_key = str(apo_key2) + '_' + str(apo_key)
    #                         self.measurement_score[ms_key] = [ms]
    #                         self.measurement_score_long[ms_key] = [ls_apo,ls_apo2,ls_holo,ls_holo2,apo_distance,holo_distance,ms]
    #                     except:
    #                         pass
    # #                if apo_key=='A7' and apo_key2=='B111':
    # #                    print("A7_B111 END")

    #     def calc_measurement_scores_single(self,chain1,chain2):
    #         logging.info("Calc single measurement score")
    # #        apo_residues = list(self.model['apo'].get_residues())
    # #        holo_residues = list(self.model['holo'].get_residues())
    #         for apo_key, ls_apo  in self.labeling_scores['apo'].items():
    # #            apo_res = apo_residues[int(apo_key)]
    #             apo_res = self.model['apo'][apo_key[0]][int(apo_key[1:])] #'A'
    #             for apo_key2, ls_apo2  in self.labeling_scores['apo'].items():
    # #                apo_res2 = apo_residues[int(apo_key2)]
    #                 apo_res2 = self.model['apo'][apo_key2[0]][int(apo_key2[1:])]
    #                 if apo_key!=apo_key2 and (chain1==apo_key[0] and chain2==apo_key2[0]):
    #                     apo_distance = apo_res['CA'] - apo_res2['CA']

    #                     try:
    # #                        holo_res = self.model['holo'][apo_key[0]][int(apo_key[1:])] #holo_residues[apo_key]
    # #                        ls_holo = self.labeling_scores['holo'][apo_key]
    # #                        holo_res2 = self.model['holo'][apo_key2[0]][int(apo_key2[1:])] #holo_residues[apo_key2]
    # #                        ls_holo2 = self.labeling_scores['holo'][apo_key2]
    # #                        holo_distance = holo_res['CA'] - holo_res2['CA']
    #                         #
    #                         ms = self.calc_measurement_score_single(float(ls_apo),float(ls_apo2),apo_distance)
    #                         if int(apo_key[1:])<int(apo_key2[1:]):
    #                             ms_key = str(apo_key) + '_' + str(apo_key2)
    #                         else:
    #                             ms_key = str(apo_key2) + '_' + str(apo_key)
    #                         self.measurement_score[ms_key] = [ms]
    #                         self.measurement_score_long[ms_key] = [ls_apo,ls_apo2,apo_distance,ms]
    #                     except:
    #                         pass

    def calc_measurement_score(self, label_scores, apo_dist, holo_dist):
        return 0

    def calc_measurement_score_single(self, label_scores, dist):
        return 0
