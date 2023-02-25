# -*- coding: utf-8 -*-
# Copyright (C) 2020, Christian Gebhardt (chrirstian.gebhardt@bio.lmu.de)
#
# This file is part of the labelizer distribution and governed by your
# choice of the "MIT license".
# Please see the LICENSE file that should have been included as part of this
# package.

"""
The fluorophore class represents a fluorphore with all important parameters and
general methods for further calculations, e.g. Förster radius determination.
"""

import csv
import logging
import math
import os

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

logger = logging.getLogger(__name__)


class Fluorophore:
    """
    The fluorophore class represents a fluorphore with all important parameters
    such as spectral properties, structural information,etc.
    It provides general methods for further calculations, e.g. Förster radius
    determination.
    """

    # name = "fluorophore-name"
    # quantum_yield = 1
    # extinction_coeff = 200000
    # abs_spectrum = None
    # em_spectrum = None
    # fps_parameter = {'ll': 23.5, 'lw': 4.5, 'R1': 8.1, 'R2': 4.2, 'R3': 1.5}

    skript_folder = os.path.dirname(os.path.realpath(__file__))
    fluo_folder = os.path.join(skript_folder, "resources", "fluorophores")

    def __init__(self):
        """
        Prepare empty fluorophore.
        """
        self.name = ""
        self.quantum_yield = 0.0
        self.extinction_coeff = 0.0
        self.abs_spectrum = None
        self.em_spectrum = None
        self.fluo_para = {}

    def load_from_drive(self, name):
        """
        Load fluorophore parameter from file.
        Method reads fluorophore parameter from file "<name>.fluo" in folder
        "resources/fluorophores" and saves attributes.
        Method reads fluorophore absorption and emission spectra from file "<name>.csv"
        in folder "resources/fluorophores" and prepares spectra to a uniform shape.

        :param name: Acceptor fluorophore for Förster radius calculation.
        :type name: string
        """
        try:
            para_file = os.path.join(Fluorophore.fluo_folder, name + ".fluo")
            with open(para_file, "r") as csv_file:
                reader = csv.reader(csv_file)
                parameters = dict(reader)
                if "name" in parameters:
                    self.name = parameters["name"]
                else:
                    self.name = name
                if "QY" in parameters:
                    self.quantum_yield = float(parameters["QY"])
                if "EC" in parameters:
                    self.extinction_coeff = float(parameters["EC"])
                if "para_ll" in parameters:
                    self.fluo_para["ll"] = float(parameters["para_ll"])
                else:
                    self.fluo_para["ll"] = 20  # educated guess
                if "para_lw" in parameters:
                    self.fluo_para["lw"] = float(parameters["para_lw"])
                else:
                    self.fluo_para["lw"] = 4.5  # educated guess
                if "para_r1" in parameters:
                    self.fluo_para["r1"] = float(parameters["para_r1"])
                else:
                    self.fluo_para["r1"] = 9.0  # educated guess
                if "para_r2" in parameters:
                    self.fluo_para["r2"] = float(parameters["para_r2"])
                else:
                    self.fluo_para["r2"] = 4.5  # educated guess
                if "para_r3" in parameters:
                    self.fluo_para["r3"] = float(parameters["para_r3"])
                else:
                    self.fluo_para["r3"] = 1.5  # educated guess
                # if all (p in parameters for p in ("para_lw", "para_r1", "para_r2", "para_r3")):
                #     self.fluo_para['lw'] = float(parameters['para_lw'])
                #     self.fluo_para['r1'] = float(parameters['para_r1'])
                #     self.fluo_para['r2'] = float(parameters['para_r2'])
                #     self.fluo_para['r3'] = float(parameters['para_r3'])
                # else:
                #     self.fluo_para['lw'] = 4.5
                #     self.fluo_para['r1'] = 9
                #     self.fluo_para['r2'] = 4.5
                #     self.fluo_para['r3'] = 1.5
        except FileNotFoundError:
            logger.warning(f"Fluorophore not found {name}")
            self.name = "unknown"
            self.quantum_yield = 0.0
            self.extinction_coeff = 0.0
        try:
            spec_file = os.path.join(Fluorophore.fluo_folder, name + ".csv")
            spec_data = pd.read_csv(spec_file).fillna(value=0)
            wl_abs, int_abs = self._interpolate_spectrum(
                np.array(spec_data.iloc[:, 0]), np.array(spec_data.iloc[:, 2])
            )
            self.abs_spectrum = np.array([wl_abs, int_abs])
            wl_em, int_em = self._interpolate_spectrum(
                np.array(spec_data.iloc[:, 0]), np.array(spec_data.iloc[:, 1])
            )
            self.em_spectrum = np.array([wl_em, int_em])
        except FileNotFoundError:
            logger.warning(f"Fluorophore spectrum not found {name}")
            self.abs_spectrum = None
            self.em_spectrum = None

    def calc_foerster_radius(self, acceptor):
        """
        Calculate Förster radius between this fluorophore (as donor) and
        second fluorphore (as acceptor).

        :param acceptor: Acceptor fluorophore for Förster radius calculation.
        :type acceptor: :class:`~labelizer.labelizer.Fluorophore`
        :return: Förster radius.
        :rtype: float
        """
        if self.em_spectrum is not None and acceptor.abs_spectrum is not None:
            assert np.alltrue(
                self.em_spectrum[0] == acceptor.abs_spectrum[0]
            ), "Spectra must be of same shape"
            overlap = self._overlap_integral(
                self.em_spectrum[0],
                self.em_spectrum[1],
                acceptor.abs_spectrum[1],
                acceptor.extinction_coeff,
            )
            kappa_2 = 0.666667
            n_ref = 1.4
            foerster_radius = 0.2108 * (
                self.quantum_yield * kappa_2 / (n_ref**4) * np.sum(overlap)
            ) ** (1 / 6.0)
        else:
            foerster_radius = -1.0
        return round(foerster_radius, 1)

    #    def load_spectrum(self, path, abs_em='abs'):
    #        """
    #        Calculate Förster radius between this fluorophore (as donor) and
    #        second fluorphore (as acceptor).
    #        :param acceptor: Acceptor fluorophore for Förster radius calculation.
    #        :type acceptor: `labelizer.labelizer.Fluorophore`
    #        :return: Förster radius.
    #        :rtype: float
    #        """
    #
    #        if abs_em == 'abs':
    #            self.abs_spectrum = path
    #
    #        elif abs_em == 'em':
    #            self.em_spectrum = path
    #        else:
    #            AttributeError("Fluorophore only has a absorption (abs) " + \
    #                           "and emission (em) spectrum")

    @classmethod
    def get_names(cls):
        names = []
        for entry in os.scandir(cls.fluo_folder):
            names.append(entry.name.split(".")[0])
        names = list(set(names))
        names.sort()
        return names

    @staticmethod
    def _interpolate_spectrum(wavelength, intensity, start=300, end=900):
        wl_min = int(math.ceil(wavelength[0]))
        wl_max = int(math.floor(wavelength[-1]))
        if wl_min <= start:
            x_low = []
        else:
            x_low = np.linspace(start, wl_min, num=(wl_min - start), endpoint=False)
        if wl_max >= end:
            x_high = []
        else:
            x_high = np.linspace(wl_max + 1, end, num=(end - wl_max), endpoint=True)
        x_interpol = np.concatenate((x_low, wavelength, x_high))
        y_interpol = np.concatenate(
            (np.zeros(len(x_low)), intensity, np.zeros(len(x_high)))
        )

        spec_interpol = interp1d(x_interpol, y_interpol, kind="cubic")
        x_new = np.linspace(start, end, num=(end - start + 1), endpoint=True)
        y_new = spec_interpol(x_new)
        return x_new, y_new

    @staticmethod
    def _overlap_integral(wavelength, em_donor, abs_acceptor, ex_coeff):
        sum_em = np.sum(em_donor)
        wl_power_4 = np.power(wavelength, 4)

        return (
            np.multiply(wl_power_4, np.multiply(em_donor, abs_acceptor))
            * ex_coeff
            / sum_em
        )
