# -*- coding: utf-8 -*-
# Copyright (C) 2020, Christian Gebhardt (chrirstian.gebhardt@bio.lmu.de)
#
# This file is part of the labelizer distribution and governed by your
# choice of the "MIT license".
# Please see the LICENSE file that should have been included as part of this
# package.

"""
The module provides global configuration options.
"""
import configparser
import os


# https://stackoverflow.com/questions/128573/using-property-on-classmethods
class classproperty(object):
    def __init__(self, fget):
        self.fget = fget

    def __get__(self, owner_self, owner_cls):
        return self.fget(owner_cls)


class DefaultConfig(object):
    # Measurement
    LS_THRESHOLD = 0.5
    # FRET
    AV_CALC_RADIUS = 20.0
    COLLISION_RADIUS = 1.7
    SPEED = "precise"
    N_KALININ_FAST = 10000
    GS_FAST = 1.2
    N_KALININ_PRECISE = 100000
    GS_PRECISE = 0.8

    # Pixel sizes are in nm, displayed units are um
    IMG_PIXELSIZE = 80

    # Optimization bounds defaults
    ENDCAP_RANGE = 20.0  # Endcap (`xl`, `xr`) are allowed to varied this many pixels from the initital guesses.

    # Plotting default values
    R_DIST_STOP = 20.0
    R_DIST_STEP = 0.5
    R_DIST_SIGMA = 0.3
    R_DIST_NORM_STOP = 2
    R_DIST_NORM_STEP = 0.05
    R_DIST_NORM_SIGMA = 0.05

    L_DIST_NBINS = 100
    L_DIST_SIGMA = 0.5

    PHI_DIST_STEP = 0.5
    PHI_DIST_SIGMA = 5

    DEBUG = True  # If True, numpy division warnings will be printed.

    # Other
    @classproperty
    def CACHE_DIR(self):
        return os.path.join(os.path.expanduser("~"), ".labelizer", "cache")


config_sections = {
    "General": ["LS_THRESHOLD", "IMG_PIXELSIZE"],
    "Fret": [
        "AV_CALC_RADIUS",
        "COLLISION_RADIUS",
        "SPEED",
        "N_KALININ_FAST",
        "GS_FAST",
        "N_KALININ_PRECISE",
        "GS_PRECISE",
    ],
    "Optimization": ["ENDCAP_RANGE"],
    "Other": ["CACHE_DIR", "DEBUG"],
}

reverse_sections = {vi: k for k, v in config_sections.items() for vi in v}


class ParsedConfig(object):
    getters = {bool: "getboolean", float: "getfloat", int: "getint", str: "get"}

    def __init__(self, config):
        self.config = config

    def __getattr__(self, name):
        _type = type(getattr(DefaultConfig, name))
        get = self.getters[_type]
        return getattr(self.config[reverse_sections[name]], get)(name)


def create_config(path=None):
    """
    Create a config file at the specified `path` using default values. If not
    path is given the file is created in the user's home directory / .labelizer.

    :param path: Path where to create the config file.
    :type path: str, optional
    """

    if not path:
        home = os.path.expanduser("~")
        path = os.path.join(home, ".labelizer")
        if not os.path.isdir(path):
            os.mkdir(path)

    config = configparser.ConfigParser()
    config.optionxform = str
    for k, v in config_sections.items():
        config[k] = {vi: getattr(DefaultConfig, vi) for vi in v}

    with open(os.path.join(path, "config.ini"), "w") as configfile:
        config.write(configfile)

    if not os.path.isdir(config["Other"]["CACHE_DIR"]):
        os.mkdir(config["Other"]["CACHE_DIR"])


def load_config(path=None):
    """
    Load the configuration file at `path`. If not file is given the default
    directory is used. If no config file is present, the default config is used.
    Parameters

    :param path: Path where to create the config file.
    :type path: str, optional
    """
    path = (
        path or os.path.join(os.path.expanduser("~"), ".labelizer", "config.ini")
    )
    global cfg
    if os.path.exists(path):
        config = configparser.ConfigParser()
        config.read(path)
        cfg = ParsedConfig(config)


cfg = DefaultConfig
load_config()

try:
    if not cfg.DEBUG:
        pass
        # np.seterr(divide='ignore', invalid='ignore')
except KeyError:
    print("Invalid config file")
