import typing
from importlib.metadata import version
from typing import List, Type

from .conservation_score import ConservationScore
from .cysteine_resemblance import CysteineResemblance
from .fluorophore import Fluorophore
from .labelizer import Labelizer, label_score_model_paper
from .methionin_exclusion import AAExclusion
from .secondary_structure import SecondaryStructure
from .solvent_exposure import SolventExposure

# Avoid cyclic imports
if typing.TYPE_CHECKING:
    from .labeling_parameter import LabelingParameter

__title__ = "labelizer"

__version__ = version(__title__)

__summary__ = "A protein analysis toolbox for fluorescent labelling and FRET assays."
__uri__ = "https://github.com/ChristianGebhardt/labelizer"

__author__ = "Christian Gebhardt"
__email__ = "christian.gebhardt@bio.lmu.de"

__license__ = "MIT"
__copyright__ = "Copyright 2020 {}".format(__author__)

available_parameter: List[Type["LabelingParameter"]] = [
    ConservationScore,
    CysteineResemblance,
    AAExclusion,
    SecondaryStructure,
    SolventExposure,
]

# Assert unique tags
assert len(set(i.file_tag for i in available_parameter)) == len(available_parameter)

__all__ = [
    "AAExclusion",
    "ConservationScore",
    "CysteineResemblance",
    "Fluorophore",
    "Labelizer",
    "SecondaryStructure",
    "SolventExposure",
    "__author__",
    "__copyright__",
    "__email__",
    "__license__",
    "__summary__",
    "__title__",
    "__uri__",
    "__version__",
    "available_parameter",
    "label_score_model_paper",
]
