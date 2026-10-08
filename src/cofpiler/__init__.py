"""Statistical stacking models of layered materials such as covalent organic frameworks."""

from importlib.metadata import version

from cofpiler.boltzmann import boltzmann_probabilities
from cofpiler.intercalate import build_intercalated_model
from cofpiler.models import StackingMode, StackingTable, ureg
from cofpiler.stacking import LayerRecord, build_stacked_model

__version__ = version("cofpiler")

__all__ = [
    "LayerRecord",
    "StackingMode",
    "StackingTable",
    "boltzmann_probabilities",
    "build_intercalated_model",
    "build_stacked_model",
    "ureg",
]
