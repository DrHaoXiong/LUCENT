"""Locus-traceable spatial cell–cell communication evidence."""
__version__ = "0.1.0"
from .spatial import SpatialData, prepare_hypothesis, test_hypothesis
from .calibration import CountCalibrationConfig
__all__ = ["SpatialData", "prepare_hypothesis", "test_hypothesis", "CountCalibrationConfig"]
