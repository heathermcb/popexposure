"""
popexposure: A package for finding the number of people residing near environmental hazards.
"""

from .pop_estimator import PopEstimator
from .utils.mask_raster_partial_pixel import mask_raster_partial_pixel
from .utils import geom_ops, geom_validator, reader

__all__ = [
    "PopEstimator",
    "geom_ops",
    "geom_validator",
    "mask_raster_partial_pixel",
    "reader",
]

# Package metadata
__version__ = "1.0.0"
__author__ = "heathermcb, joanacasey, nina-flores, lawrence-chillrud, laurenwilner"
