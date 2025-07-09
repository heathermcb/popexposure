"""
popexposure: A package for finding the number of people residing near environmental hazards.
"""

# Import the main class that users interact with
from .estimate_exposure import PopEstimator

# Optionally expose other classes/functions users might need
# from .other_module import OtherClass

# Define what gets imported with "from popexposure import *"
__all__ = ["PopEstimator"]

# Package metadata
__version__ = "1.0.0"
__author__ = "heathermcb, joanacasey, nina-flores, lawrence-chillrud, laurenwilner"
