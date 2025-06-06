#   -------------------------------------------------------------
#   Copyright (c) Micha den Heijer.
#   Licensed under the MIT License. See LICENSE in project root for information.
#   -------------------------------------------------------------
"""
Grouped Panel Data Models for Python. Extends statsmodels and linearmodels with some of the most commonly used
Panel Data Models. Implements the following papers:

- Fixed Effects:
    - Bonhomme and Manresa (2015)
    - Su, Shi and Phillips (2016)
- Interactive Effetcs"
    - Ando and Bai (2016)
    - Su and Ju (2019)

Designed to work well with Numpy and Pandas.
"""
# Import all the relevant classes and functions here
from .model import GroupedFixedEffects, GroupedInteractiveFixedEffects


__version__ = "0.0.1"
__all__ = ["GroupedFixedEffects", "GroupedInteractiveFixedEffects"]  # Name them here
