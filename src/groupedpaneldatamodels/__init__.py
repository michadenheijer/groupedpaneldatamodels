#   -------------------------------------------------------------
#   Copyright (c) Micha den Heijer.
#   Licensed under the MIT License. See LICENSE in project root for information.
#   -------------------------------------------------------------

"""
Grouped Panel Data Models for Python.

This package implements grouped fixed effects and grouped interactive fixed effects estimators for panel data analysis. It is designed to facilitate fast, flexible, and interpretable estimation in settings where unobserved heterogeneity across units can be captured by latent group structures.

Key Features:
    - Estimation of grouped fixed effects (GFE) and grouped interactive fixed effects (GIFE) models.
    - Automatic selection of the number of groups using information criteria such as BIC, AIC, and HQIC.
    - Efficient numerical implementation using vectorized operations and optional bootstrapping.
    - Clean and minimal API for integration into research pipelines.

Modules:
    - `model`: Core estimators `GroupedFixedEffects` and `GroupedInteractiveFixedEffects`.
    - `information_criteria`: Grid search utilities like `grid_search_by_ic` for model selection.

Example:
    >>> from groupedpaneldatamodels import GroupedFixedEffects
    >>> model = GroupedFixedEffects(y, x, G=3)
    >>> results = model.fit()

Author: Micha den Heijer
License: MIT
"""

# Import all the relevant classes and functions here
from .model import GroupedFixedEffects, GroupedInteractiveFixedEffects
from .information_criteria import grid_search_by_ic


__version__ = "0.1.2"  # Version of the package
__all__ = ["GroupedFixedEffects", "GroupedInteractiveFixedEffects", "grid_search_by_ic"]  # Name them here
