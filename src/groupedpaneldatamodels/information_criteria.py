from typing import Any, Literal, Type, TYPE_CHECKING
from itertools import product
from copy import deepcopy
from tqdm import trange

import numpy as np

if TYPE_CHECKING:
    from .model import _GroupedPanelModelBase

MAX_PARAM_COMBINATIONS = 1000  # Maximum number of parameter combinations to allow in ic_select


def compute_bic(n, k, var_resid):
    """
    Compute the Bayesian Information Criterion (BIC).

    Parameters:
    n (int): Number of observations.
    k (int): Number of parameters in the model.
    log_likelihood (float): Log-likelihood of the model.

    Returns:
    float: The BIC value.
    """
    return n * np.log(var_resid) + k * np.log(n)


def compute_aic(n, k, var_resid):
    """
    Compute the Akaike Information Criterion (AIC).

    Parameters:
    n (int): Number of observations.
    k (int): Number of parameters in the model.
    log_likelihood (float): Log-likelihood of the model.

    Returns:
    float: The AIC value.
    """
    return n * np.log(var_resid) + 2 * k


def compute_hqic(n, k, var_resid):
    """
    Compute the Hannan-Quinn Information Criterion (HQIC).

    Parameters:
    n (int): Number of observations.
    k (int): Number of parameters in the model.
    log_likelihood (float): Log-likelihood of the model.

    Returns:
    float: The HQIC value.
    """
    return n * np.log(var_resid) + 2 * k * np.log(np.log(n))


def compute_r_squared(y, y_hat):
    return np.corrcoef(y, y_hat)[0, 1] ** 2


# NOTE only resid computed, so R^2 is not available
def compute_statistics(n, k, resid, **kwargs):
    var_resid = np.var(resid, ddof=k)  # Residual variance with degrees of freedom correction
    # var_biased_resid = np.var(resid, ddof=0)  # Residual variance without degrees of freedom correction
    var_biased_resid = np.mean(resid**2)  # Biased residual variance (not corrected for degrees of freedom)

    return {
        "sigma^2": var_resid,
        "AIC": compute_aic(n, k, var_biased_resid) if not kwargs.get("no_aic", False) else None,
        "BIC": compute_bic(n, k, var_biased_resid) if not kwargs.get("no_bic", False) else None,
        "HQIC": compute_hqic(n, k, var_biased_resid) if kwargs.get("include_hqic", False) else None,
        # "R^2": compute_r_squared(y, y_hat) if not kwargs.get("no_r_squared", False) else None,
    }


# NOTE this is still very untested and may not work as expected
# NOTE this should probably be moved to a separate file, but for now it is here
# NOTE should disable bootstrap or any other heavy computations as they are not needed for grid search
def grid_search_by_ic(
    model_cls: Type["_GroupedPanelModelBase"],
    param_ranges: dict[str, list[Any]],
    init_params: dict[str, Any],
    fit_params: dict[str, Any] | None = None,
    ic_criterion: Literal["BIC", "AIC", "HQIC"] = "BIC",
) -> tuple["_GroupedPanelModelBase", dict[str, Any], dict[str, Any]]:
    params = param_ranges.keys()

    # Get all combinations of the parameters
    param_combinations = list(product(*param_ranges.values()))

    if len(param_combinations) > MAX_PARAM_COMBINATIONS:
        raise ValueError("Too many parameter combinations, please reduce the number of parameters or their ranges")

    best_model = None
    best_ic = float("inf")  # Start with a very high IC value
    best_params = None
    results = {}

    for combination in trange(
        len(param_combinations),
        desc=f"Selecting best model for {model_cls.__name__}@{hex(id(model_cls))}",
    ):
        params_dict = dict(zip(params, param_combinations[combination]))
        # Create a copy of the model to avoid modifying the original
        init_params = init_params or {}
        init_params.update(params_dict)
        init_params["use_bootstrap"] = False  # Disable bootstrap for grid search, as it is not needed and can be slow
        model = model_cls(**(init_params))

        # Fit the model
        fit_params = fit_params or {}
        fitted_model = model.fit(**(fit_params))

        # Store the results for this combination
        results[tuple(params_dict.items())] = {
            "IC": fitted_model.IC,
            "params": params_dict,
        }

        # Check if the IC is better than the best found so far
        if fitted_model.IC[ic_criterion] < best_ic:
            best_ic = fitted_model.IC[ic_criterion]
            best_model = fitted_model
            best_params = params_dict

    if best_model is None:
        raise ValueError("No suitable model found based on the given parameter ranges")

    # Return the best model found
    assert best_params is not None, "Best parameters should not be None"
    return best_model, results, best_params
