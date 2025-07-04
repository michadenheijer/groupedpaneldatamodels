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


def compute_statistics(n: int, k: int, resid: np.ndarray, **kwargs):
    """
    Compute residual variance and information criteria based on model residuals.

    Parameters:
        n (int): Total number of observations.
        k (int): Number of parameters in the model.
        resid (np.ndarray): Residuals from the model, shape (N, T).
        **kwargs: Additional keyword arguments to control inclusion of metrics:
            - no_aic (bool): If True, do not compute AIC.
            - no_bic (bool): If True, do not compute BIC.
            - include_hqic (bool): If True, compute HQIC.

    Returns:
        dict[str, float | None]: A dictionary containing:
            - "sigma^2": unbiased residual variance
            - "AIC": Akaike Information Criterion (if computed)
            - "BIC": Bayesian Information Criterion (if computed)
            - "HQIC": Hannan-Quinn Information Criterion (if computed)
    """
    var_resid = np.var(resid, ddof=k)
    var_biased_resid: float = np.mean(resid**2)  # type:ignore

    res = {
        "sigma^2": var_resid,
        "AIC": compute_aic(n, k, var_biased_resid) if not kwargs.get("no_aic", False) else None,
        "BIC": compute_bic(n, k, var_biased_resid) if not kwargs.get("no_bic", False) else None,
        "HQIC": compute_hqic(n, k, var_biased_resid) if kwargs.get("include_hqic", False) else None,
    }

    return res


def grid_search_by_ic(
    model_cls: Type["_GroupedPanelModelBase"],
    param_ranges: dict[str, list[Any]],
    init_params: dict[str, Any],
    fit_params: dict[str, Any] | None = None,
    ic_criterion: Literal["BIC", "AIC", "HQIC"] = "BIC",
) -> tuple["_GroupedPanelModelBase", dict[str, Any], dict[str, Any]]:
    """
    Perform a grid search over model hyperparameters using information criteria.

    This function fits the given model class over all combinations of specified parameters
    and selects the best model based on the specified information criterion (e.g., BIC, AIC, HQIC).

    Parameters:
        model_cls (Type[_GroupedPanelModelBase]): The model class to be fitted (e.g., GroupedFixedEffects, GroupedInteractiveEffects).
        param_ranges (dict[str, list[Any]]): Dictionary where keys are parameter names and values are lists of candidate values.
        init_params (dict[str, Any]): Initial parameters passed to model instantiation.
        fit_params (dict[str, Any] | None): Optional parameters passed to the `.fit()` method of the model.
        ic_criterion (Literal["BIC", "AIC", "HQIC"]): Criterion used to select the best model. Defaults to "BIC".

    Returns:
        tuple[_GroupedPanelModelBase, dict[str, Any], dict[str, Any]]:
            - The best-fitted model.
            - Dictionary of all parameter combinations with their associated IC values.
            - The parameters corresponding to the best model.
    """
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
        init_params.update({"disable_analytical_se": True})  # Disable analytical SE for grid search
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
