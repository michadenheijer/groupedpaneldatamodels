# This file has all the main code of each of the models (based on linearmodels)
# The main logic of each model is implemented in their own file, this file provides the classes
# as the main code is functional and not object oriented


# Imports
# First local imports
from .models.ando_bai import (
    grouped_interactive_effects as ando_bai,
    grouped_interactive_effects_hetrogeneous as ando_bai_heterogeneous,
)
from .models.bonhomme_manresa import (
    grouped_fixed_effects as bonhomme_manresa,
    _compute_statistics as bm_compute_statistics,
)
from .models.su_ju import interactive_effects_estimation as su_ju
from .models.su_shi_phillips import fixed_effects_estimation as su_shi_phillips
from .information_criteria import compute_statistics

# Second standard library imports
from typing import Literal, Any
from copy import deepcopy, copy
from datetime import datetime
from time import process_time

# Third party imports
from statsmodels.iolib.summary import Summary
from statsmodels.iolib.table import SimpleTable
from numpy.typing import ArrayLike
from numpy.random import default_rng, SeedSequence
from scipy.stats import norm
from tqdm import tqdm
from joblib import Parallel, delayed

import numpy as np

# Commonly used shared functions (may also put them in utility.py)
# TBD

# Errors
# TBD

# NOTE pass RNG to the models


# Base Class
class _GroupedPanelModelBase:  # type:ignore
    """
    Base class for grouped panel data models.
    This class provides the basic structure and functionality for grouped panel data models.
    It is not meant to be instantiated directly, but rather to be subclassed by specific models.
    """

    def __init__(
        self,
        dependent: ArrayLike,
        exog: ArrayLike,
        use_bootstrap: bool = False,
        random_state=None,
        **kwargs,
    ):
        """
        Initialize the base class for grouped panel data models.

        This class sets up the core structure used by all grouped panel estimators, including
        the dependent variable, exogenous variables, bootstrap settings, and general configuration options.

        Parameters:
            dependent (ArrayLike): Dependent variable array, expected shape (N, T, 1).
            exog (ArrayLike): Exogenous variable array, expected shape (N, T, K).
            use_bootstrap (bool, optional): Whether to compute bootstrap-based standard errors. Defaults to False.
            random_state (int, optional): Seed for the random number generator. Defaults to None.
            **kwargs: Optional keyword arguments for advanced configuration.
                - hide_progressbar (bool): Whether to suppress progress bars during fitting. Defaults to False.
                - disable_analytical_se (bool): Whether to skip analytical standard error calculation. Defaults to False.
        """
        # TODO Voor nu alles omzetten naar een array, maar weet niet hoe handig dat altijd is
        # want je verliest wel de namen van de kolommen, misschien net als linearmodels een
        # aparte class hiervoor maken
        self.dependent = np.asarray(dependent)
        self.exog = np.asarray(exog)
        self._N, self._T, self._K = self.exog.shape  # type:ignore
        # Parallel‑safe random number generator
        self._rng = default_rng(random_state)
        self._random_state = random_state

        # Set up relevant information that needs to be stored
        self._use_bootstrap = use_bootstrap
        # self._original_index = self.dependent.index
        self._name = self.__class__.__name__
        self._fit_datetime = None  # Time when the model was fitted
        self._fit_start = None  # Start time for fitting the model
        self._fit_duration = None  # Duration of the fitting process
        self._model_type = None  # Type of model, can be used for identification
        self._params = None
        self._IC = None
        self._params_analytical_se = None
        self._params_bootstrap_se = None
        self._hide_progressbar = kwargs.pop("hide_progressbar", False)
        self._disable_analytical_se = kwargs.pop("disable_analytical_se", False)
        self._resid = None  # Residuals of the model, to be computed after fitting

        # TODO implement self._not_null (only if neccesary)
        self._validate_data()  # TODO implement this function

        # TODO implement cov_estimators

    def __str__(self) -> str:
        return f"{self._name} ({self._model_type}) \nShape exog: {self.exog.shape}\nShape dependent: {self.dependent.shape}\n"

    def __repr__(self) -> str:
        return self.__str__() + f"\nid: {hex(id(self))}"

    def _validate_data(self) -> None:
        # TODO not that relevant for now
        pass

    @property
    def N(self) -> int:
        """
        Returns the number of observations

        Returns
        -------
        int
            The number of observations
        """
        return self._N

    @property
    def T(self) -> int:
        """
        Returns the number of time periods

        Returns
        -------
        int
            The number of time periods
        """
        return self._T

    @property
    def K(self) -> int:
        """
        Returns the number of exogenous variables

        Returns
        -------
        int
            The number of exogenous variables
        """
        return self._K

    @property
    def params(self) -> dict:
        """
        Returns the parameters of the model

        Returns
        -------
        dict
            The parameters of the model
        """
        if self._params is None:
            raise ValueError("Model has not been fitted yet")
        return self._params

    @property
    def params_bootstrap_standard_errors(self) -> dict:
        """
        Returns the bootstrap standard errors of the parameters

        Returns
        -------
        dict | None
            The bootstrap standard errors of the parameters, or None if not available
        """
        if self._params_bootstrap_se is None:
            raise ValueError("Model has not been fitted yet or no bootstrap was used")
        return self._params_bootstrap_se

    @property
    def params_analytical_standard_errors(self) -> dict:
        """
        Returns the analytical standard errors of the parameters

        Returns
        -------
        dict | None
            The analytical standard errors of the parameters, or None if not available
        """
        if self._params_analytical_se is None:
            raise ValueError("Model has not been fitted yet or no analytical se was implemented")
        return self._params_analytical_se

    @property
    def params_standard_errors(self) -> dict:
        """
        Returns the standard errors of the parameters

        Returns
        -------
        dict | None
            The standard errors of the parameters, or None if not available
        """
        if not self._use_bootstrap:
            return self.params_analytical_standard_errors

        return self.params_bootstrap_standard_errors

    @property
    def t_values(self) -> dict:
        """
        Returns the t-values of the parameters

        Returns
        -------
        dict | None
            The t-values of the parameters, or None if not available
        """
        return {param: self.params[param] / se for param, se in self.params_standard_errors.items()}

    @property
    def resid(self) -> np.ndarray:
        """
        Returns the residuals of the model

        Returns
        -------
        ArrayLike
            The residuals of the model

        Raises
        ------
        ValueError
            If the model has not been fitted yet or residuals are not available for this model.
        """
        if self._resid is None:
            raise ValueError("Model has not been fitted yet or residuals are not available for this model")
        return self._resid

    def p_values(self) -> dict:
        """
        Returns the p-values of the parameters

        Returns
        -------
        dict | None
            The p-values of the parameters, or None if not available
        """
        return {param: 2 * (1 - norm.cdf(np.abs(t))) for param, t in self.t_values.items()}

    @property
    def IC(self) -> dict:
        """
        Returns the information criteria of the model

        Returns
        -------
        dict | None
            The information criteria of the model, or None if not available
        """
        if self._IC is None:
            raise ValueError("Model has not been fitted yet or IC values are not available for this model")
        return self._IC

    # TODO: F-stat, R^2, andere dingen

    # FIXME add more pre-fit checks if needed
    # e.g. check if the data is in the right format, if the dependent variable is a 3D array, etc.
    def _pre_fit(self):
        """
        Internal method to prepare the model for fitting.
        """
        self._fit_datetime = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self._fit_start = process_time()  # Start time for fitting the model

    def _post_fit(self):
        """
        Internal method to finalize the model after fitting.
        This method should be called after the model has been fitted.
        """
        assert self._fit_start is not None, "Fit start time is not set. Did you call _pre_fit()?"
        self._fit_duration = process_time() - self._fit_start  # Calculate the time taken to fit the model

    def fit(self) -> "_GroupedPanelModelBase":
        """
        Fit the grouped panel data model to the provided dataset.

        This method should be implemented by subclasses of `_GroupedPanelModelBase`. It defines
        the main estimation routine and is responsible for setting the model's fitted parameters,
        residuals, and any diagnostics such as information criteria or standard errors.

        Returns:
            _GroupedPanelModelBase: The fitted model instance.

        Raises:
            NotImplementedError: This base method must be overridden by a subclass.
        """
        # TODO implement this function
        raise NotImplementedError("Fit function not implemented yet")

    def _get_bootstrap_standard_errors(
        self, params: tuple[str], n_boot: int = 50, require_deepcopy=False, n_jobs=-1, **kwargs
    ):
        """
        Computes bootstrap standard errors for the parameters.

        Parameters
        ----------
        params: tuple[str], required, displays which parameters to compute the bootstrap standard errors for
        n_boot: int, optional, the number of bootstrap samples to use, default is 50
        require_deepcopy: bool, optional, whether to require a deepcopy of the model, default is False
        **kwargs: Additional keyword arguments that can be passed to the model fitting function.

        Returns
        -------
        dict
            The confidence intervals for the parameters
        """

        if not self._use_bootstrap:
            return None

        # Prepare parallel bootstrap estimations using joblib
        seed_seq = SeedSequence(self._random_state)
        child_seqs = seed_seq.spawn(n_boot)
        rngs = [default_rng(s) for s in child_seqs]

        def _bootstrap_iteration(rng):
            model_copy = deepcopy(self) if require_deepcopy else copy(self)
            model_copy._use_bootstrap = False
            model_copy._disable_analytical_se = True
            sample = rng.choice(self.N, replace=True, size=self.N)
            model_copy._rng = rng
            model_copy.dependent = self.dependent[sample, :, :]
            model_copy.exog = self.exog[sample, :, :]
            return model_copy.fit(**kwargs).params

        estimations = Parallel(n_jobs=n_jobs)(
            delayed(_bootstrap_iteration)(rng) for rng in tqdm(rngs, disable=self._hide_progressbar)
        )

        self._bootstrap_estimations = estimations
        self._params_bootstrap_se = {}

        # FIXME standard errors are only correct for beta
        # a solution has to be computed for the other parameters
        for p in params:
            se = np.std([estimation[p] for estimation in estimations], axis=0, ddof=1)  # type:ignore
            self._params_bootstrap_se[p] = se

    def _get_analytical_standard_errors(self):
        """
        Computes analytical standard errors for the parameters.
        This function is a placeholder and should be implemented in subclasses.

        Returns
        -------
        dict
            The analytical standard errors for the parameters
        """
        # TODO implement this function
        raise NotImplementedError("Analytical standard errors function not implemented yet")

    def get_confidence_intervals(
        self, confidence_level: float = 0.95, conf_type: Literal["auto", "bootstrap", "analytical"] = "auto"
    ) -> dict:
        """
        Returns the confidence intervals for the parameters, prefers bootstrap if available, otherwise analytical

        Parameters
        ----------
        confidence_level: float, optional, the confidence level for the intervals, default is 0.95
        conf_type: str, optional, the type of confidence intervals to compute, can be 'auto', 'bootstrap', or 'analytical', default is 'auto'

        Returns
        -------
        dict
            The confidence intervals for the parameters
        """
        if self._params is None:
            raise ValueError("Model has not been fitted yet")

        ci = {}
        z = norm.ppf((1 + confidence_level) / 2)  # z-score for the given confidence level

        if conf_type == "auto":
            for param, se in self.params_standard_errors.items():
                ci[param] = (self._params[param] - z * se, self._params[param] + z * se)
        elif conf_type == "bootstrap":
            for param, se in self.params_bootstrap_standard_errors.items():
                ci[param] = (self._params[param] - z * se, self._params[param] + z * se)
        elif conf_type == "analytical":
            for param, se in self.params_analytical_standard_errors.items():
                ci[param] = (self._params[param] - z * se, self._params[param] + z * se)
        else:
            raise ValueError("conf_type must be one of 'auto', 'bootstrap', or 'analytical'")

        return ci

    def predict(
        self,
        params: ArrayLike,
        *,
        exog: ArrayLike | None = None,
        data: ArrayLike | None = None,
    ) -> ArrayLike:
        """
        Predicts the dependent variable based on the parameters and exogenous variables
        This function is a placeholder and should be implemented in subclasses.

        Parameters
        ----------
        params: array_like
            The parameters
        exog: array_like
            The exogenous variables

        Returns
        -------
        array_like
            The predicted dependent variable
        """
        if exog is None:
            exog = self.exog

        # TODO implement this function
        raise NotImplementedError("Predict function not implemented yet")

    def to_dict(self, store_bootstrap_iterations=False) -> dict[str, Any]:
        """
        Converts the model to a dictionary,
        which can be useful for serialization or inspection.

        Parameters
        ----------
        store_bootstrap_iterations: bool, optional, whether to store the bootstrap iterations in the dictionary, default is False, requires a lot of memory
        If set to True, the bootstrap estimations will be included in the dictionary.

        Returns
        -------
        dict
            The model as a dictionary
        """
        return {
            "name": self._name,
            "id": hex(id(self)),
            "fit_datetime": self._fit_datetime,
            "fit_duration": self._fit_duration,
            "model_type": self._model_type,
            "params": self._params,
            "IC": self._IC,
            "use_bootstrap": self._use_bootstrap,
            "bootstrap_estimations": (
                self._bootstrap_estimations
                if store_bootstrap_iterations and hasattr(self, "_bootstrap_estimations")
                else None
            ),
            "bootstrap_se": self._params_bootstrap_se if hasattr(self, "_params_bootstrap_se") else None,
            "analytical_se": self._params_analytical_se if hasattr(self, "_params_analytical_se") else None,
            "bootstrap_conf_interval": (
                self.get_confidence_intervals(conf_type="bootstrap") if self._params_bootstrap_se is not None else None
            ),
            "analytical_conf_interval": (
                self.get_confidence_intervals(conf_type="analytical")
                if self._params_analytical_se is not None
                else None
            ),
            "N": self.N,
            "T": self.T,
            "K": self.K,
        }

    @staticmethod
    def _show_float(value: float, precision: int = 4) -> str:
        """Formats a float value to a string with a specified precision.

        Args:
            value (float): The float value to format.
            precision (int, optional): The precision that is required. Defaults to 4.

        Returns:
            str: The formatted string representation of the float value.
        """
        try:
            return f"{value:.{precision}f}" if not np.isnan(value) else "N/A"
        except Exception:
            return str(value)

    # FIXME add flexibility to choose between bootstrap and analytical standard errors
    def summary(
        self,
        confidence_level: float = 0.95,
        standard_errors: Literal["auto", "bootstrap", "analytical"] = "auto",
    ):
        """
        Generates a summary of the model, including information about the model type, observations, exogenous variables,
        groups, fit time, fit duration, and standard errors.

        Args:
            confidence_level (float, optional): The preferred confidence lebel . Defaults to 0.95.
            standard_errors (Literal["auto", "bootstrap", "analytical"], optional): Which type of errors to prefer. Defaults to "auto".
            If "auto", it will use bootstrap if available, otherwise analytical standard errors. All other values will raise a NotImplementedError.
            As this is not implemented yet.

        Raises:
            ValueError: Model has not been fitted yet
            NotImplementedError: Standard errors type other than 'auto' is not implemented yet, you can manually view them using `model.params_bootstrap_standard_errors` or `model.params_analytical_standard_errors`_

        Returns:
            Summary: A summary object containing the model information, parameters, and standard errors.
        """
        # Ensure the model has been fitted
        if self._params is None:
            raise ValueError("Model has not been fitted yet")

        if standard_errors != "auto":
            raise NotImplementedError(
                "Standard errors type other than 'auto' is not implemented yet"
                + " you can manually view them using `model.params_bootstrap_standard_errors` or `model.params_analytical_standard_errors`"
            )

        # INFORMATION HEADERS
        left = [
            ["Type", f"{self._model_type}"],
            ["Observations", self.N * self.T],
            ["Exogenous Variables", self.K],
            ["Groups", self.G if hasattr(self, "G") else "N/A"],  # type:ignore
            ["Fit Time", self._fit_datetime],
            ["Fit Duration", f"{self._fit_duration:.2f} seconds"],
            [
                "Hetrogeneous Beta",
                self.heterogeneous_beta if hasattr(self, "heterogeneous_beta") else "N/A",  # type:ignore
            ],
        ]

        ic = self._IC if self._IC is not None else {}

        right = [
            ["AIC", self._show_float(ic.get("AIC", float("nan")))],
            ["BIC", self._show_float(ic.get("BIC", float("nan")))],
            ["HQIC", self._show_float(ic.get("HQIC", float("nan")))],
            ["sigma^2", self._show_float(ic.get("sigma^2", float("nan")))],
            ["Seed", self._random_state if self._random_state is not None else "N/A"],
            ["Standard Error type", "Bootstrap" if self._use_bootstrap else "Analytical"],
            ["Confidence Level", self._show_float(confidence_level)],
        ]

        top = [left + right for left, right in zip(left, right)]
        headers_top = ["Left", "Value", "Right", "Value"]

        summary = Summary()

        summary.tables.append(SimpleTable(top, title=f"{self._name} Summary"))

        # Coef.	Std.Err.	t	P>|t|	[0.025	0.975]
        headers_params = [
            "Parameter",
            "Coef.",
            "Std.Err.",
            "t",
            "P>|t|",
            f"[{(1 - confidence_level)/2:.3f}",
            f"{(1 - (1 - confidence_level)/2):.3f}]",
        ]
        # PARAMETERS TABLE

        if self._params_bootstrap_se is not None or self._params_analytical_se is not None:
            prev_first_idx = None  # Tracks the first index of the previous parameter entry

            for param in self.params_standard_errors.keys():
                params_values = []
                for idx, v in np.ndenumerate(self.params[param]):
                    se = self.params_standard_errors[param][idx]
                    t_value = self.t_values[param][idx]
                    p_value = self.p_values()[param][idx]
                    ci_lower = self.get_confidence_intervals(confidence_level)[param][0][idx]
                    ci_upper = self.get_confidence_intervals(confidence_level)[param][1][idx]

                    # Insert empty row if first index changes
                    first_idx = idx[0] if len(idx) > 0 else None
                    if prev_first_idx is not None and first_idx != prev_first_idx:
                        params_values.append([""] * len(headers_params))

                    prev_first_idx = first_idx

                    # Format row
                    row = [
                        param + str(idx).replace("(", "").replace(")", "").replace(" ", ""),
                        self._show_float(v),
                        self._show_float(se),
                        self._show_float(t_value),
                        self._show_float(p_value),
                        self._show_float(ci_lower),
                        self._show_float(ci_upper),
                    ]
                    params_values.append(row)

                params_se_table = SimpleTable(
                    params_values,
                    headers=headers_params,
                    title=f"{param}",
                )
                summary.tables.append(params_se_table)

        for param in self.params.keys():
            if (
                (self._params_bootstrap_se is not None) or (self._params_analytical_se is not None)
            ) and param in self.params_standard_errors.keys():
                continue

            if self.params[param] is None:
                continue

            if isinstance(self.params[param], dict):
                table = SimpleTable(
                    [[a, b] for a, b in self.params[param].items()], title=f"{param} coef.", headers=["Index", "Value"]
                )
                summary.tables.append(table)

            elif self.params[param].ndim == 1:
                data = self.params[param].round(4)
                table_data = [["Value"] + data.tolist()]
                headers = ["Index"] + [f"{i}" for i in range(len(data))]
                table = SimpleTable(
                    table_data,
                    title=f"{param} coef.",
                    headers=headers,
                )
                summary.tables.append(table)
            else:
                data = self.params[param].round(4).T
                index_column = [[f"{i}"] for i in range(data.shape[0])]
                table_data = [row + value for row, value in zip(index_column, data.tolist())]
                table = SimpleTable(
                    table_data,
                    title=f"{param} coef.",
                    headers=[f"{param}"] + [f"{i}" for i in range(data.shape[0])],
                )
                # If no standard errors are available, just show the parameter values
                summary.tables.append(table)
                # break  # Only show the first parameter if no standard errors are available

        return summary


class GroupedFixedEffects(_GroupedPanelModelBase):
    """
    GroupedFixedEffects

    Class for estimating grouped fixed effects in panel data models.

    Supports both Bonhomme and Manresa (2015) and Su, Shi, and Phillips (2016) estimators.

    This class clusters units into a specified number of latent groups and estimates either:
    - Group-specific slope coefficients (heterogeneous), or
    - A shared slope coefficient (homogeneous)

    depending on the specified configuration. It also optionally includes individual
    (entity-specific) fixed effects and supports bootstrap inference.

    Typical usage:
        >>> model = GroupedFixedEffects(y, x, G=3)
        >>> model.fit(max_iter=100)
        >>> model.summary()

    Attributes such as coefficients, group assignments, and residuals are made available after fitting.
    """

    def __init__(
        self,
        dependent: ArrayLike,
        exog: ArrayLike,
        G: int,
        use_bootstrap: bool = False,
        model: Literal["bonhomme_manresa", "su_shi_phillips"] = "bonhomme_manresa",
        heterogeneous_beta: bool = True,
        entity_effects: bool = False,
        kappa: float = 0.1,
        **kwargs,
    ):
        """
        Initialize the GroupedFixedEffects model for panel data analysis.

        This class estimates grouped fixed effects using either the Bonhomme and Manresa (2015)
        or Su, Shi, and Phillips (2016) estimators. The model is designed for settings where units
        can be grouped based on similarities in their fixed effects or slope coefficients.

        Parameters:
            dependent (np.ndarray): A 3D array of Y, structured as (N, T, 1), where N is the number of individuals, T is the number of time periods.
            exog (np.ndarray): A 3D array of X, structured as (N, T, K), containing the exogenous regressors.
            G (int): The (maximum) number of latent groups to estimate.
            use_bootstrap (bool, optional): Whether or not to estimate the standard errors using the Bootstrap method. Defaults to False.
            model (Literal["bonhomme_manresa", "su_shi_phillips"], optional): Which estimator to use: "bonhomme_manresa" (default) or "su_shi_phillips".
            heterogeneous_beta (bool, optional): Whether the coefficients β should be allowed to differ across groups. If False, they are homogeneous. Defaults to True.
            entity_effects (bool, optional): Whether to include individual (entity-specific) fixed effects in the estimation. Recommended for "su_shi_phillips". Defaults to False.
            kappa (float, optional): Regularization strength for SCAD penalty. Default is 0.1.

        Raises:
            ValueError: If the specified model is not supported.
        """
        super().__init__(dependent, exog, use_bootstrap, **kwargs)

        self._entity_effects = entity_effects

        self._model_type = model
        if self._model_type not in ["bonhomme_manresa", "su_shi_phillips"]:
            raise ValueError("Model must be either 'bonhomme_manresa' or 'su_shi_phillips'")

        if self._model_type == "bonhomme_manresa":
            if kappa != 0.1:  # NOTE 0.1 is the default value for kappa
                raise ValueError("kappa is not supported for the Bonhomme and Manresa model")

        self._kappa = kappa  # Regularization strength for SCAD penalty, only used in Su and Shi Phillips

        self.G = int(G)
        self.heterogeneous_beta = heterogeneous_beta

    def _fit_bm(self, n_boot: int = 50, **kwargs):
        """
        Fits the Bonhomme and Manresa model to the data

        Parameters
        ----------
        n_boot: int
            The number of bootstrap samples to use

        Returns
        -------
        self
            The fitted model
        """
        boot_n_jobs = kwargs.pop("boot_n_jobs", -1)  # type:ignore
        b, beta, g, eta, iterations, objective_value, resid = bonhomme_manresa(
            self.dependent,
            self.exog,
            self.G,
            hetrogeneous_theta=self.heterogeneous_beta,
            unit_specific_effects=self._entity_effects,
            **kwargs,
        )

        # Create dictionary mapping group number to list of individuals
        g_members = {int(group): np.where(g == group)[0].tolist() for group in np.unique(g)}
        self._params = {"beta": b.T, "alpha": beta, "g": g_members, "eta": eta}
        self._resid = resid  # Store the residuals
        assert resid is not None, "Residuals must be computed before calculating standard errors"
        num_params = self.G * self.T + self.N + self.K
        self._IC = compute_statistics(self.N * self.T, num_params, resid, include_hqic=True)
        self._get_analytical_standard_errors()
        self._get_bootstrap_standard_errors(("beta",), n_boot=n_boot, n_jobs=boot_n_jobs, **kwargs)
        self._post_fit()  # Set the fit duration and datetime
        return self

    def _fit_ssp(self, n_boot: int = 50, **kwargs):
        """
        Fits the Su and Shi Phillips model to the data

        Parameters
        ----------
        n_boot: int
            The number of bootstrap samples to use

        Returns
        -------
        self
            The fitted model
        """
        if self.heterogeneous_beta is False:
            raise ValueError("Homogeneous beta is not supported for the Su and Shi Phillips model")

        boot_n_jobs = kwargs.pop("boot_n_jobs", -1)  # type:ignore
        b, alpha, beta, resid = su_shi_phillips(
            np.squeeze(self.dependent),
            self.exog,
            self.N,
            self.T,
            self.K,
            self.G,
            use_individual_effects=self._entity_effects,
            kappa=self._kappa,
            **kwargs,
        )
        self._params = {"beta": beta.T, "b": b, "alpha": alpha}
        self._resid = resid  # Store the residuals

        # Get groupings
        b = b.T
        beta = beta.T
        dists = np.linalg.norm(b[:, None, :] - beta[None, :, :], axis=2)
        g = np.argmin(dists, axis=1)

        g_members = {int(group): np.where(g == group)[0].tolist() for group in np.unique(g)}
        self._params["g"] = g_members

        num_params = np.unique_counts(np.round(np.concat([b.ravel(), beta.ravel(), alpha.ravel()]), 2)).counts.sum()
        self._IC = compute_statistics(self.N * self.T, num_params, resid, include_hqic=True)
        self._get_analytical_standard_errors()
        self._get_bootstrap_standard_errors(("beta",), n_boot=n_boot, n_jobs=boot_n_jobs, **kwargs)
        self._post_fit()  # Set the fit duration and datetime
        return self

    def fit(self, **kwargs):
        """
        Fit the GroupedFixedEffects model to the data.

        This method estimates grouped fixed effects based on the selected model type:
        - "bonhomme_manresa" implements the algorithm by Bonhomme and Manresa (2015).
        - "su_shi_phillips" implements the algorithm by Su, Shi, and Phillips (2016).

        Parameters
        ----------
        n_boot : int, optional
            Number of bootstrap replications to compute standard errors. Default is 50.

        For model='bonhomme_manresa':
            max_iter : int, optional
                Maximum number of optimization iterations. Default is 10000.
            tol : float, optional
                Convergence tolerance for optimization. Default is 1e-6.
            gfe_iterations : int, optional
                Number of different starting values for the grouped fixed effects algorithm. Default is 100.
            enable_vns : bool, optional
                Whether to use the Variable Neighborhood Search (VNS) algorithm. Default is False.
                (Not recommended when `heterogeneous_beta=True` due to computational cost.)

        For model='su_shi_phillips':
            max_iter : int, optional
                Maximum number of optimization iterations. Default is 1000.
            tol : float, optional
                Convergence tolerance for optimization. Default is 1e-6.
            only_bfgs : bool, optional
                Whether to use only the L-BFGS optimizer. If False, alternates between L-BFGS and Nelder-Mead. Default is True.

        Returns
        -------
        self : GroupedFixedEffects
            The fitted model instance.
        """
        self._pre_fit()
        n_boot = kwargs.pop("n_boot", 50)
        if self._model_type == "bonhomme_manresa":
            return self._fit_bm(n_boot=n_boot, **kwargs)
        elif self._model_type == "su_shi_phillips":
            return self._fit_ssp(n_boot=n_boot, **kwargs)

        raise ValueError("Model must be either 'bonhomme_manresa' or 'su_shi_phillips'")

    # FIXME this could be more generalized I believe
    def _get_analytical_standard_errors_bm(self):
        """
        Computes analytical standard errors for the Bonhomme and Manresa model.
        """
        if self._params is None:
            raise ValueError("Model has not been fitted yet")

        if self._model_type != "bonhomme_manresa":
            raise ValueError("This function is only applicable for the Bonhomme and Manresa model")

        se_alpha = np.zeros((self.G, self.T))  # type:ignore
        se_beta = np.zeros((self.G, self.K))  # type:ignore

        g = self.params["g"] if self.heterogeneous_beta else {0: list(range(self.N))}
        resid = self._resid  # type:ignore
        assert resid is not None, "Residuals must be computed before calculating standard errors"

        for gamma in g.keys():
            se_alpha[gamma] = np.sqrt((resid[g[gamma]] ** 2).mean(axis=0))

        self._params_analytical_se = {"alpha": se_alpha}

        for gamma in g.keys():
            x = self.exog[g[gamma]]
            x_bar = self.exog[g[gamma]].mean(axis=0)
            x_diff = x - x_bar[np.newaxis, :, :]

            total_sum_sigma = np.zeros((self.K, self.K))
            for i in range(x.shape[0]):
                for t in range(self.T):
                    total_sum_sigma += x_diff[i, t, np.newaxis].T @ x_diff[i, t, np.newaxis]

            sigma_beta_g = total_sum_sigma / (x.shape[0] * self.T)

            resid_g = resid[g[gamma]]
            total_sum_omega = np.zeros((self.K, self.K))
            for i in range(x.shape[0]):
                for t in range(self.T):
                    for s in range(self.T):
                        diff_t = x_diff[i, t, np.newaxis]
                        diff_s = x_diff[i, s, np.newaxis]
                        total_sum_omega += resid_g[i, t] * resid_g[i, s] * (diff_t.T @ diff_s)

            omega_beta_g = total_sum_omega / (x.shape[0] * self.T)

            var_beta_g = (
                np.linalg.inv(sigma_beta_g) @ omega_beta_g @ np.linalg.inv(sigma_beta_g) / (x.shape[0] * self.T)
            )
            se_beta[gamma] = np.sqrt(np.diag(var_beta_g))

        self._params_analytical_se["beta"] = se_beta

    def _get_analytical_standard_errors_ssp(self):
        """
        Computes analytical standard errors for the Su and Shi Phillips model.
        """
        if self._params is None:
            raise ValueError("Model has not been fitted yet")

        if self._model_type != "su_shi_phillips":
            raise ValueError("This function is only applicable for the Su and Shi Phillips model")

        resid = self._resid  # type:ignore
        assert resid is not None, "Residuals must be computed before calculating standard errors"

        self._params_analytical_se = {}

        if self._entity_effects:
            se_alpha = np.atleast_2d(np.sqrt((resid**2).mean(axis=1))).T
            self._params_analytical_se["alpha"] = se_alpha

        g = self.params["g"]

        se_beta = np.zeros((self.G, self.K))  # type:ignore

        for gamma in self.params["g"].keys():
            x = self.exog[g[gamma]]
            x_bar = self.exog[g[gamma]].mean(axis=0)
            x_diff = x - x_bar[np.newaxis, :, :]

            total_sum_sigma = np.zeros((self.K, self.K))
            for i in range(x.shape[0]):
                for t in range(self.T):
                    total_sum_sigma += x_diff[i, t, np.newaxis].T @ x_diff[i, t, np.newaxis]
            sigma_beta_g = total_sum_sigma / (x.shape[0] * self.T)
            resid_g = resid[g[gamma]]
            total_sum_omega = np.zeros((self.K, self.K))
            for i in range(x.shape[0]):
                for t in range(self.T):
                    for s in range(self.T):
                        diff_t = x_diff[i, t, np.newaxis]
                        diff_s = x_diff[i, s, np.newaxis]
                        total_sum_omega += resid_g[i, t] * resid_g[i, s] * (diff_t.T @ diff_s)
            omega_beta_g = total_sum_omega / (x.shape[0] * self.T)
            var_beta_g = (
                np.linalg.inv(sigma_beta_g) @ omega_beta_g @ np.linalg.inv(sigma_beta_g) / (x.shape[0] * self.T)
            )
            se_beta[gamma] = np.sqrt(np.diag(var_beta_g))

        self._params_analytical_se["beta"] = se_beta

    def _get_analytical_standard_errors(self):
        """
        Computes analytical standard errors for the parameters.
        Note: this function will always be called, as the analytical standard errors are always computed
        Note: this function will only compute the analytical standard errors, it does not return anything

        Returns
        -------
        None
        """
        if self._disable_analytical_se:
            return
        elif self._model_type == "bonhomme_manresa":
            return self._get_analytical_standard_errors_bm()
        elif self._model_type == "su_shi_phillips":
            return self._get_analytical_standard_errors_ssp()

        raise NotImplementedError("Analytical standard errors function not implemented yet for this model")


class GroupedInteractiveFixedEffects(_GroupedPanelModelBase):
    """
    GroupedInteractiveFixedEffects

    Class for estimating grouped interactive fixed effects in panel data.

    Implements the estimators by Ando and Bai (2016) and Su and Ju (2018).

    This class extends the grouped fixed effects framework by allowing for interactive effects
    (latent factors and loadings) that vary across groups. It is suitable for capturing
    unobserved heterogeneity that follows a low-rank factor structure within each group.

    It supports both heterogeneous and homogeneous slope coefficients (only for Ando & Bai, 2016).

    Example usage:
        >>> model = GroupedInteractiveFixedEffects(y, x, G=3)
        >>> model.fit()
        >>> model.summary()

    After fitting, the model provides access to estimated coefficients, group assignments, latent
    factors, and residual diagnostics.
    """

    def __init__(
        self,
        dependent: ArrayLike,
        exog: ArrayLike,
        G: int,
        use_bootstrap: bool = False,
        model: Literal["ando_bai", "su_ju"] = "ando_bai",
        GF: ArrayLike | None = None,
        R: int | None = None,
        heterogeneous_beta: bool = True,
        **kwargs,
    ):
        """
        Initialize the GroupedInteractiveFixedEffects model for panel data analysis.

        This model estimates grouped interactive fixed effects using either the Ando and Bai (2016)
        or the Su and Ju (2018) estimators. It allows for heterogeneity in slope coefficients across
        groups and supports group-specific factor structures to model unobserved interactive effects.

        Parameters
        ----------
        dependent : ArrayLike
            A 3D array representing the dependent variable with shape (N, T, 1),
            where N is the number of individuals and T is the number of time periods.
        exog : ArrayLike
            A 3D array representing the exogenous variables with shape (N, T, K),
            where K is the number of regressors.
        G : int
            The number of latent groups to estimate.
        use_bootstrap : bool, optional
            Whether to compute bootstrap-based standard errors. Default is False.
        model : Literal["ando_bai", "su_ju"], optional
            The estimator to use. Choose between "ando_bai" (default) or "su_ju".
        GF : ArrayLike, optional
            An array specifying the number of latent factors for each group. If not provided,
            a single factor is assumed for each group.
        R : int, optional
            The total number of latent factors in the Su and Ju model. If not specified, defaults to G.
        heterogeneous_beta : bool, optional
            Whether to allow slope coefficients to vary across groups. Default is True.
        **kwargs : dict
            Additional configuration arguments, such as:
                - hide_progressbar (bool): Whether to suppress progress bars.
                - disable_analytical_se (bool): Whether to skip analytical SE calculation.
        """

        super().__init__(dependent, exog, use_bootstrap, **kwargs)

        self._model_type = model

        if self._model_type not in ["ando_bai", "su_ju"]:
            raise ValueError("Model must be either 'ando_bai' or 'su_ju'")

        self.G = int(G)
        self.GF = (
            GF if GF is not None else np.ones(G, dtype=int)
        )  # NOTE if GF is not defined, we assume all groups have one factor
        self.R = R if R is not None else self.G  # Number of factors, default to G
        self.heterogeneous_beta = heterogeneous_beta

    # FIXME best to change this into multiple functions
    def fit(self, **kwargs):
        """
        Fit the GroupedInteractiveFixedEffects model to the data.

        This method estimates grouped interactive fixed effects using the selected estimator:
        - "ando_bai" (Ando & Bai, 2016): Supports both homogeneous and heterogeneous slope coefficients across groups.
        - "su_ju" (Su & Ju, 2018): Supports only heterogeneous slopes and models a shared factor structure across units.

        It estimates group assignments, coefficients, and latent interactive components, and computes
        information criteria and standard errors.

        Parameters
        ----------
        n_boot : int, optional
            Number of bootstrap replications to compute standard errors. Default is 50.
        boot_n_jobs : int, optional
            Number of parallel jobs for the bootstrap procedure. Default is -1 (use all processors).
        kwargs : dict, optional
            Additional arguments passed to the underlying estimator routines.

        Returns
        -------
        self : GroupedInteractiveFixedEffects
            The fitted model instance.

        Raises
        ------
        ValueError
            If the selected model is not one of "ando_bai" or "su_ju", or if homogeneous beta is requested
            for "su_ju", which is not supported.
        """
        self._pre_fit()
        n_boot = kwargs.pop("n_boot", 50)
        boot_n_jobs = kwargs.pop("boot_n_jobs", -1)  # type:ignore
        assert self.GF is not None, "GF must be defined for the model"
        assert isinstance(self.GF, np.ndarray), "GF must be a numpy array"
        if self._model_type == "ando_bai":
            if self.heterogeneous_beta:
                # Use the heterogeneous version of the Ando and Bai model
                beta, g, F, Lambda, objective_value, resid = ando_bai_heterogeneous(
                    self.dependent, self.exog, self.G, self.GF, **kwargs
                )
            else:
                # Use the standard Ando and Bai model
                beta, g, F, Lambda, objective_value, resid = ando_bai(
                    self.dependent, self.exog, self.G, self.GF, **kwargs
                )

            # Create dictionary mapping group number to list of individuals
            g_members = {int(group): np.where(g == group)[0].tolist() for group in np.unique(g)}
            self._params = {"beta": beta.T, "g": g_members, "F": F.T, "Lambda": Lambda}

            num_params = self.G * self.T + self.GF.sum() + self.N + self.K
            self._resid = resid  # Store the residuals
            self._IC = compute_statistics(self.N * self.T, num_params, resid, include_hqic=True)

            self._get_analytical_standard_errors()
            self._get_bootstrap_standard_errors(
                ("beta",), n_boot=n_boot, n_jobs=boot_n_jobs, **kwargs  # type:ignore
            )

            self._post_fit()  # Set the fit duration and datetime
            return self

        elif self._model_type == "su_ju":
            if self.heterogeneous_beta == False:
                raise ValueError("Homogeneous beta is not supported for the Su and Ju model")

            b, beta, lambdas, factors, resid = su_ju(
                self.dependent, self.exog, self.N, self.T, self.K, self.G, R=self.R, **kwargs
            )

            self._params = {"b": b, "beta": beta.T, "lambdas": lambdas, "factors": factors.T}

            # Get groupings
            _b = b.T
            _beta = beta.T
            dists = np.linalg.norm(_b[:, None, :] - _beta[None, :, :], axis=2)
            g = np.argmin(dists, axis=1)
            g_members = {int(group): np.where(g == group)[0].tolist() for group in np.unique(g)}
            self._params["g"] = g_members

            self._resid = resid
            num_params = np.unique_counts(np.concat([b.ravel(), beta.ravel(), lambdas.ravel()])).counts.sum()
            self._IC = compute_statistics(self.N * self.T, num_params, resid, include_hqic=True)
            self._get_analytical_standard_errors()
            self._get_bootstrap_standard_errors(("beta",), n_boot=n_boot, n_jobs=boot_n_jobs, **kwargs)

            self._post_fit()  # Set the fit duration and datetime
            return self

        raise ValueError("Model must be either 'ando_bai' or 'su_ju'")

    # FIXME this code does not actually compute the analytical standard errors for the Ando and Bai model
    # as it computes the Bonhomme and Manresa model standard errors
    def _get_analytical_standard_errors_ab(self):
        """
        Computes analytical standard errors for the Bonhomme and Manresa model.
        """
        if self._params is None:
            raise ValueError("Model has not been fitted yet")

        if self._model_type != "ando_bai":
            raise ValueError("This function is only applicable for the Bonhomme and Manresa model")

        se_beta = np.zeros((self.G, self.K))  # type:ignore

        g = self.params["g"] if self.heterogeneous_beta else {0: list(range(self.N))}
        resid = self._resid  # type:ignore
        assert resid is not None, "Residuals must be computed before calculating standard errors"
        self._params_analytical_se = {}

        for gamma in g.keys():
            assert isinstance(self.GF, np.ndarray), "Group members must be a list of indices"

            GF = np.array(self.GF, dtype=int)

            x = self.exog[g[gamma]]
            # F = self.params["F"]
            # F_g = F[:, self.GF[:gamma].sum() : self.GF[: gamma + 1].sum()].T
            # M_F = np.identity(self.T) - F_g @ np.linalg.pinv(F_g.T @ F_g) @ F_g.T
            # Lambda = self.params["Lambda"]
            # Lambda_g = Lambda[GF[:gamma].sum() : GF[: gamma + 1].sum(), :].T
            # N_g = len(g[gamma])
            # z_init = M_F @ x
            # z_sub = np.zeros_like(z_init)

            # for i in range(N_g):
            #     for j in range(N_g):
            #         c = Lambda_g[i].T @ np.linalg.pinv(Lambda_g.T @ Lambda_g) @ Lambda_g[j]
            #         z_sub[i] = c * z_init[j]

            # z = z_init - z_sub
            # x_diff = z

            x_bar = self.exog[g[gamma]].mean(axis=0)
            x_diff = x - x_bar[np.newaxis, :, :]

            total_sum_sigma = np.zeros((self.K, self.K))
            for i in range(x.shape[0]):
                for t in range(self.T):
                    total_sum_sigma += x_diff[i, t, np.newaxis].T @ x_diff[i, t, np.newaxis]

            sigma_beta_g = total_sum_sigma / (x.shape[0] * self.T)

            resid_g = resid[g[gamma]]
            total_sum_omega = np.zeros((self.K, self.K))

            for i in range(x.shape[0]):
                for t in range(self.T):
                    for s in range(self.T):
                        diff_t = x_diff[i, t, np.newaxis]
                        diff_s = x_diff[i, s, np.newaxis]
                        total_sum_omega += resid_g[i, t] * resid_g[i, s] * (diff_t.T @ diff_s)

            omega_beta_g = total_sum_omega / (x.shape[0] * self.T)

            var_beta_g = (
                np.linalg.inv(sigma_beta_g) @ omega_beta_g @ np.linalg.inv(sigma_beta_g) / (x.shape[0] * self.T)
            )
            se_beta[gamma] = np.sqrt(np.diag(var_beta_g))

        self._params_analytical_se["beta"] = se_beta

    def _get_analytical_standard_errors_sj(self):
        """
        Computes analytical standard errors for the Su and Shi Phillips model.
        """
        if self._params is None:
            raise ValueError("Model has not been fitted yet")

        if self._model_type != "su_ju":
            raise ValueError("This function is only applicable for the Su and Shi Phillips model")

        resid = self._resid  # type:ignore
        assert resid is not None, "Residuals must be computed before calculating standard errors"

        self._params_analytical_se = {}
        g = self.params["g"]

        se_beta = np.zeros((self.G, self.K))  # type:ignore

        for gamma in self.params["g"].keys():
            x = self.exog[g[gamma]]
            x_bar = self.exog[g[gamma]].mean(axis=0)
            x_diff = x - x_bar[np.newaxis, :, :]

            total_sum_sigma = np.zeros((self.K, self.K))
            for i in range(x.shape[0]):
                for t in range(self.T):
                    total_sum_sigma += x_diff[i, t, np.newaxis].T @ x_diff[i, t, np.newaxis]
            sigma_beta_g = total_sum_sigma / (x.shape[0] * self.T)
            resid_g = resid[g[gamma]]
            total_sum_omega = np.zeros((self.K, self.K))
            for i in range(x.shape[0]):
                for t in range(self.T):
                    for s in range(self.T):
                        diff_t = x_diff[i, t, np.newaxis]
                        diff_s = x_diff[i, s, np.newaxis]
                        total_sum_omega += resid_g[i, t] * resid_g[i, s] * (diff_t.T @ diff_s)
            omega_beta_g = total_sum_omega / (x.shape[0] * self.T)
            var_beta_g = (
                np.linalg.inv(sigma_beta_g) @ omega_beta_g @ np.linalg.inv(sigma_beta_g) / (x.shape[0] * self.T)
            )
            se_beta[gamma] = np.sqrt(np.diag(var_beta_g))

        self._params_analytical_se["beta"] = se_beta

    def _get_analytical_standard_errors(self):
        """
        Computes analytical standard errors for the parameters.
        Note: this function will always be called, as the analytical standard errors are always computed

        Returns
        -------
        dict
            The analytical standard errors for the parameters
        """
        if self._disable_analytical_se:
            return
        elif self._model_type == "ando_bai":
            return self._get_analytical_standard_errors_ab()
        elif self._model_type == "su_ju":
            return self._get_analytical_standard_errors_sj()

        raise NotImplementedError("Analytical standard errors function not implemented yet for this model")
