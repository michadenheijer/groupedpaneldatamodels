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
    compute_statistics as bm_compute_statistics,
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
from tqdm import trange

import pandas as pd
import numpy as np

# Commonly used shared functions (may also put them in utility.py)
# TBD

# Errors
# TBD

# NOTE pass RNG to the models


# Base Class
class _GroupedPanelModelBase:  # type:ignore
    """
    Base class for all models

    Parameters
    ----------
    dependent: array_like
        The dependent variable
    exog: array_like
        The exogenous variables
    weights: array_like
        The weights
    check_rank: bool
        Whether to check the rank of the model, default is True, skipping may improve performance.
    """

    def __init__(
        self,
        dependent: ArrayLike,
        exog: ArrayLike,
        use_bootstrap: bool = False,
        random_state=None,
        **kwargs,
    ):
        # TODO Voor nu alles omzetten naar een array, maar weet niet hoe handig dat altijd is
        # want je verliest wel de namen van de kolommen, misschien net als linearmodels een
        # aparte class hiervoor maken
        # self.dependent = pd.DataFrame(dependent)  # type:ignore
        # self.exog = pd.DataFrame(exog)  # type:ignore
        self.dependent = np.asarray(dependent)
        self.exog = np.asarray(exog)
        self._N, self._T, self._K = self.exog.shape  # type:ignore
        self._constant = False  # Set to False and update when neccesary # FIXME not used
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
    def _has_constant(self) -> bool:
        """
        Returns whether the model has a constant

        Returns
        -------
        bool
            Whether the model has a constant
        """
        return self._constant

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
            raise ValueError("Model has not been fitted yet or no bootstrap was used")
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
        self._fit_datetime = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self._fit_start = process_time()  # Start time for fitting the model

    def _post_fit(self):
        assert self._fit_start is not None, "Fit start time is not set. Did you call _pre_fit()?"
        self._fit_duration = process_time() - self._fit_start  # Calculate the time taken to fit the model

    def fit(self) -> "_GroupedPanelModelBase":
        """
        Fits the model to the data

        Returns
        -------
        self
            The fitted model
        """
        # TODO implement this function
        raise NotImplementedError("Fit function not implemented yet")

    def _get_bootstrap_confidence_intervals(
        self, params: tuple[str], n_boot: int = 50, require_deepcopy=False, **kwargs
    ):
        """
        Computes bootstrap confidence intervals for the parameters

        Parameters
        ----------
        n_boot: int
            The number of bootstrap samples to use

        Returns
        -------
        dict
            The confidence intervals for the parameters
        """

        if not self._use_bootstrap:
            return None

        # FIXME this is possibly the worst method to implement this, but I guess this is it
        estimations = []
        # Create child RNGs so each bootstrap draw has an independent stream.
        seed_seq = SeedSequence(self._random_state)
        child_seqs = seed_seq.spawn(n_boot)
        rngs = [default_rng(s) for s in child_seqs]
        c = deepcopy(self) if require_deepcopy else copy(self)
        c._use_bootstrap = False  # Disable bootstrap for the copied model)

        for i in trange(n_boot, disable=self._hide_progressbar, desc=f"Bootstrap {self._name}@{hex(id(self))}"):
            sample = rngs[i].choice(self.N, replace=True, size=self.N)
            c._rng = rngs[i]  # ensure the copied model inherits its own generator
            c.dependent = self.dependent[sample, :, :]
            c.exog = self.exog[sample, :, :]
            estimations.append(c.fit(**kwargs).params)

        self._bootstrap_estimations = estimations
        self._params_bootstrap_se = {}

        # FIXME standard errors are only correct for beta
        # a solution has to be computed for the other parameters
        for p in params:
            se = np.std([estimation[p] for estimation in estimations], axis=0)
            self._params_bootstrap_se[p] = se

    def get_confidence_intervals(self, confidence_level: float = 0.95) -> dict:
        """
        Returns the confidence intervals for the parameters

        Returns
        -------
        dict
            The confidence intervals for the parameters
        """
        if self._params is None:
            raise ValueError("Model has not been fitted yet")

        if self._params_bootstrap_se is None:
            raise ValueError("Model has not been fitted yet or no bootstrap was used")

        ci = {}
        z = norm.ppf((1 + confidence_level) / 2)  # z-score for the given confidence level
        for param, se in self._params_bootstrap_se.items():
            ci[param] = (self._params[param] - z * se, self._params[param] + z * se)

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

    def to_dict(self) -> dict[str, Any]:
        """
        Converts the model to a dictionary

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
            "bootstrap_estimations": self._bootstrap_estimations if hasattr(self, "_bootstrap_estimations") else None,
            "bootstrap_se": self._params_bootstrap_se if hasattr(self, "_params_bootstrap_se") else None,
            "analytical_se": self._params_analytical_se if hasattr(self, "_params_analytical_se") else None,
            "bootstrap_conf_interval": (
                self.get_confidence_intervals()
                if hasattr(self, "_params_bootstrap_se") or hasattr(self, "_params_analytical_se")
                else None
            ),
            "N": self.N,
            "T": self.T,
            "K": self.K,
        }

    @staticmethod
    def _show_float(value: float, precision: int = 4) -> str:
        try:
            return f"{value:.{precision}f}" if not np.isnan(value) else "N/A"
        except Exception:
            return str(value)

    def summary(
        self,
        confidence_level: float = 0.95,
    ):
        # Ensure the model has been fitted
        if self._params is None:
            raise ValueError("Model has not been fitted yet")

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
            ["AIC", ic.get("AIC", "N/A")],
            ["BIC", ic.get("BIC", "N/A")],
            ["HQIC", ic.get("HQIC", "N/A")],
            ["sigma^2", ic.get("sigma^2", "N/A")],
            ["Seed", self._random_state if self._random_state is not None else "N/A"],
            ["Standard Error type", "Bootstrap" if self._use_bootstrap else "Analytical"],
            ["Confidence Level", confidence_level],
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
        params_values = []
        # PARAMETERS TABLE

        if self._params_bootstrap_se is not None or self._params_analytical_se is not None:
            prev_first_idx = None  # Tracks the first index of the previous parameter entry

            for param in self.params_standard_errors.keys():
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
            if self._params_bootstrap_se is not None and param in self._params_bootstrap_se.keys():
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
    def __init__(
        self,
        dependent: ArrayLike,
        exog: ArrayLike,
        G: int,
        use_bootstrap: bool = False,
        model: Literal["bonhomme_manresa", "su_shi_phillips"] = "bonhomme_manresa",
        heterogeneous_beta: bool = True,
        entity_effects: bool = False,
        **kwargs,
    ):
        super().__init__(dependent, exog, use_bootstrap, **kwargs)

        self._entity_effects = entity_effects

        self._model_type = model
        if self._model_type not in ["bonhomme_manresa", "su_shi_phillips"]:
            raise ValueError("Model must be either 'bonhomme_manresa' or 'su_shi_phillips'")

        self.G = int(G)
        self.heterogeneous_beta = heterogeneous_beta

    def fit(self, **kwargs):
        """
        Fits the model to the data

        Returns
        -------
        self
            The fitted model
        """
        self._pre_fit()
        n_boot = kwargs.pop("n_boot", 50)
        if self._model_type == "bonhomme_manresa":
            # TODO implement this function
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
            # TODO implement correct number of parameters
            num_params = self.G * self.T + self.N + self.K
            self._IC = compute_statistics(self.N * self.T, num_params, resid, include_hqic=True)
            self._get_bootstrap_confidence_intervals(("beta",), n_boot=n_boot, **kwargs)

            self._post_fit()  # Set the fit duration and datetime
            return self
        elif self._model_type == "su_shi_phillips":
            if self.heterogeneous_beta is False:
                raise ValueError("Homogeneous beta is not supported for the Su and Shi Phillips model")

            b, alpha, beta, resid = su_shi_phillips(
                np.squeeze(self.dependent),
                self.exog,
                self.N,
                self.T,
                self.K,
                self.G,
                use_individual_effects=self._entity_effects,
                **kwargs,
            )

            self._params = {"beta": beta, "b": b, "alpha": alpha}
            self._resid = resid  # Store the residuals
            num_params = np.unique_counts(np.round(np.concat([b.ravel(), beta.ravel(), alpha.ravel()]), 3)).counts.sum()
            # num_params = 0
            self._IC = compute_statistics(self.N * self.T, num_params, resid, include_hqic=True)
            self._get_bootstrap_confidence_intervals(("beta",), n_boot=n_boot, **kwargs)

            self._post_fit()  # Set the fit duration and datetime
            return self

        raise ValueError("Model must be either 'bonhomme_manresa' or 'su_shi_phillips'")


class GroupedInteractiveFixedEffects(_GroupedPanelModelBase):
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
        super().__init__(dependent, exog, use_bootstrap, **kwargs)

        self._model_type = model

        if self._model_type not in ["ando_bai", "su_ju"]:
            raise ValueError("Model must be either 'ando_bai' or 'su_ju'")

        self.G = int(G)
        self.GF = (
            GF if GF is not None else np.ones(G, dtype=int)
        )  # NOTE if GF is not defined, we assume all groups have one factor
        self.R = R if R is not None else 1  # Number of factors, default to 1
        self.heterogeneous_beta = heterogeneous_beta

    # FIXME best to change this into multiple functions
    def fit(self, **kwargs):
        """
        Fits the model to the data

        Returns
        -------
        self
            The fitted model
        """
        self._pre_fit()
        n_boot = kwargs.pop("n_boot", 50)
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
            self._params = {"beta": beta.T, "g": g_members, "F": F, "Lambda": Lambda}

            num_params = self.G * self.T + self.GF.sum() + self.N + self.K
            self._resid = resid  # Store the residuals
            self._IC = compute_statistics(self.N * self.T, num_params, resid, include_hqic=True)
            self._get_bootstrap_confidence_intervals(
                ("beta",), n_boot=n_boot, **kwargs  # type:ignore
            )

            self._post_fit()  # Set the fit duration and datetime
            return self

        elif self._model_type == "su_ju":
            if self.heterogeneous_beta == False:
                raise ValueError("Homogeneous beta is not supported for the Su and Ju model")

            b, beta, lambdas, factors, resid = su_ju(
                self.dependent, self.exog, self.N, self.T, self.K, self.G, R=self.R, **kwargs
            )

            self._params = {"b": b, "beta": beta, "lambdas": lambdas, "factors": factors}
            self._resid = resid
            num_params = np.unique_counts(np.concat([b.ravel(), beta.ravel(), lambdas.ravel()])).counts.sum()
            self._IC = compute_statistics(self.N * self.T, num_params, resid, include_hqic=True)
            self._get_bootstrap_confidence_intervals(("beta",), n_boot=n_boot, **kwargs)

            self._post_fit()  # Set the fit duration and datetime
            return self

        raise ValueError("Model must be either 'ando_bai' or 'su_ju'")
