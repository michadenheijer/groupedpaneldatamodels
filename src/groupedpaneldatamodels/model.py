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


# Second standard library imports
from typing import Literal, Any
from copy import deepcopy
from datetime import datetime
from time import process_time

# Third party imports
from numpy.typing import ArrayLike
from scipy.stats import norm
from tqdm import trange

import pandas as pd
import numpy as np

# Commonly used shared functions (may also put them in utility.py)
# TBD

# Errors
# TBD


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

    def __init__(self, dependent: ArrayLike, exog: ArrayLike, use_bootstrap: bool = False):
        # TODO Voor nu alles omzetten naar een array, maar weet niet hoe handig dat altijd is
        # want je verliest wel de namen van de kolommen, misschien net als linearmodels een
        # aparte class hiervoor maken
        # self.dependent = pd.DataFrame(dependent)  # type:ignore
        # self.exog = pd.DataFrame(exog)  # type:ignore
        self.dependent = np.asarray(dependent)
        self.exog = np.asarray(exog)
        self._N, self._T, self._K = self.exog.shape  # type:ignore
        self._constant = False  # Set to False and update when neccesary # FIXME not used

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
        self._params_standard_errors = None

        # TODO implement self._not_null (only if neccesary)
        self._validate_data()  # TODO implement this function

        # TODO implement cov_estimators

    def __str__(self) -> str:
        return f"{self._name} \nShape exog: {self.exog.shape}\nShape dependent: {self.dependent.shape}\n"

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
    def params_standard_errors(self) -> dict:
        """
        Returns the standard errors of the parameters

        Returns
        -------
        dict | None
            The standard errors of the parameters, or None if not available
        """
        if self._params_standard_errors is None:
            raise ValueError("Model has not been fitted yet or no bootstrap was used")
        return self._params_standard_errors

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

    def fit(self):
        """
        Fits the model to the data

        Returns
        -------
        self
            The fitted model
        """
        # TODO implement this function
        raise NotImplementedError("Fit function not implemented yet")

    def _get_bootstrap_confidence_intervals(self, params: tuple[str], n_boot: int = 50, **kwargs):
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

        for i in trange(n_boot):
            c = deepcopy(self)
            c._use_bootstrap = False  # Disable bootstrap for the copied model
            sample = np.random.choice(self.N, replace=True, size=self.N)
            c.dependent = c.dependent[sample, :, :]
            c.exog = c.exog[sample, :, :]
            estimations.append(c.fit(**kwargs).params)

        self._bootstrap_estimations = estimations
        self._params_standard_errors = {}

        # FIXME standard errors are only correct for beta
        # a solution has to be computed for the other parameters
        for p in params:
            se = np.std([estimation[p] for estimation in estimations], axis=0)
            self._params_standard_errors[p] = se

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

        if self._params_standard_errors is None:
            raise ValueError("Model has not been fitted yet or no bootstrap was used")

        ci = {}
        z = norm.ppf((1 + confidence_level) / 2)  # z-score for the given confidence level
        for param, se in self._params_standard_errors.items():
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
            "params_standard_errors": self._params_standard_errors,
            "use_bootstrap": self._use_bootstrap,
            "bootstrap_estimations": self._bootstrap_estimations if hasattr(self, "_bootstrap_estimations") else None,
            "N": self.N,
            "T": self.T,
            "K": self.K,
        }


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
    ):
        super().__init__(dependent, exog, use_bootstrap)

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
            beta, alpha, g, eta, iterations, objective_value = bonhomme_manresa(
                self.dependent,
                self.exog,
                self.G,
                hetrogeneous_theta=self.heterogeneous_beta,
                unit_specific_effects=self._entity_effects,
                **kwargs,
            )
            self._params = {"beta": beta, "alpha": alpha, "g": g, "eta": eta}
            sigma_squared, BIC = bm_compute_statistics(objective_value, self.N, self.T, self.K, self.G)
            self._params.update({"sigma^2": sigma_squared})
            self._IC = {"BIC": BIC}
            self._get_bootstrap_confidence_intervals(("beta",), n_boot=n_boot, **kwargs)

            self._post_fit()  # Set the fit duration and datetime
            return self
        elif self._model_type == "su_shi_phillips":
            if self.heterogeneous_beta is False:
                raise ValueError("Homogeneous beta is not supported for the Su and Shi Phillips model")

            beta, mu, alpha = su_shi_phillips(
                np.squeeze(self.dependent),
                self.exog,
                self.N,
                self.T,
                self.K,
                self.G,
                use_individual_effects=self._entity_effects,
                **kwargs,
            )

            self._params = {"beta": beta, "mu": mu, "alpha": alpha}
            self._IC = None  # TODO implement this
            self._get_bootstrap_confidence_intervals(("alpha",), n_boot=n_boot, **kwargs)

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
    ):
        super().__init__(dependent, exog, use_bootstrap)

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
        if self._model_type == "ando_bai":
            if self.heterogeneous_beta:
                # Use the heterogeneous version of the Ando and Bai model
                beta, g, F, Lambda, objective_value = ando_bai_heterogeneous(
                    self.dependent, self.exog, self.G, self.GF, **kwargs
                )
            else:
                # Use the standard Ando and Bai model
                beta, g, F, Lambda, objective_value = ando_bai(self.dependent, self.exog, self.G, self.GF, **kwargs)

            self._params = {"beta": beta, "g": g, "F": F, "Lambda": Lambda}

            sigma_squared, BIC = bm_compute_statistics(
                objective_value, self.N, self.T, self.K, self.G
            )  # TODO this needs to be updated
            self._params.update({"sigma^2": sigma_squared})  # TODO this as well
            self._IC = {"BIC": BIC}
            self._get_bootstrap_confidence_intervals(
                ("beta",), n_boot=n_boot, **kwargs  # type:ignore
            )

            self._post_fit()  # Set the fit duration and datetime
            return self

        elif self._model_type == "su_ju":
            if self.heterogeneous_beta == False:
                raise ValueError("Homogeneous beta is not supported for the Su and Ju model")

            beta, alpha, lambdas, factors = su_ju(
                self.dependent, self.exog, self.N, self.T, self.K, self.G, R=self.R, **kwargs
            )
            self._params = {"beta": beta, "alpha": alpha, "lambdas": lambdas, "factors": factors}
            self._IC = None  # TODO implement this
            self._get_bootstrap_confidence_intervals(("alpha",), n_boot=n_boot, **kwargs)

            self._post_fit()  # Set the fit duration and datetime
            return self

        raise ValueError("Model must be either 'ando_bai' or 'su_ju'")
