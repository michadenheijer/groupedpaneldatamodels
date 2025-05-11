# This file has all the main code of each of the models (based on linearmodels)


# Imports
# First local imports

# Second standard library imports

# Third party imports
from numpy.typing import ArrayLike

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

    def __init__(
        self,
        dependent: ArrayLike,
        exog: ArrayLike,
        weights: ArrayLike,
        check_rank: bool = True,
    ):
        # TODO Voor nu alles omzetten naar een array, maar weet niet hoe handig dat altijd is
        # want je verliest wel de namen van de kolommen, misschien net als linearmodels een
        # aparte class hiervoor maken
        self.dependent = pd.DataFrame(dependent)  # type:ignore
        self.exog = pd.DataFrame(exog)  # type:ignore
        self.weights = weights  # TODO implement weights, for now these are just ignored
        self.check_rank = check_rank
        self._constant = False  # Set to False and update when neccesary

        # Set up relevant information that needs to be stored
        self._original_shape = self.dependent.shape
        self._original_index = self.dependent.index
        self._constant = False
        self._name = self.__class__.__name__
        self._check_rank = bool(check_rank)
        # TODO implement self._not_null (only if neccesary)
        self._validate_data()  # TODO implement this function

        # Set weights
        self.weights: ArrayLike = self._adapt_weights(weights)  # TODO create this function

        # TODO implement cov_estimators

    def __str__(self) -> str:
        return f"{self._name} \nNum exog: {self.exog.shape}"

    def __repr__(self) -> str:
        return self.__str__() + f"\nid: {hex(id(self))}"

    def _adapt_weights(self, weights: ArrayLike | None) -> ArrayLike:
        """
        Adapts the weights to the model

        Parameters
        ----------
        weights: array_like
            The weights

        Returns
        -------
        array_like
            The adapted weights
        """
        # TODO this function is for now very simple, and should be improved
        # TODO implement weights at all
        if weights is None:
            return np.ones(self.dependent.shape[0])

        return weights

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

    # TODO: F-stat, R^2, andere dingen

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


class GroupedRandomEffects(_GroupedPanelModelBase):
    # TODO implementation
    # Should be the Lin and Ng 2012 model
    pass


class GroupedFixedEffects(_GroupedPanelModelBase):
    # TODO implementation
    # Should be the Bonhomme and Manresa (2015) implementation
    def __init__(
        self,
        dependent: ArrayLike,
        exog: ArrayLike,
        weights: ArrayLike,
        entity_effects: bool = False,
        homogeneous_coefficients: bool = False,
        drop_absorbed: bool = True,
        check_rank: bool = True,
    ):
        super().__init__(dependent, exog, weights, check_rank)

        self._drop_absorbed = drop_absorbed
        self._entity_effects = entity_effects
        self._homogeneous_coefficients = homogeneous_coefficients

    # TODO: Maybe implement a from formula, but not sure if this is needed

    def fit(self, n_sim: int = 20, enable_vns: bool = True) -> tuple[ArrayLike, ArrayLike, float]:
        """
        Fits the model to the data

        Returns
        -------
        None
            The fitted model
        """

        # TODO implement this function
        best_objective: float = np.inf
        best_labels: ArrayLike
        best_centers: ArrayLike

        # TODO implement the entire code into Cython/Numba
        # thus only use static methods and only use numpy
        for _ in range(max(n_sim, 1)):
            centers = self._generate_random_centers(
                self.dependent.values,  # type:ignore
                self.dependent.shape[0],
                self.dependent.shape[1],
                self.exog.shape[1],
                np.random.default_rng(),
            )
            labels, centers, objective = self._kmeans()

            if enable_vns:
                labels, centers, objective = self._vns(labels, centers, objective)

            if objective < best_objective:
                best_objective = objective
                best_labels = labels  # TODO maybe copy
                best_centers = centers  # TODO maybe copy

        return best_labels, best_centers, best_objective  # TODO should return a class with the results

    @staticmethod
    def _generate_random_centers(y: ArrayLike, n: int, t: int, g: int, rng: np.random.Generator) -> ArrayLike:
        """
        Generates random centers for the k-means algorithm
        Returns
        -------
        array_like
            The random centers
        """
        idx = np.random.choice(n, size=g, replace=False)
        return y[idx, :].copy()  # type:ignore

    @staticmethod
    def _kmeans(
        y: ArrayLike,
        centers: ArrayLike,
        labels: ArrayLike,
        rng: np.random.Generator,
        max_iter: int = 100,
    ) -> tuple[ArrayLike, ArrayLike, float]:
        """
        K-means algorithm

        Returns
        -------
        array_like
            The labels
        array_like
            The centers
        float
            The objective function value
        """
        for _ in range(max_iter):
            # Assign labels
            labels = np.argmin(np.linalg.norm(y[:, np.newaxis] - centers, axis=2), axis=1)

            # Update centers
            new_centers = np.array([y[labels == i].mean(axis=0) for i in range(centers.shape[0])])

            # Check convergence
            if np.all(centers == new_centers):
                break

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

    @property
    def entity_effects(self) -> bool:
        """
        Returns whether the model has entity effects

        Returns
        -------
        bool
            Whether the model has entity effects
        """
        return self._entity_effects

    @property
    def homogeneous_coefficients(self) -> bool:
        """
        Returns whether the model has homogeneous coefficients

        Returns
        -------
        bool
            Whether the model has homogeneous coefficients
        """
        return self._homogeneous_coefficients


class GroupedInteractiveEffects(_GroupedPanelModelBase):
    # TODO implementation
    # Should be the Ando and Bai (2016) model
    pass
