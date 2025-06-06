# TODO
# - Check which imports are actually needed
# - Try to remove np.squeeze arguments
# - Getting intitial clusters is quite bad, may need to be improved
# - Cache np.arange in _get_clusters
# - Lambda returns all possible factors and not only the ones that are used, which may be confusing
# - Add docstrings

from numba import njit
from numpy.linalg import lstsq
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from skglm import GeneralizedLinearEstimator
from skglm.penalties import SCAD

import numpy as np

# NOTE suppress FutureWarnings from the SCAD penalty in skglm
import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)


###### Homogeneous Case ######
def _get_factors_initial(y, GF, T):
    y_squeezed = np.squeeze(y).T
    pca = PCA(n_components=GF.sum())
    U = pca.fit_transform(y_squeezed) / pca.singular_values_
    F = np.sqrt(T) * U
    Lambda = F.T @ y_squeezed / T

    return F, Lambda


def _get_factors(y, x, beta, g, G, GF, T):
    y_squeezed = np.squeeze(y - x @ beta).T

    F = np.zeros((T, GF.sum()))
    Lambda = np.zeros((GF.sum(), y_squeezed.shape[1]))

    for i in range(G):
        y_squeezed_partial = np.atleast_2d(np.squeeze(y[g == i] - x[g == i] @ beta).T)

        pca = PCA(n_components=GF[i])
        U = pca.fit_transform(y_squeezed_partial) / pca.singular_values_

        F_partial = np.sqrt(T) * U
        F[:, GF[:i].sum() : GF[: i + 1].sum()] = F_partial

        Lambda_partial = F_partial.T @ y_squeezed / T
        Lambda[GF[:i].sum() : GF[: i + 1].sum(), :] = Lambda_partial

    return F, Lambda


def _get_clusters_initial(Lambda, G, GF):
    # FIXME this code is very bad, but it works
    while True:
        km = KMeans(n_clusters=G)
        g = km.fit_predict(Lambda.T)

        counts = np.bincount(g, minlength=G)
        if np.all(counts >= (GF + 1)):
            break

    return g


def _get_clusters(y, x, beta, Lambda, g, G, F, GF, N, T):
    y_star = y - x @ beta
    res = np.zeros((N, T, G))
    for i in range(G):
        # TODO check if this is correct
        # But I think so
        res[:, :, i] = (
            y_star.reshape(N, -1)
            - (F[:, GF[:i].sum() : GF[: i + 1].sum()] @ Lambda[GF[:i].sum() : GF[: i + 1].sum(), :]).T
        )
    res_per_grouping = (res**2).sum(axis=1)
    g = res_per_grouping.argmin(axis=1)

    # Count size of each group
    counts = np.bincount(g, minlength=G)

    # Ensure minimum group sizes (GF[i] + 1)
    min_sizes = GF + 1

    # Check if any group is below minimum size
    while np.any(counts < min_sizes):
        # Find the most deficient group
        target_group = np.argmin(counts - min_sizes)
        needed = min_sizes[target_group] - counts[target_group]

        if needed <= 0:
            continue  # This group already meets its minimum

        # Find elements not in this group
        non_target_indices = np.where(g != target_group)[0]

        # Sort by distance to target group
        distances = res_per_grouping[non_target_indices, target_group]
        closest_indices = non_target_indices[np.argsort(distances)]

        # Try to reassign closest elements
        reassigned = 0
        for idx in closest_indices:
            source_group = g[idx]

            # Only reassign if source group has enough elements to spare
            if counts[source_group] > min_sizes[source_group]:
                g[idx] = target_group
                counts[source_group] -= 1
                counts[target_group] += 1
                reassigned += 1

                if reassigned >= needed:
                    break

        # If we couldn't reassign any elements, we're stuck
        if reassigned == 0:
            raise Exception("Cannot satisfy minimum group size constraints.")

    objective_value = res_per_grouping[np.arange(N), g].sum()

    return g, objective_value


# NOTE this ignores the factor structure
def _estimate_beta_initial(y, x, K):
    beta = lstsq(x.reshape(-1, K), y.reshape(-1, 1), rcond=None)[0]
    return beta


def _estimate_beta(y, x, g, GF, F, Lambda, K, N, T, G, kappa, gamma):
    res = np.zeros((N, T, G))
    for i in range(G):
        res[:, :, i] = (
            y.reshape(N, -1) - (F[:, GF[:i].sum() : GF[: i + 1].sum()] @ Lambda[GF[:i].sum() : GF[: i + 1].sum(), :]).T
        )

    y_star = res[np.arange(N), :, g]
    # FIXME the np.arange could be cached or something
    if kappa == 0:
        beta = lstsq(x.reshape(-1, K), y_star.reshape(-1, 1), rcond=None)[0]
        return beta

    beta = np.atleast_2d(
        GeneralizedLinearEstimator(penalty=SCAD(alpha=kappa, gamma=gamma))
        .fit(x.reshape(-1, K), np.squeeze(y_star.reshape(-1, 1)))
        .coef_
    ).T
    return beta


def _grouped_interactive_effects_iteration(y, x, G, GF, N, T, K, kappa, gamma, tol, max_iterations):
    last_objective_value = np.inf
    F, Lambda = _get_factors_initial(y, GF, T)
    g = _get_clusters_initial(Lambda, G, GF)
    beta = _estimate_beta_initial(y, x, K)
    F, Lambda = _get_factors(y, x, beta, g, G, GF, T)

    obj_val_store = 5
    objective_values = np.zeros(obj_val_store)

    for i in range(max_iterations):
        g, objective_value = _get_clusters(y, x, beta, Lambda, g, G, F, GF, N, T)
        F, Lambda = _get_factors(y, x, beta, g, G, GF, T)
        beta = _estimate_beta(y, x, g, GF, F, Lambda, K, N, T, G, kappa, gamma)

        objective_values[i % obj_val_store] = objective_value
        if objective_values.max() - objective_values.min() < tol:
            break

        last_objective_value = objective_value

    return beta, g, F, Lambda, last_objective_value


def grouped_interactive_effects(
    y, x, G, GF=None, kappa=0.0, gamma=3.7, tol=1e-6, gife_iterations=100, max_iterations=1000
):
    """Runs GIFE regression"""
    N, T, K = x.shape
    if GF is None:
        GF = np.array([1] * G)
    else:
        GF = np.array(GF)  # Ensures that is an np array

    best_objective_value = np.inf
    best_g = None
    best_beta = None
    best_F = None
    best_Lambda = None

    # FIXME Lamda returns all possible factors and not only the ones that are used

    for i in range(gife_iterations):
        beta, g, F, Lambda, objective_value = _grouped_interactive_effects_iteration(
            y, x, G, GF, N, T, K, kappa, gamma, tol, max_iterations
        )

        if objective_value < best_objective_value:
            best_objective_value = objective_value
            best_g = g
            best_beta = beta
            best_F = F
            best_Lambda = Lambda

    return best_beta, best_g, best_F, best_Lambda, best_objective_value


##### Heterogeneous Case ######
# NOTE this is a copy of the homogeneous case, but with some modifications
def _get_factors_hetrogeneous(y, x, beta, g, G, GF, T, N):
    # NOTE this is not neeeded I believe
    F = np.zeros((T, GF.sum()))
    Lambda = np.zeros((GF.sum(), y.shape[0]))

    for i in range(G):
        y_squeezed_partial = np.atleast_2d(np.squeeze(y[g == i] - x[g == i] @ beta[:, i : i + 1]).T)

        pca = PCA(n_components=GF[i])
        U = pca.fit_transform(y_squeezed_partial) / pca.singular_values_

        F_partial = np.sqrt(T) * U
        F[:, GF[:i].sum() : GF[: i + 1].sum()] = F_partial

        # NOTE that here the correct beta is used
        y_squeezed = np.squeeze(y - x @ beta[:, i : i + 1]).T
        Lambda_partial = F_partial.T @ y_squeezed / T
        Lambda[GF[:i].sum() : GF[: i + 1].sum(), :] = Lambda_partial

    return F, Lambda


# FIXME this code should be modified s.t. cannot return groups less than size
# GF[i]
def _get_clusters_hetrogeneous(y, x, beta, Lambda, g, G, F, GF, N, T):
    res = np.zeros((N, T, G))
    for i in range(G):
        # TODO check if this is correct
        # But I think so
        res[:, :, i] = (
            y.reshape(N, -1)
            - x @ beta[:, i]
            - (F[:, GF[:i].sum() : GF[: i + 1].sum()] @ Lambda[GF[:i].sum() : GF[: i + 1].sum(), :]).T
        )
    res_per_grouping = (res**2).sum(axis=1)
    g = res_per_grouping.argmin(axis=1)

    # Count size of each group
    counts = np.bincount(g, minlength=G)

    # Ensure minimum group sizes (GF[i] + 1)
    min_sizes = GF + 1

    # Check if any group is below minimum size
    while np.any(counts < min_sizes):
        # Find the most deficient group
        target_group = np.argmin(counts - min_sizes)
        needed = min_sizes[target_group] - counts[target_group]

        if needed <= 0:
            continue  # This group already meets its minimum

        # Find elements not in this group
        non_target_indices = np.where(g != target_group)[0]

        # Sort by distance to target group
        distances = res_per_grouping[non_target_indices, target_group]
        closest_indices = non_target_indices[np.argsort(distances)]

        # Try to reassign closest elements
        reassigned = 0
        for idx in closest_indices:
            source_group = g[idx]

            # Only reassign if source group has enough elements to spare
            if counts[source_group] > min_sizes[source_group]:
                g[idx] = target_group
                counts[source_group] -= 1
                counts[target_group] += 1
                reassigned += 1

                if reassigned >= needed:
                    break

        # If we couldn't reassign any elements, we're stuck
        if reassigned == 0:
            raise Exception("Cannot satisfy minimum group size constraints.")

    objective_value = res_per_grouping[np.arange(N), g].sum()

    return g, objective_value


# NOTE this ignores the factor structure
def _estimate_beta_initial_hetrogeneous(y, x, g, K, G):
    beta = np.zeros((K, G))
    for i in range(G):
        y_partial = y[g == i]
        x_partial = x[g == i]

        beta[:, i] = np.squeeze(lstsq(x_partial.reshape(-1, K), y_partial.reshape(-1, 1), rcond=None)[0])
    return beta


def _estimate_beta_hetrogeneous(y, x, g, GF, F, Lambda, K, N, T, G, kappa, gamma):
    beta = np.zeros((K, G))
    res = np.zeros((N, T, G))
    for i in range(G):
        res[:, :, i] = (
            y.reshape(N, -1) - (F[:, GF[:i].sum() : GF[: i + 1].sum()] @ Lambda[GF[:i].sum() : GF[: i + 1].sum(), :]).T
        )

    y_star = res[np.arange(N), :, g]
    if kappa == 0:
        for i in range(G):
            beta[:, i] = np.squeeze(lstsq(x[g == i].reshape(-1, K), y_star[g == i].reshape(-1, 1), rcond=None)[0])
        return beta

    for i in range(G):
        beta[:, i] = (
            GeneralizedLinearEstimator(penalty=SCAD(alpha=kappa, gamma=gamma))
            .fit(x[g == i].reshape(-1, K), np.squeeze(y_star[g == i].reshape(-1, 1)))
            .coef_
        )

    return beta


def _grouped_interactive_effects_iteration_hetrogeneous(y, x, G, GF, N, T, K, kappa, gamma, tol, max_iterations):
    last_objective_value = np.inf
    F, Lambda = _get_factors_initial(y, GF, T)
    g = _get_clusters_initial(Lambda, G, GF)
    beta = _estimate_beta_initial_hetrogeneous(y, x, g, K, G)
    F, Lambda = _get_factors_hetrogeneous(y, x, beta, g, G, GF, T, K)

    obj_val_store = 5
    objective_values = np.zeros(obj_val_store)

    for i in range(max_iterations):
        g, objective_value = _get_clusters_hetrogeneous(y, x, beta, Lambda, g, G, F, GF, N, T)
        F, Lambda = _get_factors_hetrogeneous(y, x, beta, g, G, GF, T, K)
        beta = _estimate_beta_hetrogeneous(y, x, g, GF, F, Lambda, K, N, T, G, kappa, gamma)

        objective_values[i % obj_val_store] = objective_value
        if objective_values.max() - objective_values.min() < tol:
            break

        last_objective_value = objective_value

    return beta, g, F, Lambda, last_objective_value


def _reorder_groups(g, beta, F):
    """Reorders the groups based on the first value of alpha"""
    # FIXME this is not the best way to do this
    # But it works for now
    mapping = np.argsort(beta[0, :])
    ordered_g = np.argsort(mapping)[g]
    ordered_beta = beta[:, mapping]
    ordered_F = F[:, mapping]
    return ordered_g, ordered_beta, ordered_F


def grouped_interactive_effects_hetrogeneous(
    y, x, G, GF=None, kappa=0.0, gamma=3.7, tol=1e-6, gife_iterations=100, max_iterations=1000
):
    """Runs GIFE regression"""
    N, T, K = x.shape
    if GF is None:
        GF = np.array([1] * G)
    else:
        GF = np.array(GF)  # Ensures that is an np array

    best_objective_value = np.inf
    best_g = None
    best_beta = None
    best_F = None
    best_Lambda = None

    # FIXME Lamda returns all possible factors and not only the ones that are used

    for i in range(gife_iterations):
        beta, g, F, Lambda, objective_value = _grouped_interactive_effects_iteration_hetrogeneous(
            y, x, G, GF, N, T, K, kappa, gamma, tol, max_iterations
        )

        if objective_value < best_objective_value:
            best_objective_value = objective_value
            best_g = g
            best_beta = beta
            best_F = F
            best_Lambda = Lambda

    ordered_g, ordered_beta, ordered_F = _reorder_groups(best_g, best_beta, best_F)

    return ordered_beta, ordered_g, ordered_F, best_Lambda, best_objective_value
