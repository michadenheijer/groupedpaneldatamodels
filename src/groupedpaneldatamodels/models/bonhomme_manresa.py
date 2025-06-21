# TODO
# - Check imports


import numpy as np
from numpy.linalg import lstsq
from scipy.linalg import lstsq as lstsq_sp

from numba import njit

# NOTE ignore RuntimeWarning for now
import warnings

warnings.simplefilter(action="ignore", category=RuntimeWarning)


from numba import njit
import numpy as np


import numpy as np
from numba import njit


import numpy as np
from numba import njit


import numpy as np
from numba import njit


import numpy as np
from numba import njit


@njit
def fast_qr_f32(A, b, tol):
    """
    Float-32 QR least-squares with empty-matrix guard.

    Always returns a 2-D array (n, k) so Numba sees a single return type.
    """
    # Promote to float32 and make RHS 2-D -------------------------------
    A32 = A.astype(np.float32)
    b32 = b.astype(np.float32).reshape(-1, 1) if b.ndim == 1 else b.astype(np.float32)

    m, n = A32.shape
    k = b32.shape[1]  # number of RHS

    # --------- EARLY EXIT for empty design or empty sample -------------
    if m == 0 or n == 0:
        return np.zeros((n, k), dtype=np.float32), False

    # ----------------------- QR factorisation --------------------------
    Q, R = np.linalg.qr(A32)

    # Cheap full-rank test ---------------------------------------------
    rdiag_ok = True
    rcond = tol
    rmin = n if m >= n else m
    for i in range(rmin):
        if abs(R[i, i]) <= rcond:
            rdiag_ok = False
            break

    if rdiag_ok:  # fast path
        X = np.linalg.solve(R, np.ascontiguousarray(Q.T) @ np.ascontiguousarray(b32))  # (n, k)
        return X, True

    # fallback signal – same array shape/type ---------------------------
    return np.empty((n, k), dtype=np.float32), False


# NOTE disable superfast_lstsq for now, as it is untested
def superfast_lstsq(A, b, tol=1e-5):
    # X, ok = fast_qr_f32(A, b, tol)
    # if ok:  # QR succeeded
    #     return X.ravel() if b.ndim == 1 else X

    # robust fallback (also works for empty A)
    Xf, *_ = np.linalg.lstsq(A.astype(np.float32), b.astype(np.float32), rcond=-1)
    return Xf


# @njit
def _get_starting_values(y, x, G: int, N: int, K: int):
    """Generates the starting values of theta"""
    num_start_vars: int = K + G  # FIXME I believe that shape is slow in Cython
    random_draws_theta = np.random.choice(N, num_start_vars, replace=False)
    x_stack_start = x[random_draws_theta].reshape(-1, K)
    y_stack_start = y[random_draws_theta].reshape(-1, 1)

    # FIXME some errors may arise, maybe add some checks
    theta_init = lstsq(x_stack_start, y_stack_start, rcond=None)[0]

    random_draws_alpha = np.random.choice(N, size=G, replace=False)
    alpha_init = np.squeeze(y[random_draws_alpha] - x[random_draws_alpha, :, :] @ theta_init)

    return theta_init, alpha_init


@njit
def _compute_groupings(res, alpha):
    """Computes the groupings based on the residuals and alpha"""
    euclidean_distance_between_grouping = ((res[None, :, :] - alpha[:, None, :]) ** 2).sum(axis=2)
    g = np.argmin(euclidean_distance_between_grouping, axis=0)  # Closest group
    return g


# @njit
def _compute_alpha(res, g, G):
    """Computes the alpha values based on the residuals and groupings"""
    counts = np.bincount(g, minlength=G)[:, None]  # (G, 1) — number of elements in each group
    sums = np.zeros((G, res.shape[1]))  # (G, K) — sum of residuals per group
    np.add.at(sums, g, res)  # sums[i] += res[j] for all j where g[j] == i
    alpha = sums / counts  # mean = sum / count
    return alpha


# @njit
def _compute_theta(x, y, alpha, g):
    """Computes the theta values based on the x, y, alpha and groupings"""
    # FIXME check if this makes sense
    K = x.shape[2]  # FIXME I believe that shape is slow in Cython
    # theta = lstsq_sp(
    #     x.reshape(-1, K), y.reshape(-1, 1) - alpha[g].reshape(-1, 1), lapack_driver="gelsy"
    # )[  # type:ignore
    #     0
    # ]
    theta = superfast_lstsq(
        x.reshape(-1, K),
        y.reshape(-1, 1) - alpha[g].reshape(-1, 1),
    )
    return theta


# @njit
def _compute_residuals(y, x, theta):
    """Computes the residuals based on y, x and theta"""
    res = np.squeeze(y - x @ theta)
    return res


@njit
def _compute_objective_value(res, alpha, g):
    """Computes the objective value based on the residuals, alpha and groupings"""
    objective_value = ((res - alpha[g]) ** 2).sum()
    return objective_value


def _reorder_groups(g, alpha):
    """Reorders the groups based on the first value of alpha"""
    # FIXME this is not the best way to do this
    # But it works for now
    mapping = np.argsort(alpha[:, 0])
    ordered_g = np.argsort(mapping)[g]
    ordered_alpha = alpha[mapping]
    return ordered_g, ordered_alpha


def _grouped_fixed_effects_iteration(y, x, G: int, N: int, K: int, max_iter=10000, tol=1e-8):
    # FIXME possibly create some sort of array that stores these values
    # Could be used for debugging
    theta, alpha = _get_starting_values(y, x, G, N, K)
    res = _compute_residuals(y, x, theta)
    g = _compute_groupings(res, alpha)

    objective_value = np.inf

    iterations_used = 0
    for i in range(max_iter):
        alpha = _compute_alpha(res, g, G)
        theta = _compute_theta(x, y, alpha, g)
        res = _compute_residuals(y, x, theta)
        alpha = _compute_alpha(res, g, G)
        g = _compute_groupings(res, alpha)
        new_objective_value = _compute_objective_value(res, alpha, g)

        if abs(objective_value - new_objective_value) < tol:
            iterations_used = i
            objective_value = new_objective_value
            break

        objective_value = new_objective_value
    resid = _compute_residuals(y, x, theta)[np.arange(N), :, g] - alpha[g]

    return theta, alpha, g, iterations_used, objective_value, resid


def _compute_eta(y_bar, x_bar, theta):
    """Computes the eta values based on y_bar, x_bar and theta"""
    eta = np.squeeze(y_bar - x_bar @ theta)
    return eta


def _hkmeans(y, x, theta, alpha, g, G, max_iter=1000, tol=1e-6):
    objective_value = np.inf

    for i in range(max_iter):
        res = _compute_residuals(y, x, theta)
        alpha = _compute_alpha(res, g, G)
        theta = _compute_theta(x, y, alpha, g)
        res = _compute_residuals(y, x, theta)
        alpha = _compute_alpha(res, g, G)
        g = _compute_groupings(res, alpha)
        alpha = _compute_alpha(res, g, G)
        new_objective_value = _compute_objective_value(res, alpha, g)

        if abs(objective_value - new_objective_value) < tol:
            return theta, alpha, g, i, new_objective_value

        objective_value = new_objective_value

    return theta, alpha, g, max_iter, objective_value


def _run_vns(y, x, g, G, N, alpha, theta, init_objective_value, max_vns_iter=10, tol=1e-8, max_alg1_iter=20):
    # FIXME this is not the best way to do this
    # But it works for now
    best_objective_value = init_objective_value
    objective_value = np.inf
    g = g.copy()

    i = 1
    while i <= max_vns_iter:
        # Randomly change a few groupings
        g_new = g.copy()

        # FIXME this should check if there are empty groups
        g_new[np.random.choice(N, size=i, replace=False)] = np.random.choice(G, size=i, replace=True)

        # Apply algorithm 1
        theta_new, alpha_new, g_new, iterations_used, objective_value = _hkmeans(y, x, theta, alpha, g_new, G)

        # Local 1 step search
        changed = True
        while changed:
            changed = False
            for j in range(N):
                for k in range(G):
                    if g_new[j] == k:
                        continue

                    # FIXME add check to not leave any group empty

                    g_local = g_new.copy()
                    g_local[j] = k

                    # FIXME this may be very slow
                    if (np.bincount(g_local, minlength=G) == 0).any():
                        continue

                    theta_local, alpha_local, g_local, iterations_used, objective_value_local = _hkmeans(
                        y, x, theta_new, alpha_new, g_local, G, max_alg1_iter
                    )

                    if objective_value_local < objective_value:
                        g_new = g_local
                        theta_new = theta_local
                        alpha_new = alpha_local
                        objective_value = objective_value_local
                        changed = True
                        # print(f"Changed group {j} to {k} in iteration {i} with objective value {objective_value}")

        if objective_value < best_objective_value:
            g = g_new
            theta = theta_new
            alpha = alpha_new
            best_objective_value = objective_value
            i = 1

        else:
            i += 1

    return g, theta, alpha, objective_value


def _grouped_fixed_effects_iteration_vns(y, x, G: int, N: int, K: int, max_iter=10000, tol=1e-8, neighbor_max=10):
    # FIXME possibly create some sort of array that stores these values
    # Could be used for debugging
    theta, alpha = _get_starting_values(y, x, G, N, K)
    res = _compute_residuals(y, x, theta)
    g = _compute_groupings(res, alpha)

    objective_value = np.inf

    for i in range(max_iter):
        alpha = _compute_alpha(res, g, G)
        theta = _compute_theta(x, y, alpha, g)
        res = _compute_residuals(y, x, theta)
        alpha = _compute_alpha(res, g, G)
        g = _compute_groupings(res, alpha)
        alpha = _compute_alpha(res, g, G)
        new_objective_value = _compute_objective_value(res, alpha, g)

        # TODO run VNS
        g, theta, alpha, new_objective_value = _run_vns(y, x, g, G, N, alpha, theta, new_objective_value)

        if abs(objective_value - new_objective_value) < tol:
            iterations_used = i
            objective_value = new_objective_value
            break

        objective_value = new_objective_value
    resid = _compute_residuals(y, x, theta)[np.arange(N), :, g] - alpha[g]

    return theta, alpha, g, max_iter, objective_value, resid


# FIXME not used right now but still neccesary
def compute_statistics(objective_value, N, T, K, G):
    """Computes the statistics based on the objective value, N, T, K and G"""
    # FIXME this is not the best way to do this
    # But it works for now
    sigma_squared = 1 / (N * T - G * T - N - K) * objective_value
    BIC = 1 / (N * T) * objective_value + sigma_squared * (G * T + N + K) / (N * T)
    return sigma_squared, BIC


# FIXME remove this variable
def _get_starting_values_hetrogeneous(y, x, G: int, N: int, K: int):
    """Generates the starting values of theta"""
    num_start_vars: int = x.shape[1] + G  # FIXME I believe that shape is slow in Cython
    random_draws_theta = np.random.choice(N, num_start_vars, replace=False)
    x_stack_start = x[random_draws_theta].reshape(-1, K)
    y_stack_start = y[random_draws_theta].reshape(-1, 1)

    # FIXME some errors may arise, maybe add some checks
    theta_init = lstsq(x_stack_start, y_stack_start, rcond=None)[0]

    random_draws_alpha = np.random.choice(N, size=G, replace=False)
    alpha_init = np.squeeze(y[random_draws_alpha] - x[random_draws_alpha, :, :] @ theta_init)

    # FIXME this should be changed such that
    return np.tile(theta_init, (1, G)), alpha_init


# TODO create function for hetrogeneous theta


# @njit
# def _compute_groupings_hetrogeneous(res, alpha):
#     """Computes the groupings based on the residuals and alpha"""
#     euclidean_distance_between_grouping = np.squeeze((res - alpha.T[None, :, :]) ** 2).sum(axis=1)
#     g = np.argmin(euclidean_distance_between_grouping, axis=1)  # Closest group
#     return g


def _compute_groupings_hetrogeneous(res, alpha):
    """Assign each unit to the nearest group (least-squares sense)."""
    dists = np.square(res - alpha.T[None, :, :]).sum(axis=1)  # (N, G)
    return np.argmin(dists, axis=1).astype(np.uint8)


# # @njit
def _compute_alpha_hetrogeneous(res, g, G):
    """Computes the alpha values based on the residuals and groupings"""
    n, T, _ = res.shape
    counts = np.bincount(g, minlength=G)  # shape: (G,)
    R = res[np.arange(n), :, g]  # shape: (n, T)
    sums = np.zeros((G, T))
    np.add.at(sums, g, R)  # sums[k] += R[i] whenever g[i]==k
    alphas = sums / counts[:, None]  # broadcast counts over T
    return alphas


# # @njit
def _compute_theta_hetrogeneous(x, y, alpha, g, G, K):
    """Computes the theta values based on the x, y, alpha and groupings"""
    # FIXME check if this makes sense
    theta = np.zeros((K, G))
    for i in range(G):
        # theta[:, i] = lstsq_sp(  # type:ignore
        #     x[g == i].reshape(-1, K),
        #     np.squeeze((np.squeeze(y[g == i]) - alpha[i]).reshape(-1, 1)),
        #     lapack_driver="gelsy",
        # )[0]
        theta[:, i] = superfast_lstsq(  # type:ignore
            x[g == i].reshape(-1, K),
            np.squeeze((np.squeeze(y[g == i]) - alpha[i]).reshape(-1, 1)),
        )

    return theta


# @njit
def _compute_residuals_hetrogeneous(y, x, theta, G):
    """Computes the residuals based on y, x and theta"""
    res = np.tile(y, (1, G)) - x @ theta
    return res


@njit
def _compute_objective_value_hetrogeneous(res, alpha, g, N):
    """Computes the objective value based on the residuals, alpha and groupings"""
    s = 0
    for i in range(N):
        s += ((res[i, :, g[i]] - alpha[g[i]]) ** 2).sum()

    return s


def _reorder_groups_hetrogeneous(g, alpha, theta, hetrogeneous_theta):
    """Reorders the groups based on the first value of alpha"""
    # FIXME this is not the best way to do this
    # But it works for now
    mapping = np.argsort(theta[0])
    ordered_g = np.argsort(mapping)[g]
    ordered_alpha = alpha[mapping]
    if not hetrogeneous_theta:
        return ordered_g, ordered_alpha, theta
    ordered_theta = theta[:, mapping]
    return ordered_g, ordered_alpha, ordered_theta


def _hkmeans_hetrogeneous(y, x, theta, alpha, g, G, K, N, max_iter=1000, tol=1e-6):
    objective_value = np.inf

    for i in range(max_iter):
        res = _compute_residuals_hetrogeneous(y, x, theta, G)
        alpha = _compute_alpha_hetrogeneous(res, g, G)
        theta = _compute_theta_hetrogeneous(x, y, alpha, g, G, K)
        res = _compute_residuals_hetrogeneous(y, x, theta, G)
        alpha = _compute_alpha_hetrogeneous(res, g, G)
        g = _compute_groupings_hetrogeneous(res, alpha)
        alpha = _compute_alpha_hetrogeneous(res, g, G)
        new_objective_value = _compute_objective_value_hetrogeneous(res, alpha, g, N)

        if abs(objective_value - new_objective_value) < tol:
            return theta, alpha, g, i, new_objective_value

        objective_value = new_objective_value

    return theta, alpha, g, max_iter, objective_value


def _run_vns_hetrogeneous(
    y, x, g, G, N, alpha, theta, init_objective_value, K, max_vns_iter=10, tol=1e-8, max_alg1_iter=20
):
    # FIXME this is not the best way to do this
    # But it works for now
    best_objective_value = init_objective_value
    objective_value = np.inf
    g = g.copy()

    i = 1
    while i <= max_vns_iter:
        # print(f"vns iter: {i}")
        # Randomly change a few groupings
        g_new = g.copy()

        # FIXME this should check if there are empty groups
        g_new[np.random.choice(N, size=i, replace=False)] = np.random.choice(G, size=i, replace=True)

        # Apply algorithm 1
        theta_new, alpha_new, g_new, iterations_used, objective_value = _hkmeans_hetrogeneous(
            y, x, theta, alpha, g_new, G, K, N
        )

        # Local 1 step search
        changed = True
        while changed:
            changed = False
            for j in range(N):
                for k in range(G):
                    if g_new[j] == k:
                        continue

                    # FIXME add check to not leave any group empty

                    g_local = g_new.copy()
                    g_local[j] = k

                    # FIXME this may be very slow
                    if (np.bincount(g_local, minlength=G) == 0).any():
                        continue

                    theta_local, alpha_local, g_local, iterations_used, objective_value_local = _hkmeans_hetrogeneous(
                        y, x, theta_new, alpha_new, g_local, G, K, N, max_alg1_iter
                    )

                    if objective_value_local < objective_value:
                        g_new = g_local
                        theta_new = theta_local
                        alpha_new = alpha_local
                        objective_value = objective_value_local
                        changed = True
                        # print(f"Changed group {j} to {k} in iteration {i} with objective value {objective_value}")

        if objective_value < best_objective_value:
            g = g_new
            theta = theta_new
            alpha = alpha_new
            best_objective_value = objective_value
            i = 1
            # print(f"Succes at iteration {i} with objective value {best_objective_value}")

        else:
            i += 1

    return g, theta, alpha, objective_value


# FIXME no stopping condition set up
def _grouped_fixed_effects_iteration_vns_hetrogeneous(
    y, x, G: int, N: int, K: int, max_iter=10000, tol=1e-8, neighbor_max=10
):
    # FIXME possibly create some sort of array that stores these values
    # Could be used for debugging
    theta, alpha = _get_starting_values_hetrogeneous(y, x, G, N, K)
    res = _compute_residuals_hetrogeneous(y, x, theta, G)
    g = _compute_groupings_hetrogeneous(res, alpha)

    objective_value = np.inf

    for i in range(max_iter):
        # print(f"ITERATION: {i}")
        alpha = _compute_alpha_hetrogeneous(res, g, G)
        theta = _compute_theta_hetrogeneous(x, y, alpha, g, G, K)
        res = _compute_residuals_hetrogeneous(y, x, theta, G)
        alpha = _compute_alpha_hetrogeneous(res, g, G)
        g = _compute_groupings_hetrogeneous(res, alpha)
        alpha = _compute_alpha_hetrogeneous(res, g, G)
        new_objective_value = _compute_objective_value_hetrogeneous(res, alpha, g, N)

        # TODO run VNS
        g, theta, alpha, new_objective_value = _run_vns_hetrogeneous(
            y, x, g, G, N, alpha, theta, new_objective_value, K
        )

        objective_value = new_objective_value
    resid = _compute_residuals_hetrogeneous(y, x, theta, G)[np.arange(N), :, g] - alpha[g]

    return theta, alpha, g, max_iter, objective_value, resid


def _grouped_fixed_effects_iteration_hetrogeneous(y, x, G: int, N: int, K: int, max_iter=10000, tol=1e-8):
    # FIXME possibly create some sort of array that stores these values
    # Could be used for debugging
    theta, alpha = _get_starting_values_hetrogeneous(y, x, G, N, K)
    res = _compute_residuals_hetrogeneous(y, x, theta, G)
    g = _compute_groupings_hetrogeneous(res, alpha)

    objective_value = np.inf

    iterations_used = 0
    for i in range(max_iter):
        alpha = _compute_alpha_hetrogeneous(res, g, G)
        theta = _compute_theta_hetrogeneous(x, y, alpha, g, G, K)
        res = _compute_residuals_hetrogeneous(y, x, theta, G)
        alpha = _compute_alpha_hetrogeneous(res, g, G)
        g = _compute_groupings_hetrogeneous(res, alpha)
        new_objective_value = _compute_objective_value_hetrogeneous(res, alpha, g, N)

        if abs(objective_value - new_objective_value) < tol:
            objective_value = new_objective_value
            break

        objective_value = new_objective_value
        iterations_used = i

    resid = _compute_residuals_hetrogeneous(y, x, theta, G)[np.arange(N), :, g] - alpha[g]

    return theta, alpha, g, iterations_used, objective_value, resid


def _compute_eta_hetrogeneous(y_bar, x_bar, theta):
    """Computes the eta values based on y_bar, x_bar and theta"""
    eta = np.squeeze(y_bar - x_bar @ theta.mean(axis=1))
    return eta


def grouped_fixed_effects(
    y,
    x,
    G,
    max_iter=10000,
    tol=1e-6,
    gfe_iterations=100,
    unit_specific_effects=False,
    enable_vns=False,
    hetrogeneous_theta=True,
):
    """
    Computes the grouped fixed effects using the algorithm described in the paper.
    """
    # FIXME not really required if unit_specific_effects is False
    # But this seems like the easiest implementation
    x_bar = np.mean(x, axis=1, keepdims=True)
    y_bar = np.mean(y, axis=1, keepdims=True)

    best_theta = None
    best_alpha = None
    best_g = None
    best_iterations_used = None
    best_objective_value = np.inf
    best_resid = None

    N = x.shape[0]
    T = x.shape[1]
    K = x.shape[2]

    # FIXME this code does not drop constant binary variables by itself
    # The input class should be able to do this
    if unit_specific_effects:
        # Demean x and y
        x = x - x_bar
        y = y - y_bar

    for _ in range(gfe_iterations):
        # FIXME should be better way to do this
        if hetrogeneous_theta and enable_vns:
            theta, alpha, g, iterations_used, objective_value, resid = (
                _grouped_fixed_effects_iteration_vns_hetrogeneous(
                    y,
                    x,
                    G,
                    N,
                    K,
                    max_iter,
                    tol,
                )
            )
        elif enable_vns:
            theta, alpha, g, iterations_used, objective_value, resid = _grouped_fixed_effects_iteration_vns(
                y,
                x,
                G,
                N,
                K,
                max_iter,
                tol,
            )

        elif hetrogeneous_theta:
            theta, alpha, g, iterations_used, objective_value, resid = _grouped_fixed_effects_iteration_hetrogeneous(
                y,
                x,
                G,
                N,
                K,
                max_iter,
                tol,
            )

        else:
            theta, alpha, g, iterations_used, objective_value, resid = _grouped_fixed_effects_iteration(
                y, x, G, N, K, max_iter, tol
            )

        if objective_value < best_objective_value:
            best_objective_value = objective_value
            best_g, best_alpha, best_theta = _reorder_groups_hetrogeneous(g, alpha, theta, hetrogeneous_theta)
            # best_g, best_alpha = g, alpha
            best_iterations_used = iterations_used
            best_resid = resid.reshape(-1)  # NOTE always reshape to 1d array

    # FIXME this should not be true for hetrogeneous theta
    # Because does not work
    assert best_theta is not None, "No theta found, something went wrong in the estimation process."
    assert best_g is not None, "No groupings found, something went wrong in the estimation process."
    if unit_specific_effects:
        eta = _compute_eta_hetrogeneous(y_bar, x_bar, best_theta)
        return best_theta, best_alpha, best_g, eta, best_iterations_used, best_objective_value, best_resid

    # NOTE None is for lack of unit specific effects
    return best_theta, best_alpha, best_g, None, best_iterations_used, best_objective_value, best_resid
