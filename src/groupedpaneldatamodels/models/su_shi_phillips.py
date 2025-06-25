import numpy as np
import scipy.optimize as opt

from numba import njit
from numpy.linalg import lstsq
from sklearn.cluster import KMeans


@njit(fastmath=True)
def _norm(beta, alpha):
    alpha_norm2 = np.sum(alpha * alpha, axis=0)  # (n,)
    beta_norm2 = np.sum(beta * beta, axis=0)  # (m,)

    # all pairwise dot products
    dot = beta.T @ alpha  # (m, n)

    # expand to (m,1) and (1,n) for broadcasting
    d2 = beta_norm2.reshape((-1, 1)) + alpha_norm2.reshape((1, -1)) - 2.0 * dot

    # clamp tiny negatives from round-off
    d2 = np.maximum(d2, 0.0)

    return np.sqrt(d2)


@njit(fastmath=True)
def _row_prod(a):
    """Row-wise product for 2-D array `a` – returns shape (a.shape[0],)."""
    m, n = a.shape
    out = np.empty(m, dtype=a.dtype)
    for i in range(m):
        p = 1.0
        for j in range(n):
            p *= a[i, j]
        out[i] = p
    return out


@njit(fastmath=True)
def _objective_value(y, x, beta, alpha, mu, kappa):
    base = ((y - np.sum(x * beta.T[:, None, :], axis=2) - mu) ** 2).mean()
    penalty = np.mean(_row_prod(_norm(beta, alpha))) * kappa
    return base + penalty


@njit(fastmath=True)
def _objective_value_without_individual_effects(y, x, beta, alpha, kappa):
    base = ((y - np.sum(x * beta.T[:, None, :], axis=2)) ** 2).mean()
    penalty = np.mean(_row_prod(_norm(beta, alpha))) * kappa
    return base + penalty


def _generate_initial_estimates(y, x, N, T, K, G):
    beta_init = np.zeros((K, N))

    for i in range(N):
        beta_init[:, i : i + 1] = lstsq(x[i].reshape(T, K), y[i].reshape(T, 1))[0]
    alpha_init = KMeans(n_clusters=G).fit(beta_init.T).cluster_centers_.T

    for j in range(G):
        if np.abs(beta_init.T - alpha_init[:, j]).min() < 1e-2:
            alpha_init[:, j] += 1e-1 * np.sign(alpha_init[:, j])

    mu_init = np.mean(y, axis=1)

    return beta_init, alpha_init, mu_init


def _order_alpha(alpha):
    """Reorders the groups based on the first value of alpha"""
    # FIXME this is not the best way to do this
    # But it works for now
    mapping = np.argsort(alpha[0])
    ordered_alpha = alpha[:, mapping]
    return ordered_alpha


def _compute_resid(y, x, beta, alpha, mu):
    diff = beta[:, :, None] - alpha[:, None, :]
    dist2 = np.sum(diff**2, axis=0)
    nearest = np.argmin(dist2, axis=1)
    beta_rounded = alpha[:, nearest]

    return y - np.sum(x * beta_rounded.T[:, None, :], axis=2) - mu


def fixed_effects_estimation(
    y, x, N, T, K, G, max_iter=1000, only_bfgs=True, tol=1e-6, use_individual_effects=True, kappa=0.1
):
    """Internal estimation function for Su, Shi and Phillips (2016) model.

    Args:
        y (np.ndarray): Dependent variable, shape (N, T).
        x (np.array): Explanatory variables, shape (N, T, K).
        N (int): Number of individuals (cross-sectional units).
        T (int): Number of time periods.
        K (int): Number of explanatory variables.
        G (int): Number of groups.
        max_iter (int, optional): Maximum number of iterations. Defaults to 1000.
        only_bfgs (bool, optional): Only uses BFGS, instead of Nelder-Mead. Defaults to True.
        tol (float, optional): Acceptable tolerance before stopping the estimation. Defaults to 1e-6.
        use_individual_effects (bool, optional): Enables indvidual effects. Defaults to True.
        kappa (float, optional): Penalization parameter. Defaults to 0.1.

    Returns:
        tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
            - beta: Estimated coefficients for each individual, shape (K, N)
            - mu: Estimated individual fixed effects, shape (N, 1)
            - alpha: Group-level representative coefficients, shape (K, G)
            - resid: Residuals of the model, shape (N, T)
    """
    # FIXME this works, though the code is not very clean
    if use_individual_effects:
        y_demeaned = y - np.mean(y, axis=1, keepdims=True)
        x_demeaned = x - np.mean(x, axis=1, keepdims=True)

        beta, mu, alpha, resid = fixed_effects_estimation(
            y_demeaned, x_demeaned, N, T, K, G, max_iter, only_bfgs, tol, False, kappa
        )

        mu = np.mean(y - np.sum(x * beta.T[:, None, :], axis=2), axis=1, keepdims=True)
        return beta, mu, _order_alpha(alpha), resid

    beta, alpha, mu = _generate_initial_estimates(y, x, N, T, K, G)

    obj_value = np.inf
    last_obj_value = np.inf

    for i in range(max_iter):
        for j in range(G):
            alpha_fixed = alpha.copy()

            def unpack_local(theta):
                if use_individual_effects:
                    beta = theta[: K * N].reshape(K, N)
                    mu = theta[K * N : K * N + N].reshape(N, 1)
                    alpha = alpha_fixed.copy()
                    alpha[:, j : j + 1] = theta[K * N + N :].reshape(K, 1)
                    return beta, mu, alpha
                else:
                    beta = theta[: K * N].reshape(K, N)
                    alpha = alpha_fixed.copy()
                    alpha[:, j : j + 1] = theta[K * N :].reshape(K, 1)
                    return beta, None, alpha

            def obj_local(theta):
                beta, mu, alpha = unpack_local(theta)

                if use_individual_effects:
                    return _objective_value(y, x, beta, alpha, mu, kappa)

                else:
                    return _objective_value_without_individual_effects(y, x, beta, alpha, kappa)

            def pack_local(beta, mu, alpha):
                if use_individual_effects:
                    return np.concatenate((beta.flatten(), mu.flatten(), alpha[:, j].flatten()), axis=0)
                else:
                    return np.concatenate((beta.flatten(), alpha[:, j].flatten()), axis=0)

            if i % 2 == 0 or only_bfgs:
                minimizer = opt.minimize(
                    obj_local,
                    pack_local(beta, mu, alpha),
                    method="L-BFGS-B",  #  FIXME BFGS may be better, but slower
                    options={"maxiter": 100 if only_bfgs else 10},
                    tol=tol,
                )
                beta, mu, alpha = unpack_local(minimizer.x)
                obj_value = minimizer.fun
                # print(f"BFGS Iteration {i}, Group {j}, Objective Value: {obj_value:.6f}")
            else:
                minimizer = opt.minimize(
                    obj_local,
                    pack_local(beta, mu, alpha),
                    method="Nelder-Mead",
                    options={"adaptive": True, "maxiter": 100},
                    tol=tol,
                )
                beta, mu, alpha = unpack_local(minimizer.x)
                obj_value = minimizer.fun
                # print(f"Nelder-Mead Iteration {i}, Group {j}, Objective Value: {obj_value:.6f}")

        # TODO fix this, because does not fully function
        # NOTE added i % 2 == 0 to avoid convergence too early
        if np.abs(obj_value - last_obj_value) < tol and (i % 2 == 0 or only_bfgs):
            # print("Convergence reached.")
            break

        last_obj_value = obj_value

    if use_individual_effects:
        resid = _compute_resid(y, x, beta, alpha, mu)
        assert mu is not None, "Mu should not be None when use_individual_effects is True"
        return beta, mu, _order_alpha(alpha), resid
    else:
        mu = np.zeros((N, 1))
        resid = _compute_resid(y, x, beta, alpha, mu)
        return beta, mu, _order_alpha(alpha), resid
