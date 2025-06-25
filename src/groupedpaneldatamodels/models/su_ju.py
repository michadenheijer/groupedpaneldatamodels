import numpy as np
import scipy.optimize as opt

from numpy.linalg import lstsq, eigh, eigvalsh
from numba import njit
from sklearn.cluster import KMeans


@njit(fastmath=True)
def _base(y, x, beta, N, T, R):
    res = (y - np.sum(x * beta.T[:, None, :], axis=2)).T
    v_res = (res @ res.T) / N
    return eigvalsh(v_res)[:-R].sum() / T


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
def _objective_value_without_individual_effects(y, x, beta, alpha, kappa, N, R, T):
    base = _base(y, x, beta, N, T, R)
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


def _compute_resid(y, x, beta, alpha, factors, lambdas, N, T):
    diff = beta[:, :, None] - alpha[:, None, :]
    dist2 = np.sum(diff**2, axis=0)
    nearest = np.argmin(dist2, axis=1)
    beta_rounded = alpha[:, nearest]

    return y - np.sum(x * beta_rounded.T[:, None, :], axis=2) - (factors @ lambdas).T


def interactive_effects_estimation(
    y: np.ndarray,
    x: np.ndarray,
    N: int,
    T: int,
    K: int,
    G: int,
    R: int,
    max_iter: int = 1000,
    only_bfgs: bool = True,
    tol: float = 1e-6,
    kappa: float = 0.1,
):
    """Internal function to estimate the interactive effects model as described by Su and Ju (2018).

    Args:
        y (np.ndarray): Dependent variable, shape (N, T, 1).
        x (np.ndarray): Explanatory variables, shape (N, T, K).
        N (int): Number of individuals (cross-sectional units).
        T (int): Number of time periods.
        K (int): Number of explanatory variables.
        G (int): Number of groups.
        R (int): Number of common factors (global).
        max_iter (int, optional): Maximum number of acceptable iterations, may be too large. Defaults to 1000.
        only_bfgs (bool, optional): Only uses BFGS. Defaults to True.
        tol (float, optional): Acceptable tolerance for the stopping condition. Defaults to 1e-6.
        kappa (float, optional): Kappa penalty parameter. Defaults to 0.1.

    Returns:
        tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
            - beta: Estimated coefficients for each individual, shape (K, N)
            - alpha: Ordered group-level representative coefficients, shape (K, G)
            - lambdas: Factor loadings, shape (R, N)
            - factors: Common factors, shape (T, R)
            - resid: Residuals of the model, shape (N, T)
    """
    y = np.squeeze(y, axis=2)
    beta, alpha, _ = _generate_initial_estimates(y, x, N, T, K, G)

    alpha_prev = alpha.copy()
    obj_value = np.inf
    last_obj_value = np.inf
    for i in range(max_iter):
        for j in range(G):
            alpha_fixed = alpha.copy()

            def unpack_local(theta):
                beta = theta[: K * N].reshape(K, N)
                alpha = alpha_fixed.copy()
                alpha[:, j : j + 1] = theta[K * N :].reshape(K, 1)
                return beta, alpha

            def obj_local(theta):
                beta, alpha = unpack_local(theta)
                return _objective_value_without_individual_effects(y, x, beta, alpha, kappa, N, R, T)

            def pack_local(beta, alpha):
                return np.concatenate((beta.flatten(), alpha[:, j].flatten()), axis=0)

            if i % 2 == 0 or only_bfgs:
                minimizer = opt.minimize(
                    obj_local,
                    pack_local(beta, alpha),
                    method="L-BFGS-B",
                    options={"maxiter": 100 if only_bfgs else 10},
                    tol=1e-4,
                )
                beta, alpha = unpack_local(minimizer.x)
                obj_value = minimizer.fun
                # print(f"BFGS Iteration {i}, Group {j}, Objective Value: {obj_value:.6f}")
            else:
                minimizer = opt.minimize(
                    obj_local,
                    pack_local(beta, alpha),
                    method="Nelder-Mead",
                    options={"adaptive": True, "maxiter": 100},
                    tol=1e-6,
                )
                beta, alpha = unpack_local(minimizer.x)
                obj_value = minimizer.fun
                # print(f"Nelder-Mead Iteration {i}, Group {j}, Objective Value: {obj_value:.6f}")

        if (
            np.abs(obj_value - last_obj_value) < tol
            and np.linalg.norm(alpha - alpha_prev) / np.linalg.norm(alpha_prev) < tol
        ):
            # print("Convergence reached.")
            break

        last_obj_value = obj_value
        alpha_prev = alpha.copy()

    res = (np.squeeze(y) - np.sum(x * beta.T[:, None, :], axis=2)).T
    res_var = (res @ res.T) / N
    factors = eigh(res_var).eigenvectors[:, -R:]
    factors = factors[:, ::-1]  # Reverse to have descending order

    lambdas = np.zeros((R, N), dtype=np.float32)

    for i in range(R):
        lambdas[i, :] = factors[:, i].T @ res

    resid = _compute_resid(y, x, beta, alpha, factors, lambdas, N, T)

    return beta, _order_alpha(alpha), lambdas, factors, resid
