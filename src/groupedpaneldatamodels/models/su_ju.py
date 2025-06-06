import numpy as np
import scipy.optimize as opt

from numpy.linalg import lstsq, eigh, eigvalsh
from sklearn.cluster import KMeans


def _objective_value_without_individual_effects(y, x, beta, alpha, kappa, N, R, T):
    y = np.squeeze(y, axis=2)  # FIXME preferably this should be done outside
    res = (y - np.sum(x * beta.T[:, None, :], axis=2)).T
    v_res = (res @ res.T) / N
    base = eigvalsh(v_res)[:-R].sum() / T
    penalty = np.mean(np.prod(np.linalg.norm(beta[:, :, None] - alpha[:, None, :], axis=0), axis=1)) * kappa
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


def order_alpha(alpha):
    """Reorders the groups based on the first value of alpha"""
    # FIXME this is not the best way to do this
    # But it works for now
    mapping = np.argsort(alpha[0])
    ordered_alpha = alpha[:, mapping]
    return ordered_alpha


def interactive_effects_estimation(y, x, N, T, K, G, R, max_iter=1000, only_bfgs=False, tol=1e-6, kappa=10):
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
                    obj_local, pack_local(beta, alpha), method="BFGS", options={"maxiter": 10}, tol=1e-6
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

    lambdas = np.zeros((R, N))

    for i in range(R):
        lambdas[i, :] = factors[:, i].T @ res

    return beta, order_alpha(alpha), lambdas, factors
