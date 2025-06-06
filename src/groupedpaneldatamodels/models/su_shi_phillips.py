import numpy as np
import scipy.optimize as opt

from numpy.linalg import lstsq
from sklearn.cluster import KMeans


def _objective_value(y, x, beta, alpha, mu, kappa):
    base = ((y - np.sum(x * beta.T[:, None, :], axis=2) - mu) ** 2).mean()
    penalty = np.mean(np.prod(np.linalg.norm(beta[:, :, None] - alpha[:, None, :], axis=0), axis=1)) * kappa
    return base + penalty


def _objective_value_without_individual_effects(y, x, beta, alpha, kappa):
    base = ((y - np.sum(x * beta.T[:, None, :], axis=2)) ** 2).mean()
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


def fixed_effects_estimation(
    y, x, N, T, K, G, max_iter=1000, only_bfgs=True, tol=1e-6, use_individual_effects=True, kappa=10
):
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
                    tol=1e-4,
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
                    tol=1e-4,
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
        return beta, mu, order_alpha(alpha)
    else:
        mu = np.zeros((N, 1))
        return beta, mu, order_alpha(alpha)
