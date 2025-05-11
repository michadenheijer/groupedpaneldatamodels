# bonhomme_manresa.pyx
# cython: profile=False
# cython: language_level=3
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True

"""
Grouped-fixed-effects (GFE) estimator – Cython version.
The core numerical work is still done with vectorised NumPy routines,
but the tight Python loops and attribute look-ups are eliminated.
"""

import numpy as _np
cimport numpy as np
cimport cython
from libc.math cimport INFINITY

ctypedef np.float64_t DTYPE_t     # element type we will use everywhere


# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------
cdef inline tuple _get_starting_values(
        np.ndarray[DTYPE_t, ndim=3] y,
        np.ndarray[DTYPE_t, ndim=3] x,
        int G, int N, int K):
    """
    Draws random rows and solves a tiny OLS to obtain a first-step θ̂ and α̂.
    """
    cdef int num_start_vars = x.shape[1] + G    # same logic as original code

    # Random draws for θ
    cdef np.ndarray[np.intp_t, ndim=1] idx_theta = _np.random.choice(
        N, num_start_vars, replace=False).astype(_np.intp)

    cdef np.ndarray[DTYPE_t, ndim=2] X0 = x[idx_theta].reshape(-1, K)
    cdef np.ndarray[DTYPE_t, ndim=2] Y0 = y[idx_theta].reshape(-1, 1)
    cdef np.ndarray[DTYPE_t, ndim=2] theta_init = _np.linalg.lstsq(X0, Y0, rcond=None)[0]

    # Random draws for α
    cdef np.ndarray[np.intp_t, ndim=1] idx_alpha = _np.random.choice(
        N, G, replace=False).astype(_np.intp)

    # Y[idx] − X[idx] @ θ
    cdef np.ndarray[DTYPE_t, ndim=2] alpha_init = _np.squeeze(
        y[idx_alpha] - _np.matmul(x[idx_alpha], theta_init))

    return theta_init, alpha_init


cdef inline np.ndarray[np.intp_t, ndim=1] _compute_groupings(
        np.ndarray[DTYPE_t, ndim=2] res,
        np.ndarray[DTYPE_t, ndim=2] alpha):
    """
    For every observation choose the closest α (Euclidean distance).
    """
    # res shape: (N, K)  -> add singleton axis so broadcasting works
    # alpha shape: (G, K)
    cdef np.ndarray[DTYPE_t, ndim=2] dist = (
        (res[None, :, :] - alpha[:, None, :]) ** 2).sum(axis=2)

    return _np.argmin(dist, axis=0)      # (N,) int array


cdef inline np.ndarray[DTYPE_t, ndim=2] _compute_alpha(
        np.ndarray[DTYPE_t, ndim=2] res,
        np.ndarray[np.intp_t, ndim=1] g,
        int G):
    """
    α̂_g = average residual in group g.
    """
    cdef np.ndarray[np.intp_t, ndim=1] counts = _np.bincount(
        g, minlength=G).astype(_np.intp)          # (G,)

    cdef np.ndarray[DTYPE_t, ndim=2] sums = _np.zeros(
        (G, res.shape[1]), dtype=_np.float64)              # (G, K)

    _np.add.at(sums, g, res)                               # in-place accumulation
    return sums / counts[:, None]                          # broadcast along columns

# ---------------------------------------------------------------------------
# Fast α̂ update using explicit C loops (no np.add.at)
# ---------------------------------------------------------------------------
@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline np.ndarray[DTYPE_t, ndim=2] _compute_alpha_fast(
        np.ndarray[DTYPE_t, ndim=2] res,
        np.ndarray[np.intp_t, ndim=1] g,
        int G):
    """
    Faster α̂ update using tight nested loops to accumulate group means.
    """
    cdef int N = res.shape[0]
    cdef int L = res.shape[1]
    cdef np.ndarray[DTYPE_t, ndim=2] sums = _np.zeros((G, L), dtype=res.dtype)
    cdef np.ndarray[np.int64_t, ndim=1] counts = _np.zeros(G, dtype=_np.int64)
    cdef int i, l, gi
    for i in range(N):
        gi = <int>g[i]
        counts[gi] += 1
        for l in range(L):
            sums[gi, l] += res[i, l]
    for gi in range(G):
        if counts[gi]:
            sums[gi, :] /= counts[gi]
    return sums

cdef inline np.ndarray[DTYPE_t, ndim=2] _compute_theta(
        np.ndarray[DTYPE_t, ndim=3] x,
        np.ndarray[DTYPE_t, ndim=3] y,
        np.ndarray[DTYPE_t, ndim=2] alpha,
        np.ndarray[np.intp_t, ndim=1] g):
    """
    θ̂ from pooled OLS of (y − α_g) on x.
    """
    cdef int K = x.shape[2]
    # broadcast α_g so that (N,L,1) − (N,L,1) is well‑defined
    cdef np.ndarray[DTYPE_t, ndim=2] theta = _np.linalg.lstsq(
        x.reshape(-1, K),
        (y - alpha[g][:, :, None]).reshape(-1, 1),
        rcond=None
    )[0]
    return theta


# ---------------------------------------------------------------------------
# Fast θ̂ update using precomputed X'X and X'
# ---------------------------------------------------------------------------
cdef inline np.ndarray[DTYPE_t, ndim=2] _compute_theta_fast(
        np.ndarray[DTYPE_t, ndim=2] XtX,
        np.ndarray[DTYPE_t, ndim=2] Xt,
        np.ndarray[DTYPE_t, ndim=3] y,
        np.ndarray[DTYPE_t, ndim=2] alpha,
        np.ndarray[np.intp_t, ndim=1] g):
    """
    Faster θ̂ update using precomputed X'X and X'.
    """
    cdef np.ndarray[DTYPE_t, ndim=2] XtY = _np.matmul(
        Xt, (y - alpha[g][:, :, None]).reshape(-1, 1))
    return _np.linalg.solve(XtX, XtY)


cdef inline np.ndarray[DTYPE_t, ndim=2] _compute_residuals(
        np.ndarray[DTYPE_t, ndim=3] y,
        np.ndarray[DTYPE_t, ndim=3] x,
        np.ndarray[DTYPE_t, ndim=2] theta):
    """
    r_i = y_i − x_i θ̂
    """
    return _np.squeeze(y - _np.matmul(x, theta))


cdef inline double _compute_objective_value(
        np.ndarray[DTYPE_t, ndim=2] res,
        np.ndarray[DTYPE_t, ndim=2] alpha,
        np.ndarray[np.intp_t, ndim=1] g):
    """
    Sum of squared residuals.
    """
    return float(((res - alpha[g]) ** 2).sum())


# -----------------------------------------------------------------------------
# One full inner cycle of the algorithm
# -----------------------------------------------------------------------------

cdef tuple _gfe_iteration(
        np.ndarray[DTYPE_t, ndim=3] y,
        np.ndarray[DTYPE_t, ndim=3] x,
        int G, int N, int K,
        int max_iter, double tol):
    """
    Core fixed-point iterations (one random initialisation).
    """
    cdef np.ndarray[DTYPE_t, ndim=2] theta, alpha
    cdef np.ndarray[np.intp_t, ndim=1] g
    cdef np.ndarray[DTYPE_t, ndim=2] res
    cdef double old_obj = INFINITY, new_obj
    cdef int it_used = 0

    # Pre‑compute X'X and X' once: they depend only on x
    cdef np.ndarray[DTYPE_t, ndim=2] X_flat = x.reshape(-1, K)
    cdef np.ndarray[DTYPE_t, ndim=2] Xt = X_flat.T
    cdef np.ndarray[DTYPE_t, ndim=2] XtX = _np.matmul(Xt, X_flat)

    theta, alpha = _get_starting_values(y, x, G, N, K)
    res = _compute_residuals(y, x, theta)
    g = _compute_groupings(res, alpha)

    for it_used in range(max_iter):
        alpha = _compute_alpha_fast(res, g, G)
        theta = _compute_theta_fast(XtX, Xt, y, alpha, g)
        res = _compute_residuals(y, x, theta)
        alpha = _compute_alpha_fast(res, g, G)
        g = _compute_groupings(res, alpha)

        new_obj = _compute_objective_value(res, alpha, g)
        if abs(old_obj - new_obj) < tol:
            break
        old_obj = new_obj

    return theta, alpha, g, it_used, new_obj


# -----------------------------------------------------------------------------
# Public function – exactly the same signature you had before
# -----------------------------------------------------------------------------
def grouped_fixed_effects(np.ndarray[DTYPE_t, ndim=3] y,
                          np.ndarray[DTYPE_t, ndim=3] x,
                          int G,
                          int max_iter=10000,
                          double tol=1e-8,
                          int gfe_iterations=20):
    """
    Run several random starts and return the best solution.
    """
    cdef int N = x.shape[0]
    cdef int K = x.shape[2]

    cdef np.ndarray[DTYPE_t, ndim=2] best_theta = None
    cdef np.ndarray[DTYPE_t, ndim=2] best_alpha = None
    cdef np.ndarray[np.intp_t, ndim=1] best_g = None
    cdef int best_it = -1
    cdef double best_obj = INFINITY

    cdef tuple res_tup
    for _ in range(gfe_iterations):
        res_tup = _gfe_iteration(y, x, G, N, K, max_iter, tol)
        if res_tup[4] < best_obj:
            best_theta, best_alpha, best_g, best_it, best_obj = res_tup

    return best_theta, best_alpha, best_g, best_it, best_obj