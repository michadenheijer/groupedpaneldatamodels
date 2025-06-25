from numba import njit
import numpy as np


@njit
def _fast_qr_f32(A, b, tol):
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
def superfast_lstsq(A, b, tol=1e-5, disable_superfast=True):
    """An internal function to compute least-squares solutions using a faster QR factorization,
    falls back on regular `lstsq` if QR factorization .

    Args:
        A (ArrayLike): Left-hand side matrix (design matrix).
        b (ArrayLike): Right-hand side vector or matrix (observations).
        tol (float, optional): Acceptable tolerance for QR-factorization. Defaults to 1e-5.
        disable_superfast (bool, optional): Disables QR-factorization. Defaults to True.

    Returns:
        ArrayLike: Least-squares solution to the equation Ax = b.
    """
    if not disable_superfast:
        X, ok = _fast_qr_f32(A, b, tol)
        if ok:  # QR succeeded
            return X.ravel() if b.ndim == 1 else X

    # robust fallback (also works for empty A)
    Xf, *_ = np.linalg.lstsq(A.astype(np.float32), b.astype(np.float32), rcond=-1)
    return Xf
