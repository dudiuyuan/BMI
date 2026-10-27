"""Self-implemented squared-cost DTW, accelerated by Numba, not a DTW library."""
import numpy as np
from numba import njit


@njit(cache=True)
def _dtw(x, y, window):
    n, m = len(x), len(y)
    prev = np.full(m + 1, np.inf)
    prev[0] = 0.
    for i in range(1, n + 1):
        curr = np.full(m + 1, np.inf)
        for j in range(max(1, i - window), min(m, i + window) + 1):
            curr[j] = (x[i-1] - y[j-1]) ** 2 + min(prev[j], curr[j-1], prev[j-1])
        prev = curr
    return np.sqrt(prev[m])


def dtw(x, y, window=None):
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if x.ndim != 1 or y.ndim != 1 or not len(x) or not len(y):
        raise ValueError('Nonempty one-dimensional sequences required')
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError('DTW requires finite data')
    w = max(len(x), len(y)) if window is None else int(window)
    if w < abs(len(x)-len(y)) or w < 0:
        raise ValueError('Window cannot connect the endpoints')
    return float(_dtw(x, y, w))


@njit(cache=True)
def dtw_matrix(X, window):
    n = len(X)
    D = np.zeros((n, n), dtype=np.float32)
    for i in range(n):
        for j in range(i+1, n):
            D[i, j] = _dtw(X[i], X[j], window)
            D[j, i] = D[i, j]
    return D
