import numpy as np
from scipy.spatial.distance import pdist, squareform
from dtw import dtw_matrix


def distance_matrix(X, method='dtw', window=1):
    X = np.asarray(X, dtype=float)
    if X.ndim != 2 or not np.isfinite(X).all():
        raise ValueError('A finite gene-by-time matrix is required')
    if window < 0:
        raise ValueError('window must be nonnegative')
    if method == 'dtw':
        D = dtw_matrix(X, window)
    elif method == 'euclidean':
        D = squareform(pdist(X, 'euclidean')).astype(np.float32)
    elif method == 'pearson':
        C = X-X.mean(axis=1, keepdims=True)
        norm = np.linalg.norm(C, axis=1)
        C /= np.where(norm > 1e-12, norm, 1)[:, None]
        D = np.clip(1-C@C.T, 0, 2).astype(np.float32)
        # Undefined correlation: constant vs nonconstant = 1; both constants = 0.
        constant = norm <= 1e-12
        D[np.ix_(constant, constant)] = 0
    else:
        raise ValueError('Unknown distance method')
    np.fill_diagonal(D, 0)
    return D
