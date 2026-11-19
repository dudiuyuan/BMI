"""Alternating k-medoids; multiple starts, deterministic tie and empty-cluster handling."""
import numpy as np


def cluster(D, k, seed=11, n_init=5, max_iter=100):
    D = np.asarray(D)
    n = len(D)
    if D.shape != (n, n) or not np.isfinite(D).all() or np.any(D < 0):
        raise ValueError('Invalid distance matrix')
    if not np.allclose(D, D.T) or not np.allclose(np.diag(D), 0):
        raise ValueError('Distance matrix must be symmetric with zero diagonal')
    if not 1 <= k <= n or n_init < 1 or max_iter < 1:
        raise ValueError('Invalid k/n_init/max_iter')
    rng = np.random.default_rng(seed)
    best = None
    for _ in range(n_init):
        medoids = np.sort(rng.choice(n, k, replace=False))
        for iteration in range(max_iter):
            labels = D[:, medoids].argmin(axis=1)
            labels[medoids] = np.arange(k)  # each medoid owns itself even at zero-distance ties
            updated = medoids.copy()
            for c in range(k):
                idx = np.flatnonzero(labels == c)
                costs = D[np.ix_(idx, idx)].sum(axis=1, dtype=np.float64)
                updated[c] = idx[np.argmin(costs)]
            if np.array_equal(updated, medoids):
                break
            medoids = updated
        labels = D[:, medoids].argmin(axis=1)
        labels[medoids] = np.arange(k)
        objective = float(D[np.arange(n), medoids[labels]].sum(dtype=np.float64))
        if best is None or objective < best[2]:
            best = (labels.copy(), medoids.copy(), objective, iteration+1)
    return best
