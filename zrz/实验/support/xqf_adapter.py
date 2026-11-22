"""Adapter: run Xu Qiufan's k-medoids (xqf/xqf_clustering/k_medoids.py) on a saved DTW matrix.

Her module is loaded from its file in the shared repo rather than copied,
so updates to her code are picked up automatically. Output is converted to
this project's convention: (assignment, clusters), clusters ordered by
cluster index and each a sorted list of gene positions.
"""
import importlib.util
from pathlib import Path

from support.diagnostic_tools import _validate_matrix

XQF_PATH = Path(__file__).resolve().parents[3] / 'xqf' / 'xqf_clustering' / 'k_medoids.py'
DEFAULT_SEED = 42


def _load_xqf():
    if not XQF_PATH.is_file():
        raise ValueError(f'Cannot find xqf k-medoids code at {XQF_PATH}')
    spec = importlib.util.spec_from_file_location('xqf_k_medoids', XQF_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def xqf_k_medoids(matrix, k, seed=DEFAULT_SEED, n_init=5, max_iter=100):
    """Random multi-start alternating k-medoids (xqf). Labels are never read."""
    _validate_matrix(matrix, k)
    n = len(matrix)
    # Her code demands exact symmetry and an exact zero diagonal; values read
    # back from CSV pass our 1e-10 tolerance check, so make them exact here.
    d = [[0.0 if i == j else (matrix[i][j] if i < j else matrix[j][i]) for j in range(n)]
         for i in range(n)]
    labels, medoids = _load_xqf().k_medoids(d, k, seed=seed, n_init=n_init, max_iter=max_iter)
    if len(labels) != n or any(not 0 <= label < k for label in labels):
        unassigned = sum(1 for label in labels if not 0 <= label < k)
        raise ValueError(f'xqf k-medoids left {unassigned} gene(s) without a valid cluster; '
                         'check assign_clusters in xqf/xqf_clustering/k_medoids.py')
    clusters = [[] for _ in range(k)]
    for gene, label in enumerate(labels):
        clusters[label].append(gene)
    return list(labels), clusters
