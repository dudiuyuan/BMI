#!/usr/bin/env python3
"""Euclidean-distance baseline: same input and output as dtw.py.

    python3 euclidean.py --input centered.csv --output outputs/name

distance = sqrt( sum over time points of (x_t - y_t)^2 ), comparing equal time
points directly (no time warping). Same scale as DTW: it is DTW restricted to
the diagonal path. Input values are used as they are (no per-gene scaling).
"""
from math import fsum, sqrt

from dtw import run


def euclidean_distance(x, y, times=None):
    return sqrt(fsum((a - b)**2 for a, b in zip(x, y)))


if __name__ == '__main__':
    run('euclidean', euclidean_distance, 'sqrt(sum of squared differences at equal time points)',
        description=__doc__)
