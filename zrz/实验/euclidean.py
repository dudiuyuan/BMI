#!/usr/bin/env python3
"""Euclidean-distance baseline for gene-expression time series.

Same input, preprocessing and output as dtw.py; only the distance differs:
    distance = sqrt( sum over time points of (x_t - y_t)^2 )
Equal time points are compared directly, with no time warping. This is the
same scale as DTW (square root of summed squared differences): DTW restricted
to the diagonal path gives exactly this value.

    python3 euclidean.py --input centered.csv --output outputs/name
"""
from math import fsum, sqrt

from support.baseline_runner import run_baseline


def euclidean_distance(x, y):
    if len(x) != len(y):
        raise ValueError('Euclidean distance needs sequences of equal length')
    return sqrt(fsum((a-b)**2 for a, b in zip(x, y)))


if __name__ == '__main__':
    run_baseline('euclidean', euclidean_distance,
                 'sqrt(sum of squared differences at equal time points)', __doc__)
