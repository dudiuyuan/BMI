#!/usr/bin/env python3
"""Pearson-correlation baseline for gene-expression time series.

Same input, preprocessing and output as dtw.py; only the distance differs:
    distance = 1 - r,  r = Pearson correlation over equal time points
so 0 = same shape, 1 = uncorrelated, 2 = mirror image. Only shape counts:
baseline and amplitude are ignored, and there is no time warping.
A flat (constant) curve has no defined r: two flat curves get distance 0,
a flat and a non-flat curve get distance 1 (treated as uncorrelated).

    python3 pearson.py --input centered.csv --output outputs/name
"""
from math import fsum, sqrt

from support.baseline_runner import run_baseline


def pearson_distance(x, y):
    if len(x) != len(y):
        raise ValueError('Pearson distance needs sequences of equal length')
    mx, my = fsum(x)/len(x), fsum(y)/len(y)
    dx, dy = [a-mx for a in x], [b-my for b in y]
    sxx, syy = fsum(a*a for a in dx), fsum(b*b for b in dy)
    # 'Flat' allows rounding residue left after centring a constant curve.
    tiny = lambda v, ss: sqrt(ss/len(v)) <= 1e-12*max(1.0, max(abs(e) for e in v))
    flat_x, flat_y = tiny(x, sxx), tiny(y, syy)
    if flat_x or flat_y:
        return 0.0 if flat_x and flat_y else 1.0
    r = fsum(a*b for a, b in zip(dx, dy)) / sqrt(sxx*syy)
    return 1.0 - max(-1.0, min(1.0, r))  # clamp rounding just outside [-1, 1]


if __name__ == '__main__':
    run_baseline('pearson', pearson_distance,
                 '1 - Pearson correlation at equal time points; flat/flat = 0, flat/non-flat = 1', __doc__)
