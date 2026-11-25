#!/usr/bin/env python3
"""Pearson-correlation baseline: same input and output as dtw.py.

    python3 pearson.py --input centered.csv --output outputs/name

distance = 1 - r (r = Pearson correlation over equal time points):
0 = same shape, 1 = uncorrelated, 2 = mirror image. Only the shape counts,
not baseline or amplitude; no time warping. A flat curve has no defined r:
two flat curves get 0, a flat and a non-flat curve get 1.
"""
from math import fsum, sqrt

from dtw import run


def pearson_distance(x, y, times=None):
    mx, my = fsum(x) / len(x), fsum(y) / len(y)
    dx, dy = [a - mx for a in x], [b - my for b in y]
    sxx, syy = fsum(a * a for a in dx), fsum(b * b for b in dy)
    flat = lambda v, ss: sqrt(ss / len(v)) <= 1e-12 * max(1.0, max(abs(e) for e in v))
    fx, fy = flat(x, sxx), flat(y, syy)
    if fx or fy:
        return 0.0 if fx and fy else 1.0
    r = fsum(a * b for a, b in zip(dx, dy)) / sqrt(sxx * syy)
    return 1.0 - max(-1.0, min(1.0, r))


if __name__ == '__main__':
    run('pearson', pearson_distance, '1 - Pearson correlation; flat/flat = 0, flat/non-flat = 1',
        description=__doc__)
