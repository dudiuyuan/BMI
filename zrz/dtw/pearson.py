
from math import fsum, sqrt

from dtw import get_args, run


def is_flat(values, sum_squares):
    biggest = 1.0
    for v in values:
        biggest = max(biggest, abs(v))
    return sqrt(sum_squares / len(values)) <= 1e-12 * biggest


def pearson_distance(x, y, times):
    mean_x = fsum(x) / len(x)
    mean_y = fsum(y) / len(y)
    dx = []
    dy = []
    for a in x: dx.append(a - mean_x)
    for b in y: dy.append(b - mean_y)

    sxx = fsum([a * a for a in dx])
    syy = fsum([b * b for b in dy])
    flat_x = is_flat(x, sxx)
    flat_y = is_flat(y, syy)
    if flat_x and flat_y: return 0.0
    if flat_x or flat_y: return 1.0

    sxy = fsum([a * b for a, b in zip(dx, dy)])
    r = sxy / sqrt(sxx * syy)
    r = max(-1.0, min(1.0, r))    
    return 1.0 - r    # 0:same, 1:unrelated, 2:opposite


def main():
    args = get_args()
    run('pearson', pearson_distance, args.input, args.output)


if __name__ == '__main__':
    main()
