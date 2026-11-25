from math import fsum, sqrt

from dtw import get_args, run


def euclidean_distance(x, y, times):
    squares = []
    for a, b in zip(x, y):
        squares.append((a - b) ** 2)
    return sqrt(fsum(squares))


def main():
    args = get_args()
    run('euclidean', euclidean_distance, args.input, args.output)


if __name__ == '__main__':
    main()
