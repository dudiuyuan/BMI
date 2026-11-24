#!/usr/bin/env python3
"""DTW distances between gene-expression time series, with five fixed modes.

    python3 dtw.py --mode diagonal-time --input centered.csv --output outputs/name
    python3 dtw.py --list

Modes (parameters fixed; chosen by grid search, see README.md):
    standard                 plain DTW, no penalties
    time                     time penalty
    time-dynamics            time penalty + expression-change rate
    diagonal-time-dynamics   diagonal window + time penalty + expression-change rate
    diagonal-time            diagonal window + time penalty

Input values are used as they are (no per-gene scaling). Python 3.10+, no
third-party packages. Writes distance_matrix.csv and run.json to --output.
"""
import argparse
import csv
import hashlib
import json
import re
import time
from array import array
from bisect import bisect_left, bisect_right
from dataclasses import asdict, dataclass
from math import inf, isfinite, sqrt
from pathlib import Path

MODES = {
    'standard': ('标准 DTW，不加任何惩罚', {}),
    'time': ('时间惩罚', {'time_penalty': 2}),
    'time-dynamics': ('时间惩罚 + 表达变化速率', {'extend_penalty': 0.05, 'time_penalty': 0.5, 'dynamics': 0.3}),
    'diagonal-time-dynamics': ('对角线窗口 + 时间惩罚 + 表达变化速率',
                               {'window': 0.25, 'extend_penalty': 0.05, 'time_penalty': 0.5, 'dynamics': 0.3}),
    'diagonal-time': ('对角线窗口 + 时间惩罚', {'window': 0.1667, 'extend_penalty': 0.05, 'time_penalty': 2}),
}


@dataclass(frozen=True)
class Config:
    """DTW parameters; every penalty off (0) unless a mode sets it.

    window          max time difference of matched points, as a fraction of the total time span (0 = no limit)
    step_penalty    extra cost of every horizontal/vertical (stretching) step
    open_penalty    extra cost of the first step of a stretch
    extend_penalty  extra cost of each further step of the same stretch
    time_penalty    weight of |t_i - t_j|**power added to each matched pair
    power           exponent of the time difference
    dynamics        raises the time penalty where expression changes fast
    """
    window: float = 0.0
    step_penalty: float = 0.0
    open_penalty: float = 0.0
    extend_penalty: float = 0.0
    time_penalty: float = 0.0
    power: float = 2.0
    dynamics: float = 0.0


def _speed(x, t):
    """Local rate of change |dx/dt| at each point (central difference, one-sided at the ends)."""
    if len(x) == 1:
        return [0.0]
    out = []
    for i in range(len(x)):
        lo, hi = max(0, i - 1), min(len(x) - 1, i + 1)
        out.append(abs((x[hi] - x[lo]) / (t[hi] - t[lo])))
    return out


def dtw_distance(x, y, times, cfg):
    """DTW distance between two curves sampled at the same times.

    Minimizes the summed local cost along an endpoint-to-endpoint path:
    (x_i - y_j)^2 + time_penalty * |t_i - t_j|^power * (1 + dynamics * (v_i + v_j) / 2),
    plus stretch penalties for horizontal/vertical steps; returns sqrt of the
    minimum. Times are rescaled to [0, 1]. Only cells inside the window are visited.
    """
    span = (times[-1] - times[0]) or 1.0
    t = [(v - times[0]) / span for v in times]
    vx, vy = _speed(x, t), _speed(y, t)
    a = cfg.step_penalty + cfg.open_penalty      # first step of a stretch
    b = cfg.step_penalty + cfg.extend_penalty    # further steps of the same stretch
    n = len(t)
    # Three states per cell: reached diagonally, horizontally, vertically.
    previous = [[inf] * 3 for _ in range(n + 1)]
    previous[0][0] = 0.0
    for i in range(1, n + 1):
        current = [[inf] * 3 for _ in range(n + 1)]
        if cfg.window > 0:
            lo = bisect_left(t, t[i-1] - cfg.window - 1e-9)
            hi = bisect_right(t, t[i-1] + cfg.window + 1e-9)
        else:
            lo, hi = 0, n
        for j in range(max(lo + 1, 1), hi + 1):
            delta = abs(t[i-1] - t[j-1])
            if cfg.window > 0 and delta > cfg.window + 1e-12:
                continue
            local = ((x[i-1] - y[j-1])**2
                     + cfg.time_penalty * delta**cfg.power * (1 + cfg.dynamics * (vx[i-1] + vy[j-1]) / 2))
            d, h, v = previous[j-1], current[j-1], previous[j]
            current[j][0] = local + min(d)
            current[j][1] = local + min(h[0] + a, h[1] + b, h[2] + a)
            current[j][2] = local + min(v[0] + a, v[1] + a, v[2] + b)
        previous = current
    cost = min(previous[n])
    if not isfinite(cost):
        raise ValueError('No alignment fits inside the window')
    return sqrt(cost)


def read_matrix(path):
    """Read gene_id,t0,t10,... ; returns (ids, times, values). No missing values allowed."""
    with open(path, encoding='utf-8-sig', newline='') as f:
        reader = csv.reader(f)
        header = next(reader, [])
        if len(header) < 3 or header[0] != 'gene_id':
            raise ValueError('Expected header gene_id,t0,t10,...')
        if not all(re.fullmatch(r't[-+]?\d+(?:\.\d+)?', h) for h in header[1:]):
            raise ValueError('Time columns must be named like t0, t10, t20')
        times = [float(h[1:]) for h in header[1:]]
        if any(b <= a for a, b in zip(times, times[1:])):
            raise ValueError('Time columns must strictly increase')
        ids, values, seen = [], [], set()
        for line, row in enumerate(reader, 2):
            if len(row) != len(header):
                raise ValueError(f'Line {line}: wrong number of fields')
            if not row[0].strip() or row[0] in seen:
                raise ValueError(f'Line {line}: empty or duplicate gene_id')
            try:
                x = [float(v) for v in row[1:]]
            except ValueError as exc:
                raise ValueError(f'Line {line}: missing or non-numeric value; use processed/centered.csv') from exc
            if not all(isfinite(v) for v in x):
                raise ValueError(f'Line {line}: non-finite value')
            ids.append(row[0]); values.append(x); seen.add(row[0])
    if len(ids) < 2:
        raise ValueError('Need at least two genes')
    return ids, times, values


def run(method, distance, definition, argv=None, description=None, extra=None):
    """Shared command line: read input, compute all pairwise distances, write the matrix and run.json."""
    parser = argparse.ArgumentParser(description=description or __doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    if method == 'dtw':
        parser.add_argument('--mode', choices=list(MODES), help='Which DTW mode (see --list)')
        parser.add_argument('--list', action='store_true', help='Show the modes and their parameters')
    parser.add_argument('--input', type=Path, help='CSV: gene_id,t0,t10,... (e.g. processed/centered.csv)')
    parser.add_argument('--output', type=Path, help='Result folder')
    parser.add_argument('--limit', type=int, help='Only use the first N genes (quick test)')
    args = parser.parse_args(argv)
    if method == 'dtw' and args.list:
        for name, (text, params) in MODES.items():
            shown = '，'.join(f'{k}={v:g}' for k, v in params.items()) or '全部关闭'
            print(f'{name:<24} {text}（{shown}）')
        return
    if method == 'dtw' and not args.mode:
        parser.error('choose --mode (see --list)')
    if not args.input or not args.output:
        parser.error('--input and --output are required')
    try:
        ids, times, values = read_matrix(args.input)
        if args.limit is not None:
            if args.limit < 2:
                raise ValueError('limit must be at least 2')
            ids, values = ids[:args.limit], values[:args.limit]
        f = distance(args) if method == 'dtw' else distance
        args.output.mkdir(parents=True, exist_ok=True)
        started = time.perf_counter()
        n = len(ids)
        d = array('d', bytes(8 * n * n))
        for i in range(n):
            for j in range(i + 1, n):
                d[i*n + j] = d[j*n + i] = f(values[i], values[j], times)
            if i % 500 == 0 and n > 500:
                print(f'Rows: {i+1}/{n}', flush=True)
        with (args.output/'distance_matrix.csv').open('w', newline='') as out:
            w = csv.writer(out)
            w.writerow(['gene_id'] + ids)
            for i in range(n):
                w.writerow([ids[i]] + list(d[i*n:(i+1)*n]))
        meta = {'method': method, **(extra(args) if extra else {}),
                'distance_definition': definition, 'preprocessing': 'none (input values used as they are)',
                'input': str(args.input.resolve()), 'input_sha256': hashlib.sha256(args.input.read_bytes()).hexdigest(),
                'genes': n, 'timepoints': times, 'elapsed_seconds': time.perf_counter() - started}
        (args.output/'run.json').write_text(json.dumps(meta, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
        print(f'Completed: {args.output} ({meta["elapsed_seconds"]:.1f}s)')
    except (ValueError, OSError) as exc:
        parser.error(str(exc))


def main(argv=None):
    run('dtw', lambda args: (lambda x, y, t, cfg=Config(**MODES[args.mode][1]): dtw_distance(x, y, t, cfg)),
        'sqrt(minimum cumulative cost over endpoint-to-endpoint warping paths)', argv,
        extra=lambda args: {'mode': args.mode, 'mode_description': MODES[args.mode][0],
                            'params': asdict(Config(**MODES[args.mode][1]))})


if __name__ == '__main__':
    main()
