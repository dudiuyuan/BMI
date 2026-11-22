"""Grid search helpers for tune_dtw.ipynb.

Kept in a module (not in the notebook) so the worker function can run in
parallel processes on macOS, where notebook-defined functions cannot be sent
to worker processes.
"""
import csv
import json
import time
from collections import Counter
from dataclasses import asdict
from itertools import product
from math import inf, log

from dtw import Config, align, _normalize_times, _speed, _times
from euclidean import euclidean_distance
from pearson import pearson_distance
from support.xqf_adapter import xqf_k_medoids

PARAMS = ('window', 'step_penalty', 'open_penalty', 'extend_penalty', 'time_penalty', 'power', 'dynamics')
FIELDS = (['method', 'name'] + list(PARAMS) +
          ['n_params', 'scheme', 'k', 'ari_mean', 'ari_min', 'nmi_mean', 'silhouette_true',
           'early_late_cocluster', 'seconds'])
FLAGS = {'window': '--window', 'step_penalty': '--step-penalty', 'open_penalty': '--open-penalty',
         'extend_penalty': '--extend-penalty', 'time_penalty': '--time-penalty', 'power': '--power',
         'dynamics': '--dynamics'}


# --------------------------------------------------------------------------- metrics

def ari(labels, assignment):
    n = len(labels)
    choose2 = lambda m: m*(m-1)/2
    cell = Counter(zip(labels, assignment))
    observed = sum(map(choose2, cell.values()))
    r = sum(map(choose2, Counter(labels).values())); c = sum(map(choose2, Counter(assignment).values()))
    expected = r*c/choose2(n); maximum = (r+c)/2
    return (observed-expected)/(maximum-expected) if maximum != expected else 1.0


def nmi(labels, assignment):
    """Normalized mutual information, arithmetic-mean normalization (sklearn default)."""
    n = len(labels)
    a, b, ab = Counter(labels), Counter(assignment), Counter(zip(labels, assignment))
    h = lambda counts: -sum(v/n*log(v/n) for v in counts.values())
    ha, hb = h(a), h(b)
    mi = sum(v/n*log((v/n)/((a[x]/n)*(b[y]/n))) for (x, y), v in ab.items())
    return 1.0 if ha == hb == 0 else (mi/((ha+hb)/2) if ha+hb > 0 else 0.0)


def silhouette(matrix, labels):
    """Mean silhouette of the TRUE labels: how well the distances separate the real
    types, independent of any clustering run. Ranges from -1 to 1."""
    groups = {}
    for i, l in enumerate(labels):
        groups.setdefault(l, []).append(i)
    if len(groups) < 2:
        return 0.0
    total = 0.0
    for i, l in enumerate(labels):
        own = groups[l]
        if len(own) == 1:
            continue
        a = sum(matrix[i][j] for j in own if j != i)/(len(own)-1)
        b = min(sum(matrix[i][j] for j in g)/len(g) for m, g in groups.items() if m != l)
        if max(a, b) > 0:
            total += (b-a)/max(a, b)
    return total/len(labels)


# --------------------------------------------------------------------------- grid

def describe(params):
    """Short name listing only the parameters that are switched on."""
    on = [f'{k}={v:g}' for k, v in params.items() if (v and k != 'power') or (k == 'power' and v != 2)]
    return ', '.join(on) if on else 'standard DTW'


def effective_key(params):
    """Settings that give the same distances share one key:
    - power and dynamics do nothing when time_penalty is 0;
    - the DP only sees step+open and step+extend, so e.g. step=0.1 is the same
      as open=0.1, extend=0.1."""
    q = dict(params)
    if not q['time_penalty']:
        q['power'], q['dynamics'] = 2.0, 0.0
    a, b = q['step_penalty'] + q['open_penalty'], q['step_penalty'] + q['extend_penalty']
    return json.dumps([round(q['window'], 12), round(a, 12), round(b, 12), round(q['time_penalty'], 12),
                       round(q['power'], 12), round(q['dynamics'], 12)])


def expand_grid(grid):
    """All combinations of the candidate lists, with equivalent settings merged.
    Returns (unique parameter dicts, number of raw combinations)."""
    missing = set(PARAMS) - set(grid)
    if missing:
        raise ValueError(f'GRID is missing: {sorted(missing)}')
    names = list(PARAMS)
    raw = 0
    unique = {}
    for values in product(*(grid[n] for n in names)):
        raw += 1
        params = {n: float(v) for n, v in zip(names, values)}
        Config(**params)  # validates ranges
        if not params['time_penalty']:
            params['power'], params['dynamics'] = 2.0, 0.0
        key = effective_key(params)
        simpler = lambda q: (sum(1 for k, v in q.items() if k != 'power' and v), q['open_penalty'], q['extend_penalty'])
        if key not in unique or simpler(params) < simpler(unique[key]):
            unique[key] = params
    return list(unique.values()), raw


def command(params, scale_floor=0.1):
    """dtw.py command line for these parameters; --scale-floor is added when it is not
    dtw.py's default (0.1), so the hand-over reproduces the tuning conditions."""
    parts = [f'{FLAGS[k]} {v:g}' for k, v in params.items() if (v and k != 'power') or (k == 'power' and v != 2)]
    if scale_floor != 0.1:
        parts.append(f'--scale-floor {scale_floor:g}')
    return 'python3 dtw.py --input <centered.csv> --output outputs/<name> ' + ' '.join(parts)


# --------------------------------------------------------------------------- fast DTW matrix

def dtw_matrix(values, times, cfg, chunk=50000):
    """Full DTW distance matrix; same numbers as dtw.align, computed faster.

    With numpy, all gene pairs advance through the DP table together: each of
    the T x T cells is one array operation over every pair, instead of a Python
    loop per pair. The recurrence, penalties, window test and operation order
    are those of dtw.align (equal time grid for all genes), so the distances
    are identical. Without numpy it falls back to calling align pair by pair.
    """
    n = len(values)
    matrix = [[0.0]*n for _ in range(n)]
    try:
        import numpy as np
    except ImportError:
        for i in range(n):
            for j in range(i+1, n):
                matrix[i][j] = matrix[j][i] = align(values[i], values[j], cfg, times_x=times, times_y=times,
                                                    return_path=False).distance
        return matrix
    T = len(times)
    t, _ = _normalize_times(_times(times, T), _times(times, T))
    X = np.asarray(values, dtype=float)
    V = np.asarray([_speed([float(v) for v in row], t) for row in values], dtype=float)
    a = cfg.step_penalty + cfg.open_penalty
    b = cfg.step_penalty + cfg.extend_penalty
    # Per cell: time-penalty factor (same scalar expression as align) or None outside the window.
    cells = {}
    for i in range(T):
        for j in range(T):
            delta = abs(t[i] - t[j])
            if not (cfg.window > 0 and delta > cfg.window + 1e-12):
                cells[i, j] = cfg.time_penalty * delta**cfg.power
    rows, cols = np.triu_indices(n, 1)
    for start in range(0, len(rows), chunk):
        ri, ci = rows[start:start+chunk], cols[start:start+chunk]
        x, y, vx, vy = X[ri], X[ci], V[ri], V[ci]
        P = len(ri)
        inf_row = np.full(P, inf)
        previous = [[inf_row]*3 for _ in range(T+1)]
        previous[0] = [np.zeros(P), inf_row, inf_row]
        for i in range(1, T+1):
            current = [[inf_row]*3 for _ in range(T+1)]
            for j in range(1, T+1):
                factor = cells.get((i-1, j-1))
                if factor is None:
                    continue
                local = (x[:, i-1]-y[:, j-1])**2 + factor * (1 + cfg.dynamics * (vx[:, i-1]+vy[:, j-1])/2)
                d, h, v = previous[j-1], current[j-1], previous[j]
                current[j] = [local + np.minimum(np.minimum(d[0], d[1]), d[2]),
                              local + np.minimum(np.minimum(h[0]+a, h[1]+b), h[2]+a),
                              local + np.minimum(np.minimum(v[0]+a, v[1]+a), v[2]+b)]
            previous = current
        last = previous[T]
        dist = np.sqrt(np.minimum(np.minimum(last[0], last[1]), last[2]))
        if not np.all(np.isfinite(dist)):
            raise ValueError('No finite endpoint-to-endpoint alignment; check window or numeric scale')
        for i, j, d in zip(ri.tolist(), ci.tolist(), dist.tolist()):
            matrix[i][j] = matrix[j][i] = d
    return matrix


# --------------------------------------------------------------------------- one configuration

def evaluate(job):
    """Compute one distance matrix and score it under every label scheme."""
    method, params, values, times, truth, schemes, seeds = job
    started = time.perf_counter()
    n = len(values)
    if method == 'dtw':
        matrix = dtw_matrix(values, times, Config(**params))
    else:
        f = euclidean_distance if method == 'euclidean' else pearson_distance
        matrix = [[0.0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                matrix[i][j] = matrix[j][i] = f(values[i], values[j])
    seconds = time.perf_counter() - started
    early = [i for i, l in enumerate(truth['true_cluster']) if l == 'early_peak']
    late = [i for i, l in enumerate(truth['true_cluster']) if l == 'late_peak']
    key = method if method != 'dtw' else effective_key(params)
    name = method if method != 'dtw' else describe(params)
    rows = []
    for scheme, k in schemes.items():
        labels = truth[scheme]
        runs = [xqf_k_medoids(matrix, k, seed=s)[0] for s in seeds]
        aris = [ari(labels, r) for r in runs]
        together = [sum(r[i] == r[j] for i in early for j in late)/(len(early)*len(late)) if early and late else 0.0
                    for r in runs]
        rows.append({'key': key, 'name': name, 'method': method, 'params': params,
                     'n_params': sum(1 for kk, v in params.items() if kk != 'power' and v),
                     'scheme': scheme, 'k': k, 'ari_mean': sum(aris)/len(aris), 'ari_min': min(aris),
                     'nmi_mean': sum(nmi(labels, r) for r in runs)/len(runs),
                     'silhouette_true': silhouette(matrix, labels),
                     'early_late_cocluster': sum(together)/len(together), 'seconds': seconds})
    return rows


def rank_key(row):
    """Higher ARI first (mean, then worst seed), then NMI, then how well the distances
    separate the true types, then fewer parameters switched on, then faster."""
    return (-row['ari_mean'], -row['ari_min'], -row['nmi_mean'], -row['silhouette_true'],
            row['n_params'], row['seconds'])


# --------------------------------------------------------------------------- saving / resuming

def load_results(path):
    """Read results.csv: one column per parameter (empty for Euclidean/Pearson).
    The resume key is rebuilt from those columns."""
    if not path.exists() or path.stat().st_size == 0:
        return []
    rows = []
    with path.open(newline='') as f:
        for r in csv.DictReader(f):
            for key in ('ari_mean', 'ari_min', 'nmi_mean', 'silhouette_true', 'early_late_cocluster', 'seconds'):
                r[key] = float(r[key])
            r['n_params'], r['k'] = int(r['n_params']), int(r['k'])
            r['params'] = {} if r['method'] != 'dtw' else {n: float(r[n]) for n in PARAMS}
            r['key'] = r['method'] if r['method'] != 'dtw' else effective_key(r['params'])
            rows.append(r)
    return rows


def append_results(path, rows):
    new = not path.exists() or path.stat().st_size == 0
    with path.open('a', newline='') as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, extrasaction='ignore')
        if new:
            w.writeheader()
        w.writerows({**r, **{n: r['params'].get(n, '') for n in PARAMS}} for r in rows)


def parameter_effects(results, scheme, grid):
    """For each parameter value: how many configurations used it, their mean ARI and the
    best ARI among them. Shows which values tend to help, beyond the single winner."""
    dtw_rows = [r for r in results if r['scheme'] == scheme and r['method'] == 'dtw']
    out = []
    for name in PARAMS:
        for v in grid[name]:
            rows = [r for r in dtw_rows if r['params'][name] == float(v)]
            if name in ('power', 'dynamics'):
                rows = [r for r in rows if r['params']['time_penalty']]
            if rows:
                out.append({'parameter': name, 'value': v, 'configs': len(rows),
                            'ari_mean': sum(r['ari_mean'] for r in rows)/len(rows),
                            'ari_best': max(r['ari_mean'] for r in rows)})
    return out
