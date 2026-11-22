#!/usr/bin/env python3
"""DTW for gene-expression time series: core algorithm plus command-line runner.

Exact, endpoint-constrained DTW with composable penalties (Python 3.10+).
No third-party packages. Horizontal/vertical states repeat a matched sample;
they are not sequence-deletion gaps. All penalties are nonnegative.

Run experiments with explicitly specified parameters; no method presets:
    python3 dtw.py --input centered.csv --output outputs/name [options]
"""
import argparse
import csv
import hashlib
import json
import random
import re
import time
from bisect import bisect_left, bisect_right
from dataclasses import asdict, dataclass
from html import escape
from math import ceil, fsum, inf, isfinite, sqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parent


# ---------------------------------------------------------------------------
# Core algorithm
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Config:
    window: float | None = 0.0  # 0/None disables the elapsed-time constraint
    step_penalty: float = 0.0
    time_penalty: float = 0.0
    power: float = 2.0
    dynamics: float = 0.0
    open_penalty: float | None = 0.0
    extend_penalty: float | None = 0.0

    def __post_init__(self):
        # Accept nulls in old experiment records; all switches now use numeric zero.
        for key in ('window', 'open_penalty', 'extend_penalty'):
            if getattr(self, key) is None:
                object.__setattr__(self, key, 0.0)
        for key in ('step_penalty', 'time_penalty', 'dynamics', 'open_penalty', 'extend_penalty'):
            value = getattr(self, key)
            if not isfinite(value) or value < 0:
                raise ValueError(f'{key} must be finite and nonnegative')
        if not isfinite(self.power) or self.power <= 0:
            raise ValueError('power is an exponent and must be finite and positive')
        if not isfinite(self.window) or not 0 <= self.window <= 1:
            raise ValueError('window must be in [0, 1]; zero disables the constraint')



@dataclass
class Result:
    cost: float
    distance: float
    path: list[tuple[int, int]] | None
    diagnostics: dict | None


def _values(values):
    values = [float(v) for v in values]
    if not values or not all(isfinite(v) for v in values):
        raise ValueError('Sequences must be nonempty and finite; impute missing values first')
    return values


def _times(times, n):
    if times is None:
        return [i / max(1, n - 1) for i in range(n)]
    t = _values(times)
    if len(t) != n or any(b <= a for a, b in zip(t, t[1:])):
        raise ValueError('Time coordinates must match sequence length and strictly increase')
    return t


def _speed(x, t):
    if len(x) == 1:
        return [0.0]
    result = []
    for i in range(len(x)):
        lo, hi = max(0, i - 1), min(len(x) - 1, i + 1)
        result.append(abs((x[hi] - x[lo]) / (t[hi] - t[lo])))
    return result


def _normalize_times(tx, ty):
    """Map both time axes onto [0, 1] using their shared span (as align does)."""
    start, end = min(tx[0], ty[0]), max(tx[-1], ty[-1])
    span = end - start or 1.0
    return [(v-start)/span for v in tx], [(v-start)/span for v in ty]


def _window_columns(ty, t, window):
    """0-based half-open range of columns that can lie inside the window around t.

    Slightly wider than the window; align still applies its exact window test to
    every visited cell, so the result is unchanged and only out-of-window cells
    are no longer visited one by one.
    """
    if window <= 0:
        return 0, len(ty)
    margin = window + 1e-9
    return bisect_left(ty, t - margin), bisect_right(ty, t + margin)


def envelope(y, ty, window, tx=None):
    """LB_Keogh envelope: for each row time, max and min of y over window columns.

    ty (and tx, default ty) must be normalized like align does. With window 0
    every column is allowed, so the envelope is the curve's overall max/min.
    """
    upper, lower = [], []
    for t in (ty if tx is None else tx):
        lo, hi = _window_columns(ty, t, window)
        cols = [y[j] for j in range(lo, hi) if not (window > 0 and abs(t-ty[j]) > window + 1e-12)]
        upper.append(max(cols) if cols else inf)
        lower.append(min(cols) if cols else -inf)
    return upper, lower


def lb_keogh_terms(x, env):
    """Per-row LB_Keogh terms: squared distance from x_i to the envelope band of y."""
    upper, lower = env
    return [(v-u)**2 if v > u else ((v-l)**2 if v < l else 0.0) for v, u, l in zip(x, upper, lower)]


def lb_keogh(x, env):
    """LB_Keogh lower bound on the DTW distance (any penalties only add cost)."""
    return sqrt(fsum(lb_keogh_terms(x, env)))


def align(x, y, config=None, *, times_x=None, times_y=None, return_path=True,
          prune_cells=False, max_cost=None, envelope_y=None):
    """Minimize sum of local costs plus transition penalties, then take sqrt.

    Time coordinates are normalized using the shared time span. Supply both
    coordinate arrays or neither. No post-hoc path-length normalization.
    Ties prefer diagonal, horizontal, vertical states in that order.
    Cost-only mode uses O(len(y)) working memory; traceback uses O(m*n).
    Only columns inside the window are visited.

    Optional speed-ups (both off by default):
    prune_cells -- exact cell skipping. The diagonal path is always a valid
        alignment, so its cost is an upper bound on the optimum. A cell whose
        cumulative cost plus the LB_Keogh bound for the remaining rows exceeds
        it cannot be on the optimal path and is dropped. Cost and path are
        identical with or without it. Needs equal time axes; otherwise ignored.
    max_cost -- early abandoning. Cells that can only lead to a total cost above
        max_cost are dropped; when a whole row is gone, no alignment can cost
        <= max_cost and Result(inf, inf, None, {'abandoned': True}) is returned.
    envelope_y -- precomputed envelope(y, ...) on the same time grid, to avoid
        recomputing it for every pair.
    """
    cfg = config or Config()
    x, y = _values(x), _values(y)
    if (times_x is None) != (times_y is None):
        raise ValueError('Supply both time arrays or neither')
    tx, ty = _times(times_x, len(x)), _times(times_y, len(y))
    tx, ty = _normalize_times(tx, ty)
    vx, vy = _speed(x, tx), _speed(y, ty)
    # Each non-diagonal move pays step_penalty plus its run opening/extension.
    a = cfg.step_penalty + cfg.open_penalty
    b = cfg.step_penalty + cfg.extend_penalty
    m, n = len(x), len(y)
    same_grid = m == n and tx == ty
    limit = inf
    if prune_cells and same_grid:
        # Diagonal path: zero time difference, no stretch moves, inside any window.
        limit = fsum((x[i]-y[i])**2 for i in range(m))
    if max_cost is not None:
        limit = min(limit, max_cost)
    suffix = [0.0] * (m + 1)  # suffix[i]: lower bound for rows i+1..m (1-based)
    if limit < inf:
        env = envelope_y if (envelope_y is not None and same_grid) else envelope(y, ty, cfg.window, tx)
        terms = lb_keogh_terms(x, env)
        for i in range(m - 1, -1, -1):
            suffix[i] = suffix[i+1] + terms[i]
        keep = limit + 1e-9 * max(1.0, limit)  # tolerance: never drop an optimal cell
    else:
        keep = inf
    previous = [[inf] * 3 for _ in range(n + 1)]
    previous[0][0] = 0.0
    pointers = {} if return_path else None
    lo = hi = 0  # columns of the previous row that hold a finite cost
    for i in range(1, m + 1):
        current = [[inf] * 3 for _ in range(n + 1)]
        row_keep = keep - suffix[i]
        wlo, whi = _window_columns(ty, tx[i-1], cfg.window)
        new_lo, new_hi = None, -1
        j = max(wlo + 1, lo, 1)
        while j <= whi:
            if j > hi + 1 and min(current[j-1]) == inf:
                break  # nothing further right in this row is reachable
            delta = abs(tx[i-1] - ty[j-1])
            if cfg.window > 0 and delta > cfg.window + 1e-12:
                j += 1
                continue
            local = (x[i-1]-y[j-1])**2 + cfg.time_penalty * delta**cfg.power * (1 + cfg.dynamics * (vx[i-1]+vy[j-1])/2)
            candidates = [previous[j-1],
                          [current[j-1][s] + (b if s == 1 else a) for s in range(3)],
                          [previous[j][s] + (b if s == 2 else a) for s in range(3)]]
            alive = False
            for state, options in enumerate(candidates):
                source = min(range(3), key=options.__getitem__)
                value = local + options[source]
                if value <= row_keep:
                    current[j][state] = value
                    alive = alive or value < inf
                if return_path:
                    pointers[i, j, state] = source
            if alive:
                if new_lo is None:
                    new_lo = j
                new_hi = j
            j += 1
        if new_lo is None:
            if max_cost is not None:
                return Result(inf, inf, None, {'abandoned': True})
            raise ValueError('No finite endpoint-to-endpoint alignment; check window or numeric scale')
        lo, hi = new_lo, new_hi
        previous = current
    state = min(range(3), key=previous[n].__getitem__)
    cost = previous[n][state]
    if not isfinite(cost):
        if max_cost is not None:
            return Result(inf, inf, None, {'abandoned': True})
        raise ValueError('No finite endpoint-to-endpoint alignment; check window or numeric scale')
    if not return_path:
        return Result(cost, sqrt(cost), None, None)
    i, j, path = m, n, []
    while i and j:
        path.append((i-1, j-1))
        source = pointers[i, j, state]
        if state == 0:
            i, j = i-1, j-1
        elif state == 1:
            j -= 1
        else:
            i -= 1
        state = source
    path.reverse()
    moves = [(q[0]-p[0], q[1]-p[1]) for p, q in zip(path, path[1:])]
    runs, length, last = [], 0, None
    for move in moves:
        if move != (1, 1):
            if move == last:
                length += 1
            else:
                if length:
                    runs.append(length)
                length = 1
        else:
            if length:
                runs.append(length)
            length = 0
        last = move
    if length:
        runs.append(length)
    diagnostics = {'path_length': len(path), 'non_diagonal_steps': sum(runs),
                   'stretch_runs': len(runs), 'stretch_run_lengths': runs,
                   'max_time_displacement': max(abs(tx[i]-ty[j]) for i, j in path)}
    return Result(cost, sqrt(cost), path, diagnostics)


# ---------------------------------------------------------------------------
# Input, preprocessing and command-line runner
# ---------------------------------------------------------------------------

def read_matrix(path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        reader = csv.reader(f)
        header = next(reader, [])
        if len(header) < 2 or header[0] != 'gene_id':
            raise ValueError('Expected gene_id,t0,t10,... header')
        if not all(re.fullmatch(r't[-+]?\d+(?:\.\d+)?', h) for h in header[1:]):
            raise ValueError('Invalid time-column names')
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
                raise ValueError(f'Line {line}: nonnumeric/missing data; use processed/centered.csv') from exc
            if not all(isfinite(v) for v in x):
                raise ValueError(f'Line {line}: nonfinite data')
            ids.append(row[0]); values.append(x); seen.add(row[0])
    if not ids:
        raise ValueError('No genes in input')
    return ids, times, values


def preprocessing_config(scale_floor):
    """A positive floor enables guarded, per-gene amplitude normalization."""
    if not isfinite(scale_floor) or scale_floor < 0:
        raise ValueError('scale_floor must be finite and nonnegative; 0 disables normalization')
    if scale_floor == 0:
        return {'method': 'none', 'version': 1}
    return {'method': 'center_scale', 'scale_floor': scale_floor, 'version': 1}


def prepare_expression(values, preprocessing):
    """Replay the saved transform before both DTW and diagnostic plotting.

    No labels or other genes are used. This makes pair, subset, and full-matrix
    runs agree. Missing preprocessing metadata means legacy unchanged input.
    """
    spec = preprocessing or {'method': 'none', 'version': 1}
    if spec.get('version', 1) != 1:
        raise ValueError('Unsupported expression preprocessing version')
    if spec.get('method') == 'none':
        return [list(row) for row in values]
    if spec.get('method') != 'center_scale':
        raise ValueError('Unknown expression preprocessing method')
    floor = spec['scale_floor']
    if not isfinite(floor) or floor <= 0:
        raise ValueError('Saved center_scale preprocessing requires a positive scale_floor')
    result = []
    for row in values:
        mean = fsum(row) / len(row)
        centered = [v - mean for v in row]
        sd = sqrt(fsum(v*v for v in centered) / len(centered))
        result.append([v / max(sd, floor) for v in centered])
    return result


def plot_pair(path, ids, times, x, y, result):
    """Standalone SVG: original time positions and the exact warping links."""
    lo, hi = min(x+y), max(x+y)
    span = hi-lo or 1
    tspan = times[-1]-times[0] or 1
    px = lambda i: 65+850*(times[i]-times[0])/tspan
    py = lambda v, base: base-105*(v-lo)/span
    parts = ['<svg xmlns="http://www.w3.org/2000/svg" width="980" height="470" viewBox="0 0 980 470">',
             '<rect width="980" height="470" fill="#faf9f6"/>',
             '<g font-family="sans-serif" font-size="14" fill="#202c35">',
             '<text x="45" y="32" font-size="22">DTW alignment</text>']
    for i, j in result.path:
        parts.append(f'<line x1="{px(i)}" y1="{py(x[i],190)}" x2="{px(j)}" y2="{py(y[j],375)}" stroke="#c5cbd0"/>')
    for series, base, color, label in [(x,190,'#007f82',ids[0]),(y,375,'#c35a35',ids[1])]:
        points = ' '.join(f'{px(i)},{py(v,base)}' for i,v in enumerate(series))
        parts.append(f'<polyline points="{points}" fill="none" stroke="{color}" stroke-width="3"/>')
        parts.append(f'<text x="45" y="{base-120}" fill="{color}">{escape(label)}</text>')
        for i,v in enumerate(series):
            parts.append(f'<circle cx="{px(i)}" cy="{py(v,base)}" r="3" fill="{color}"/>')
    for i,t in enumerate(times):
        parts.append(f'<text x="{px(i)}" y="410" text-anchor="middle">{t:g}</text>')
    parts.append(f'<text x="45" y="448">Time (input coordinates) · distance = {result.distance:.6g} · max normalized lag = {result.diagnostics["max_time_displacement"]:.3g}</text></g></svg>')
    path.write_text('\n'.join(parts), encoding='utf-8')


def write_distance_matrix(path, ids, values, distance):
    """Write the full symmetric gene x gene matrix; shared by DTW and the baselines.

    distance(i, j) is called once per unordered pair (gene positions, i < j). Rows are streamed to disk:
    O(N) matrix-output memory. The lower triangle is kept until its row is
    written, using a temporary disk store.
    """
    import sqlite3
    import tempfile
    with tempfile.TemporaryDirectory(prefix='dtw-') as temp:
        db=sqlite3.connect(str(Path(temp)/'distances.db'))
        db.execute('CREATE TABLE d (i INTEGER,j INTEGER,v REAL,PRIMARY KEY(i,j))')
        try:
            with path.open('w',newline='') as f:
                w=csv.writer(f); w.writerow(['gene_id']+ids)
                for i in range(len(ids)):
                    lower=dict(db.execute('SELECT j,v FROM d WHERE i=?',(i,)))
                    row=[lower[j] for j in range(i)]+[0.0]
                    db.execute('DELETE FROM d WHERE i=?',(i,))
                    for j in range(i+1,len(ids)):
                        d=distance(i,j)
                        row.append(d); db.execute('INSERT INTO d VALUES (?,?,?)',(j,i,d))
                    w.writerow([ids[i]]+row)
                    db.commit()
                    if i % 100 == 0:
                        print(f'Rows: {i+1}/{len(ids)}',flush=True)
        finally:
            db.close()


def pruning_config(mode='none', threshold=None, quantile=None, sample=2000):
    """Validate the optional speed-up settings for full-matrix runs.

    cells: exact cell skipping inside each pair's DP table (distances unchanged).
    pairs: whole-pair skipping with LB_Keogh + early abandoning. Distances above
        the threshold are recorded as the threshold: min(DTW distance, threshold).
        The threshold is a distance, or the given quantile of exact DTW distances
        over a fixed-seed random sample of pairs (default quantile 0.9).
    """
    if mode not in ('none', 'cells', 'pairs', 'both'):
        raise ValueError('prune must be none, cells, pairs or both')
    pairs = mode in ('pairs', 'both')
    if not pairs and (threshold is not None or quantile is not None):
        raise ValueError('--prune-threshold/--prune-quantile need --prune pairs or both')
    if threshold is not None and quantile is not None:
        raise ValueError('Give either --prune-threshold or --prune-quantile, not both')
    if threshold is not None and (not isfinite(threshold) or threshold <= 0):
        raise ValueError('prune threshold must be a positive distance')
    if pairs and threshold is None:
        quantile = 0.9 if quantile is None else quantile
        if not 0 < quantile <= 1:
            raise ValueError('prune quantile must be in (0, 1]')
        if sample < 1:
            raise ValueError('prune sample must be positive')
    spec = {'cells': mode in ('cells', 'both'), 'pairs': pairs}
    if pairs:
        spec.update(threshold=threshold, quantile=quantile,
                    sample_pairs=sample if threshold is None else None, sample_seed=0)
    return spec


class MatrixDistance:
    """DTW distance between genes i and j for the full matrix, with optional pruning."""

    def __init__(self, values, times, cfg, cells=False):
        self.values, self.times, self.cfg, self.cells = values, times, cfg, cells
        self.threshold = None
        tn, _ = _normalize_times(list(map(float, times)), list(map(float, times)))
        self.envelopes = [envelope(v, tn, cfg.window) for v in values]
        self.cache = {}
        self.stats = {'pairs_skipped_by_lower_bound': 0, 'pairs_abandoned_early': 0,
                      'pairs_computed_fully': 0}

    def exact(self, i, j):
        if (i, j) not in self.cache:
            self.cache[i, j] = align(self.values[i], self.values[j], self.cfg,
                                     times_x=self.times, times_y=self.times, return_path=False,
                                     prune_cells=self.cells, envelope_y=self.envelopes[j]).distance
        return self.cache[i, j]

    def choose_threshold(self, quantile, sample, seed=0):
        """Quantile of exact distances over a fixed-seed random sample of pairs."""
        n = len(self.values)
        if n < 2:
            raise ValueError('Need at least two genes')
        if n*(n-1)//2 <= sample:
            pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
        else:
            rng, chosen = random.Random(seed), set()
            while len(chosen) < sample:
                chosen.add(tuple(sorted(rng.sample(range(n), 2))))
            pairs = sorted(chosen)
        distances = sorted(self.exact(i, j) for i, j in pairs)
        return distances[max(1, ceil(quantile*len(distances))) - 1], len(pairs)

    def __call__(self, i, j):
        tau = self.threshold
        if (i, j) in self.cache:  # already computed exactly while choosing the threshold
            self.stats['pairs_computed_fully'] += 1
            d = self.cache.pop((i, j))
            return d if tau is None else min(d, tau)
        if tau is None:
            self.stats['pairs_computed_fully'] += 1
            d = self.exact(i, j)
            del self.cache[i, j]  # keep memory flat over the whole matrix
            return d
        x, y = self.values[i], self.values[j]
        # DTW here is symmetric, so both one-sided LB_Keogh bounds are valid; use the larger.
        if max(lb_keogh(x, self.envelopes[j]), lb_keogh(y, self.envelopes[i])) > tau:
            self.stats['pairs_skipped_by_lower_bound'] += 1
            return tau
        result = align(x, y, self.cfg, times_x=self.times, times_y=self.times, return_path=False,
                       prune_cells=self.cells, max_cost=tau*tau, envelope_y=self.envelopes[j])
        if result.diagnostics and result.diagnostics.get('abandoned'):
            self.stats['pairs_abandoned_early'] += 1
            return tau
        self.stats['pairs_computed_fully'] += 1
        return min(result.distance, tau)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--scale-floor', type=float, default=0.1,
                        help='Per-gene SD floor for centering/scaling (default 0.1); 0 uses input unchanged')
    parser.add_argument('--clustering', choices=['kmedoids','average','xqf'], default='xqf',
                        help='Preferred downstream clustering for diagnose; does not change DTW distances. '
                             'xqf = Xu Qiufan k-medoids from xqf/xqf_clustering')
    parser.add_argument('--clustering-seed', type=int, default=42,
                        help='Random seed saved for --clustering xqf (default 42)')
    parser.add_argument('--window', type=float, default=0.0, help='Normalized time window; 0 disables constraint')
    parser.add_argument('--step-penalty', type=float, default=0.0)
    parser.add_argument('--time-penalty', type=float, default=0.0)
    parser.add_argument('--power', type=float, default=2.0)
    parser.add_argument('--dynamics', type=float, default=0.0)
    parser.add_argument('--open-penalty', type=float, default=0.0)
    parser.add_argument('--extend-penalty', type=float, default=0.0)
    parser.add_argument('--pair', nargs=2, metavar=('GENE_A','GENE_B'), help='Only align this pair of gene IDs')
    parser.add_argument('--limit', type=int, help='Explicitly use first N genes (demonstration only)')
    parser.add_argument('--config-file', type=Path, help='Config JSON or run.json; run.json also reuses saved preprocessing, clustering and pruning')
    parser.add_argument('--prune', choices=['none','cells','pairs','both'], default='none',
                        help='Optional speed-ups for the full matrix (default none). cells: exact cell skipping, '
                             'distances unchanged. pairs: LB_Keogh + early abandoning, distances above the '
                             'threshold recorded as the threshold. both: cells and pairs')
    parser.add_argument('--prune-threshold', type=float,
                        help='Distance threshold for --prune pairs/both')
    parser.add_argument('--prune-quantile', type=float,
                        help='Threshold = this quantile of exact distances on a random sample of pairs '
                             '(default 0.9 when no threshold is given)')
    parser.add_argument('--prune-sample', type=int, default=2000,
                        help='Number of random pairs used for --prune-quantile (default 2000, fixed seed)')
    args = parser.parse_args(argv)
    try:
        ids, times, values = read_matrix(args.input)
        original_count = len(ids)
        if args.limit is not None:
            if args.limit < 1:
                raise ValueError('limit must be positive')
            ids, values = ids[:args.limit], values[:args.limit]
        saved = json.loads(args.config_file.read_text(encoding='utf-8')) if args.config_file else None
        if saved is not None and not isinstance(saved, dict):
            raise ValueError('Config file must contain a JSON object')
        cfg = Config(**saved.get('config', saved)) if saved is not None else Config(
            window=args.window, step_penalty=args.step_penalty,
            time_penalty=args.time_penalty, power=args.power, dynamics=args.dynamics,
            open_penalty=args.open_penalty, extend_penalty=args.extend_penalty)
        # Old run.json files explicitly reproduce their old, unscaled inputs.
        saved_run = saved is not None and 'config' in saved
        preprocessing = (saved.get('preprocessing', {'method': 'none', 'version': 1})
                         if saved_run else preprocessing_config(args.scale_floor))
        clustering = saved.get('clustering', 'average') if saved_run else args.clustering
        if clustering not in ('average', 'kmedoids', 'xqf'):
            raise ValueError('Unknown saved clustering algorithm')
        clustering_seed = (saved.get('clustering_seed', 42) if saved_run else args.clustering_seed)
        if saved_run and 'pruning' in saved:
            # Replay the saved settings exactly, including the threshold value actually used.
            old = saved['pruning']
            mode = {(False, False): 'none', (True, False): 'cells',
                    (False, True): 'pairs', (True, True): 'both'}[bool(old.get('cells')), bool(old.get('pairs'))]
            pruning = pruning_config(mode, threshold=old.get('threshold') if old.get('pairs') else None)
        else:
            pruning = pruning_config(args.prune, args.prune_threshold, args.prune_quantile, args.prune_sample)
        values = prepare_expression(values, preprocessing)
        if args.pair and any(g not in ids for g in args.pair):
            raise ValueError('Pair gene ID not present in selected input')
        args.output.mkdir(parents=True, exist_ok=True)
        started = time.perf_counter()
        metadata = {'config':asdict(cfg), 'input':str(args.input.resolve()),
                    'input_sha256':hashlib.sha256(args.input.read_bytes()).hexdigest(),
                    'original_genes':original_count, 'selected_genes':len(ids), 'timepoints':times,
                    'distance_definition':'sqrt(minimum cumulative local and transition cost)',
                    'kernel_sha256':hashlib.sha256((ROOT/'dtw.py').read_bytes()).hexdigest(),
                    'preprocessing':preprocessing, 'clustering':clustering,
                    **({'clustering_seed':clustering_seed} if clustering == 'xqf' else {}),
                    'expression_transform':('centered per gene, divided by max(population SD, scale_floor)'
                                            if preprocessing['method']=='center_scale' else 'none; input values used unchanged'),
                    'labels_read':False, 'mode':'pair' if args.pair else 'matrix'}
        if args.pair:
            a, b = [ids.index(g) for g in args.pair]
            result = align(values[a],values[b],cfg,times_x=times,times_y=times)
            with (args.output/'alignment.csv').open('w',newline='') as f:
                w=csv.writer(f); w.writerow(['index_a','index_b','time_a','time_b','value_a','value_b'])
                for i,j in result.path:
                    w.writerow([i,j,times[i],times[j],values[a][i],values[b][j]])
            plot_pair(args.output/'alignment.svg',args.pair,times,values[a],values[b],result)
            metadata.update(pair=args.pair, cost=result.cost, distance=result.distance, diagnostics=result.diagnostics)
        else:
            distance = MatrixDistance(values, times, cfg, cells=pruning['cells'])
            if pruning['pairs']:
                if pruning['threshold'] is None:
                    pruning['threshold'], pruning['sample_pairs'] = distance.choose_threshold(
                        pruning['quantile'], pruning['sample_pairs'], pruning['sample_seed'])
                distance.threshold = pruning['threshold']
                metadata['distance_definition'] = 'min(sqrt(minimum cumulative local and transition cost), pruning.threshold)'
            write_distance_matrix(args.output/'distance_matrix.csv', ids, values, distance)
            pruning.update(distance.stats)
            metadata['pruning'] = pruning
            metadata['pairs_computed']=len(ids)*(len(ids)-1)//2
        metadata['elapsed_seconds']=time.perf_counter()-started
        (args.output/'run.json').write_text(json.dumps(metadata,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
        print(f'Completed: {args.output.resolve()} ({metadata["elapsed_seconds"]:.3f}s)')
    except (ValueError,OSError,OverflowError,TypeError) as exc:
        parser.error(str(exc))


if __name__ == '__main__':
    main()
