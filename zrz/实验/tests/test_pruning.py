"""Window-range rewrite, exact cell skipping, and whole-pair skipping."""
import csv
import json
import random
import re
import subprocess
import sys
import tempfile
import unittest
from math import inf, isfinite, sqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import dtw
from dtw import (Config, Result, align, envelope, lb_keogh, _normalize_times, _values, _times,
                 _speed, MatrixDistance, pruning_config, read_matrix, prepare_expression,
                 preprocessing_config)
from diagnose import build_report


# The align implementation before the window-range rewrite, kept verbatim as the reference.
if True:
    def reference_align(x, y, config=None, *, times_x=None, times_y=None, return_path=True):
        """Minimize sum of local costs plus transition penalties, then take sqrt.

        Time coordinates are normalized using the shared time span. Supply both
        coordinate arrays or neither. No post-hoc path-length normalization.
        Ties prefer diagonal, horizontal, vertical states in that order.
        Cost-only mode uses O(len(y)) working memory; traceback uses O(m*n).
        """
        cfg = config or Config()
        x, y = _values(x), _values(y)
        if (times_x is None) != (times_y is None):
            raise ValueError('Supply both time arrays or neither')
        tx, ty = _times(times_x, len(x)), _times(times_y, len(y))
        start, end = min(tx[0], ty[0]), max(tx[-1], ty[-1])
        span = end - start or 1.0
        tx, ty = [(v-start)/span for v in tx], [(v-start)/span for v in ty]
        vx, vy = _speed(x, tx), _speed(y, ty)
        # Each non-diagonal move pays step_penalty plus its run opening/extension.
        a = cfg.step_penalty + cfg.open_penalty
        b = cfg.step_penalty + cfg.extend_penalty
        n = len(y)
        previous = [[inf] * 3 for _ in range(n + 1)]
        previous[0][0] = 0.0
        pointers = {} if return_path else None
        for i in range(1, len(x) + 1):
            current = [[inf] * 3 for _ in range(n + 1)]
            for j in range(1, n + 1):
                delta = abs(tx[i-1] - ty[j-1])
                if cfg.window > 0 and delta > cfg.window + 1e-12:
                    continue
                local = (x[i-1]-y[j-1])**2 + cfg.time_penalty * delta**cfg.power * (1 + cfg.dynamics * (vx[i-1]+vy[j-1])/2)
                candidates = [previous[j-1],
                              [current[j-1][s] + (b if s == 1 else a) for s in range(3)],
                              [previous[j][s] + (b if s == 2 else a) for s in range(3)]]
                for state, options in enumerate(candidates):
                    source = min(range(3), key=options.__getitem__)
                    current[j][state] = local + options[source]
                    if return_path:
                        pointers[i, j, state] = source
            previous = current
        state = min(range(3), key=previous[n].__getitem__)
        cost = previous[n][state]
        if not isfinite(cost):
            raise ValueError('No finite endpoint-to-endpoint alignment; check window or numeric scale')
        if not return_path:
            return Result(cost, sqrt(cost), None, None)
        i, j, path = len(x), n, []
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

def random_config(rng):
    return Config(window=rng.choice([0, 0, 0.0834, 0.1667, 0.25, 0.5, 1]),
                  step_penalty=rng.choice([0, 0, 0.1, 0.5]),
                  time_penalty=rng.choice([0, 0, 0.5, 2, 10]),
                  power=rng.choice([1, 2, 3]),
                  dynamics=rng.choice([0, 0, 0.3]),
                  open_penalty=rng.choice([0, 0, 0.2, 2]),
                  extend_penalty=rng.choice([0, 0, 0.05, 0.5]))


def random_pair(rng, same_grid=True):
    m = rng.randint(2, 13)
    n = m if same_grid else rng.randint(2, 13)
    x = [rng.gauss(0, 1) for _ in range(m)]
    y = [rng.gauss(0, 1) for _ in range(n)]
    if same_grid:
        t = sorted(rng.sample(range(0, 200), m))
        return x, y, t, t
    return x, y, sorted(rng.sample(range(0, 200), m)), sorted(rng.sample(range(0, 200), n))


def outcome(f):
    try:
        r = f()
        return ('ok', r.cost, r.distance, r.path, r.diagnostics)
    except ValueError:
        return ('error',)


class WindowRewriteTests(unittest.TestCase):
    def test_identical_to_previous_implementation(self):
        rng = random.Random(11)
        for k in range(1500):
            cfg = random_config(rng)
            x, y, tx, ty = random_pair(rng, same_grid=k % 3 != 0)
            for path in (True, False):
                new = outcome(lambda: align(x, y, cfg, times_x=tx, times_y=ty, return_path=path))
                old = outcome(lambda: reference_align(x, y, cfg, times_x=tx, times_y=ty, return_path=path))
                self.assertEqual(new, old)


class CellSkippingTests(unittest.TestCase):
    def test_cost_and_path_unchanged(self):
        rng = random.Random(12)
        for k in range(1500):
            cfg = random_config(rng)
            x, y, tx, ty = random_pair(rng, same_grid=k % 4 != 0)
            for path in (True, False):
                plain = outcome(lambda: align(x, y, cfg, times_x=tx, times_y=ty, return_path=path))
                pruned = outcome(lambda: align(x, y, cfg, times_x=tx, times_y=ty, return_path=path,
                                               prune_cells=True))
                self.assertEqual(pruned, plain)

    def test_lower_bound_never_exceeds_dtw(self):
        rng = random.Random(13)
        for _ in range(1500):
            cfg = random_config(rng)
            x, y, t, _ = random_pair(rng)
            tn, _ = _normalize_times([float(v) for v in t], [float(v) for v in t])
            d = align(x, y, cfg, times_x=t, times_y=t, return_path=False).distance
            self.assertLessEqual(lb_keogh(x, envelope(y, tn, cfg.window)), d + 1e-12)
            self.assertLessEqual(lb_keogh(y, envelope(x, tn, cfg.window)), d + 1e-12)


class EarlyAbandonTests(unittest.TestCase):
    def test_abandons_only_above_max_cost(self):
        rng = random.Random(14)
        for _ in range(1500):
            cfg = random_config(rng)
            x, y, t, _ = random_pair(rng)
            exact = align(x, y, cfg, times_x=t, times_y=t, return_path=False)
            for factor in (0.3, 0.9, 1.0, 1.1, 3):
                cap = exact.cost * factor
                r = align(x, y, cfg, times_x=t, times_y=t, return_path=False, max_cost=cap,
                          prune_cells=rng.random() < .5)
                if exact.cost <= cap:
                    self.assertEqual(r.cost, exact.cost)
                elif r.diagnostics is None:
                    self.assertLessEqual(r.cost, cap * (1 + 1e-8) + 1e-12)
                else:
                    self.assertTrue(r.diagnostics['abandoned'])
                    self.assertEqual(r.distance, inf)


class MatrixPruningTests(unittest.TestCase):
    def setUp(self):
        rng = random.Random(15)
        self.times = list(range(0, 130, 10))
        self.values = [[rng.gauss(0, 1) for _ in self.times] for _ in range(25)]
        self.cfg = Config(window=0.1667, open_penalty=0.2, extend_penalty=0.05)

    def full(self):
        n = len(self.values)
        return {(i, j): align(self.values[i], self.values[j], self.cfg, times_x=self.times,
                              times_y=self.times, return_path=False).distance
                for i in range(n) for j in range(i+1, n)}

    def test_pair_skipping_gives_truncated_exact_distances(self):
        full = self.full()
        tau = sorted(full.values())[len(full)//2]
        for cells in (False, True):
            md = MatrixDistance(self.values, self.times, self.cfg, cells=cells)
            md.threshold = tau
            got = {p: md(*p) for p in full}
            for p in full:
                self.assertEqual(got[p], min(full[p], tau))
            self.assertEqual(sum(md.stats.values()), len(full))
            self.assertGreater(md.stats['pairs_skipped_by_lower_bound'] + md.stats['pairs_abandoned_early'], 0)

    def test_quantile_threshold_and_cells_only_mode(self):
        full = self.full()
        md = MatrixDistance(self.values, self.times, self.cfg, cells=True)
        tau, used = md.choose_threshold(0.5, 10**6)
        ordered = sorted(full.values())
        self.assertEqual(used, len(full))
        self.assertEqual(tau, ordered[-(-len(ordered)//2) - 1])
        md2 = MatrixDistance(self.values, self.times, self.cfg, cells=True)
        self.assertEqual({p: md2(*p) for p in full}, full)
        self.assertEqual(md2.cache, {})

    def test_settings_validation(self):
        self.assertEqual(pruning_config(), {'cells': False, 'pairs': False})
        self.assertEqual(pruning_config('pairs')['quantile'], 0.9)
        for args in [('cells', 1.0), ('pairs', 1.0, 0.5), ('pairs', -1), ('pairs', None, 0),
                     ('pairs', None, 1.5), ('bogus',)]:
            with self.assertRaises(ValueError):
                pruning_config(*args)


class CommandLineTests(unittest.TestCase):
    def test_prune_options_recorded_replayed_and_diagnosed(self):
        rng = random.Random(16)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); data = root/'data'; (data/'processed').mkdir(parents=True)
            source = data/'processed'/'centered.csv'
            genes = [f'g{i}' for i in range(12)]
            with source.open('w', newline='') as f:
                w = csv.writer(f); w.writerow(['gene_id'] + [f't{t}' for t in range(0, 130, 10)])
                for g in genes:
                    w.writerow([g] + [round(rng.gauss(0, 1), 4) for _ in range(13)])
            (data/'ground_truth.csv').write_text('gene_id,true_cluster\n' +
                ''.join(f'{g},{"a" if i < 6 else "b"}\n' for i, g in enumerate(genes)))

            def run(name, *options):
                out = root/name
                subprocess.run([sys.executable, '-B', str(ROOT/'dtw.py'), '--input', str(source),
                                '--output', str(out), '--window', '0.1667', *options],
                               check=True, capture_output=True, text=True)
                with (out/'distance_matrix.csv').open() as f:
                    rows = [[float(v) for v in r[1:]] for r in list(csv.reader(f))[1:]]
                return json.loads((out/'run.json').read_text()), rows

            base_meta, base = run('plain')
            self.assertEqual(base_meta['pruning'], {'cells': False, 'pairs': False,
                'pairs_skipped_by_lower_bound': 0, 'pairs_abandoned_early': 0, 'pairs_computed_fully': 66})
            meta, cells = run('cells', '--prune', 'cells')
            self.assertEqual(cells, base)
            meta, both = run('both', '--prune', 'both', '--prune-quantile', '0.5')
            tau = meta['pruning']['threshold']
            self.assertEqual(meta['pruning']['sample_pairs'], 66)
            self.assertIn('pruning.threshold', meta['distance_definition'])
            self.assertEqual(both, [[min(v, tau) for v in row] for row in base])
            replay_meta, replay = run('replay', '--config-file', str(root/'both'/'run.json'))
            self.assertEqual(replay_meta['pruning']['threshold'], tau)
            self.assertEqual(replay, both)
            report = build_report(root/'both', 2)
            payload = json.loads(re.search(r'<script id="report-data" type="application/json">(.*?)</script>',
                                           report.read_text(), re.S).group(1))
            self.assertEqual(payload['warning'], '')
            self.assertEqual(payload['pruning']['threshold'], tau)


if __name__ == '__main__':
    unittest.main()
