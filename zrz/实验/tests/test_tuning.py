"""Grid-search tuning: merging of equivalent settings, metrics, and a full notebook run."""
import csv
import json
import random
import sys
import tempfile
import types
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from dtw import Config, align
from support.tuning import (ari, nmi, silhouette, expand_grid, effective_key, command, describe, dtw_matrix,
                            load_results, PARAMS)

NOTEBOOK = ROOT/'tune_dtw.ipynb'


def run_notebook(settings, notebook=None):
    """Execute a notebook's code cells in order, replacing its settings cell (no Jupyter needed)."""
    if 'IPython.display' not in sys.modules:
        display = types.ModuleType('IPython.display')
        display.Markdown, display.display = str, lambda *a, **k: None
        ipython = sys.modules.setdefault('IPython', types.ModuleType('IPython'))
        # matplotlib checks these whenever it draws: report 'no running IPython'.
        ipython.get_ipython = getattr(ipython, 'get_ipython', lambda: None)
        ipython.version_info = getattr(ipython, 'version_info', (0, 0, 0))
        sys.modules['IPython.display'] = display
    nb = notebook or NOTEBOOK
    cells = [''.join(c['source']) for c in json.loads(nb.read_text())['cells'] if c['cell_type'] == 'code']
    assert cells[0].startswith('# 设置'), 'first code cell must be the settings cell'
    ns = {}
    exec(settings, ns)
    for src in cells[1:]:
        exec(compile(src, str(nb), 'exec'), ns)
    return ns


def tiny_dataset(folder):
    rng = random.Random(5)
    shapes = {'early_peak': [0, 1, 3, 1, 0, 0, 0, 0], 'late_peak': [0, 0, 0, 0, 1, 3, 1, 0],
              'oscillatory': [0, 2, 0, -2, 0, 2, 0, -2], 'transient_induction': [0, 0, 2, 2, 2, 0, 0, 0],
              'transient_repression': [0, 0, -2, -2, -2, 0, 0, 0], 'flat_noise': [0]*8}
    (folder/'processed').mkdir(parents=True)
    with (folder/'processed'/'centered.csv').open('w', newline='') as f, \
            (folder/'ground_truth.csv').open('w', newline='') as g:
        w = csv.writer(f); w.writerow(['gene_id'] + [f't{10*i}' for i in range(8)])
        t = csv.writer(g); t.writerow(['gene_id', 'true_cluster', 'shape_cluster'])
        for kind, shape in shapes.items():
            for r in range(2):
                w.writerow([f'{kind}_{r}'] + [v + rng.gauss(0, .05) for v in shape])
                t.writerow([f'{kind}_{r}', kind, 'broad_peak' if 'peak' in kind else kind])


SMALL_GRID = ("GRID = {'window': [0, 0.25], 'step_penalty': [0, 0.5], 'open_penalty': [0, 0.5], "
              "'extend_penalty': [0, 0.5], 'time_penalty': [0, 2], 'power': [1, 2], 'dynamics': [0, 1]}\n")


class GridTests(unittest.TestCase):
    def test_merged_settings_really_give_identical_distances(self):
        grid = {'window': [0, 0.1667], 'step_penalty': [0, 0.1, 0.5], 'open_penalty': [0, 0.1, 0.4],
                'extend_penalty': [0, 0.1, 0.5], 'time_penalty': [0, 2], 'power': [1, 2], 'dynamics': [0, 1]}
        from itertools import product
        rng = random.Random(1)
        x = [rng.gauss(0, 1) for _ in range(9)]; y = [rng.gauss(0, 1) for _ in range(9)]
        by_key = {}
        for combo in product(*grid.values()):
            q = dict(zip(grid, map(float, combo)))
            by_key.setdefault(effective_key(q), set()).add(round(align(x, y, Config(**q)).cost, 12))
        self.assertTrue(all(len(costs) == 1 for costs in by_key.values()))
        unique, raw = expand_grid(grid)
        self.assertEqual(raw, 2*3*3*3*2*2*2)
        self.assertEqual(len(unique), len(by_key))

    def test_merging_rules_and_names(self):
        same = lambda **kw: effective_key({**dict(window=0., step_penalty=0., open_penalty=0., extend_penalty=0.,
                                                  time_penalty=0., power=2., dynamics=0.), **kw})
        self.assertEqual(same(power=1, dynamics=1), same())                     # no time penalty
        self.assertNotEqual(same(time_penalty=2, power=1), same(time_penalty=2))
        self.assertEqual(same(step_penalty=.1), same(open_penalty=.1, extend_penalty=.1))
        unique, _ = expand_grid({'window': [0], 'step_penalty': [0, .1], 'open_penalty': [0, .1],
                                 'extend_penalty': [0, .1], 'time_penalty': [0], 'power': [2], 'dynamics': [0]})
        self.assertIn({'window': 0., 'step_penalty': .1, 'open_penalty': 0., 'extend_penalty': 0.,
                       'time_penalty': 0., 'power': 2., 'dynamics': 0.}, unique)  # simplest form kept
        self.assertEqual(describe(unique[0]), 'standard DTW')
        with self.assertRaises(ValueError):
            expand_grid({'window': [0]})

    def test_metrics_hand_values(self):
        self.assertEqual(ari(list('aabb'), [1, 1, 0, 0]), 1)
        self.assertAlmostEqual(ari(list('aabb'), [0, 1, 0, 1]), -0.5)
        self.assertAlmostEqual(nmi(list('aabb'), [0, 1, 0, 1]), 0)
        self.assertAlmostEqual(nmi([0, 0, 1, 1], [0, 0, 1, 2]), 0.8)  # sklearn reference value
        d = [[abs(a-b) for b in (0, 1, 10, 11)] for a in (0, 1, 10, 11)]
        self.assertAlmostEqual(silhouette(d, list('aabb')), 0.9, places=2)
        self.assertEqual(command({'window': .25, 'step_penalty': 0, 'open_penalty': 0, 'extend_penalty': 0,
                                  'time_penalty': 2, 'power': 3, 'dynamics': 0}).split('<name> ')[1],
                         '--window 0.25 --time-penalty 2 --power 3')
        q = dict(window=0, step_penalty=0, open_penalty=0, extend_penalty=.5, time_penalty=10, power=2, dynamics=0)
        self.assertNotIn('--scale-floor', command(q))
        self.assertTrue(command(q, 0).endswith('--extend-penalty 0.5 --time-penalty 10 --scale-floor 0'))


class FastMatrixTests(unittest.TestCase):
    def test_numpy_matrix_matches_align(self):
        rng = random.Random(21)
        times = sorted(rng.sample(range(0, 200), 11))
        values = [[rng.gauss(0, 1) for _ in times] for _ in range(12)]
        for _ in range(60):
            cfg = Config(window=rng.choice([0, 0.0834, 0.1667, 0.3, 1]), step_penalty=rng.choice([0, .1, .5]),
                         open_penalty=rng.choice([0, .2, 2]), extend_penalty=rng.choice([0, .05, .5]),
                         time_penalty=rng.choice([0, .5, 2, 10]), power=rng.choice([1, 2, 3]),
                         dynamics=rng.choice([0, .3, 1]))
            fast = dtw_matrix(values, times, cfg, chunk=7)   # small chunks exercise the chunk loop
            for i in range(12):
                self.assertEqual(fast[i][i], 0.0)
                for j in range(i+1, 12):
                    slow = align(values[i], values[j], cfg, times_x=times, times_y=times, return_path=False).distance
                    self.assertAlmostEqual(fast[i][j], slow, delta=1e-12)
                    self.assertEqual(fast[i][j], fast[j][i])


class NotebookTests(unittest.TestCase):
    def test_full_run_parallel_resume_and_reports(self):
        with tempfile.TemporaryDirectory() as tmp:
            data, out = Path(tmp)/'d', Path(tmp)/'o'
            tiny_dataset(data)
            settings = (f"DATA={str(data)!r}\nOUTPUT={str(out)!r}\nSCALE_FLOOR=0.1\nN_SEEDS=3\nWORKERS=2\nTOP=5\n"
                        "SCHEMES={'true_cluster': 6, 'shape_cluster': 5}\n" + SMALL_GRID)
            ns = run_notebook(settings)
            with (out/'results.csv').open() as f:
                rows = list(csv.DictReader(f))
            self.assertNotIn('key', rows[0]); self.assertNotIn('params', rows[0])
            self.assertTrue(all(name in rows[0] for name in PARAMS))  # one column per parameter
            keys = {r['key'] for r in load_results(out/'results.csv')}
            self.assertEqual(len(keys), len(ns['configs']) + 2)          # every merged config + 2 references
            self.assertEqual(len(rows), 2*len(keys))                        # one row per label scheme
            self.assertTrue({'euclidean', 'pearson', 'standard DTW'} <= {r['name'] for r in rows})
            rec = json.loads((out/'recommended.json').read_text())
            self.assertEqual(set(rec['schemes']), {'true_cluster', 'shape_cluster'})
            self.assertEqual(len(rec['schemes']['true_cluster']['runners_up']), 5)
            self.assertIn('python3 dtw.py', rec['schemes']['true_cluster']['best']['command'])
            self.assertIn('平均 ARI', (out/'summary.md').read_text())
            again = run_notebook(settings)                                   # resume: nothing left to do
            self.assertEqual(again['jobs'], [])
            with (out/'results.csv').open() as f:
                self.assertEqual(len(list(csv.DictReader(f))), len(rows))


class NoScaleTests(unittest.TestCase):
    def test_notebook_without_scaling_records_it_in_the_commands(self):
        with tempfile.TemporaryDirectory() as tmp:
            data, out = Path(tmp)/'d', Path(tmp)/'o'
            tiny_dataset(data)
            run_notebook(f"DATA={str(data)!r}\nOUTPUT={str(out)!r}\nSCALE_FLOOR=0\nN_SEEDS=2\nWORKERS=2\nTOP=5\n"
                         "SCHEMES={'true_cluster': 6, 'shape_cluster': 5}\n" + SMALL_GRID)
            rec = json.loads((out/'recommended.json').read_text())
            self.assertEqual(rec['tuned_on']['scale_floor'], 0)
            self.assertTrue(rec['schemes']['true_cluster']['best']['command'].endswith('--scale-floor 0'))
            self.assertIn('scale_floor=0', (out/'summary.md').read_text())

if __name__ == '__main__':
    unittest.main()
