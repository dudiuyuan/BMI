"""final_selection.ipynb runs end to end on real-format tuning results (plots not shown)."""
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT/'tests'))
from test_tuning import run_notebook, tiny_dataset, SMALL_GRID

SELECT = ROOT/'final_selection.ipynb'


class SelectionNotebookTests(unittest.TestCase):
    def test_runs_and_draws(self):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        shown = []
        plt.show = lambda *a, **k: (shown.append(plt.gcf()), plt.close('all'))
        nb = json.loads(SELECT.read_text())
        for c in nb['cells']:  # IPython magics are not Python
            c['source'] = [l for l in c['source'] if not l.lstrip().startswith('%')]
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp); tiny_dataset(tmp/'d')
            for floor, out in ((0, 'a'), (0.1, 'b')):
                run_notebook(f"DATA={str(tmp/'d')!r}\nOUTPUT={str(tmp/out)!r}\nSCALE_FLOOR={floor}\nN_SEEDS=2\n"
                             "WORKERS=2\nTOP=5\nSCHEMES={'true_cluster': 6, 'shape_cluster': 5}\n" + SMALL_GRID)
            clean = tmp/'sel.ipynb'; clean.write_text(json.dumps(nb))
            for compare in (repr(str(tmp/'b'/'results.csv')), 'None'):
                ns = run_notebook(f"RESULTS={str(tmp/'a'/'results.csv')!r}\nCOMPARE={compare}\n"
                                  "SCHEME='true_cluster'\nMETRIC='ari_mean'\n"
                                  f"NEW_RESULTS={str(tmp/'b'/'results.csv')!r}\nPICK_ABOVE=-1\n", clean)
                self.assertEqual(len(ns['picked']), len(ns['dtw']))     # threshold -1 picks everything
                self.assertEqual(ns['missing'], [])
                self.assertEqual(len(ns['sel']), len(ns['dtw']))
                self.assertGreater(len(ns['dtw']), 0)
                self.assertEqual(ns['y'], sorted(ns['y']))
                self.assertEqual(ns['value_at'](ns['y'], 100), max(ns['y']))
        self.assertEqual(len(shown), 4)  # two plots per run


if __name__ == '__main__':
    unittest.main()
