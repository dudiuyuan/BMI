"""Validation notebook: re-scores the tuning winners on a second dataset."""
import csv
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT/'tests'))
from test_tuning import run_notebook, tiny_dataset, SMALL_GRID

VALIDATE = ROOT/'validate_dtw.ipynb'


class ValidationNotebookTests(unittest.TestCase):
    def test_rescores_chosen_configs_on_new_data(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            tiny_dataset(tmp/'base')
            run_notebook(f"DATA={str(tmp/'base')!r}\nOUTPUT={str(tmp/'tune')!r}\nSCALE_FLOOR=0.1\nN_SEEDS=3\n"
                         "WORKERS=2\nTOP=5\nSCHEMES={'true_cluster': 6, 'shape_cluster': 5}\n" + SMALL_GRID)
            (tmp/'new').mkdir(); tiny_dataset(tmp/'new'/'d')
            settings = (f"TUNING={str(tmp/'tune'/'results.csv')!r}\nSELECT_SCHEME='true_cluster'\nMIN_ARI=0.0\n"
                        f"DATA={str(tmp/'new'/'d')!r}\nOUTPUT={str(tmp/'val')!r}\nSCALE_FLOOR=0.1\nN_SEEDS=3\n"
                        "SCHEMES={'true_cluster': 6, 'shape_cluster': 5}\nWORKERS=2\n")
            ns = run_notebook(settings, VALIDATE)
            self.assertEqual(len(ns['chosen']), len([r for r in ns['tuned'] if r['method'] == 'dtw']))
            with (tmp/'val'/'results.csv').open() as f:
                rows = list(csv.DictReader(f))
            self.assertEqual(len(rows), 2*(len(ns['chosen']) + len(ns['references'])))
            text = (tmp/'val'/'summary.md').read_text()
            self.assertIn('Spearman', text); self.assertIn('baseline 排名', text)
            again = run_notebook(settings, VALIDATE)          # resume: nothing recomputed
            self.assertEqual(again['jobs'], [])

    def test_spearman(self):
        import json
        cells = [''.join(c['source']) for c in json.loads(VALIDATE.read_text())['cells'] if c['cell_type'] == 'code']
        src = next(c for c in cells if 'def spearman' in c)
        ns = {}
        exec(src.split('\nnew = ')[0], ns)
        self.assertAlmostEqual(ns['spearman']([1, 2, 3, 4], [10, 20, 30, 40]), 1)
        self.assertAlmostEqual(ns['spearman']([1, 2, 3, 4], [4, 3, 2, 1]), -1)
        self.assertAlmostEqual(ns['spearman']([1, 2, 2, 3], [1, 3, 2, 4]), 4.5/22.5**.5)  # tie -> average ranks


if __name__ == '__main__':
    unittest.main()
