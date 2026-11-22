"""Euclidean / Pearson baselines: distances, CLI output, and the shared downstream."""
import csv
import json
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from dtw import Config, align, read_matrix, prepare_expression, preprocessing_config
from euclidean import euclidean_distance
from pearson import pearson_distance
from diagnose import build_report


class DistanceTests(unittest.TestCase):
    def test_euclidean_hand_value_and_dtw_diagonal(self):
        self.assertAlmostEqual(euclidean_distance([0, 0], [3, 4]), 5)
        import random
        rng = random.Random(3)
        for _ in range(50):
            n = rng.randint(2, 13)
            x = [rng.gauss(0, 1) for _ in range(n)]
            y = [rng.gauss(0, 1) for _ in range(n)]
            # DTW forced onto the diagonal path is exactly the Euclidean distance.
            diagonal = align(x, y, Config(window=1e-9)).distance
            self.assertAlmostEqual(euclidean_distance(x, y), diagonal, places=12)

    def test_pearson_values_and_flat_curves(self):
        self.assertAlmostEqual(pearson_distance([1, 2, 3], [10, 20, 30]), 0)
        self.assertAlmostEqual(pearson_distance([1, 2, 3], [3, 2, 1]), 2)
        self.assertAlmostEqual(pearson_distance([1, 0, -1, 0], [0, 1, 0, -1]), 1)
        self.assertEqual(pearson_distance([5, 5, 5], [1, 1, 1]), 0)
        self.assertEqual(pearson_distance([5, 5, 5], [1, 2, 3]), 1)
        flat = prepare_expression([[0.1, 0.1, 0.1]], preprocessing_config(.1))[0]
        self.assertEqual(pearson_distance(flat, [1, 2, 3]), 1)

    def test_unequal_lengths_rejected(self):
        for f in (euclidean_distance, pearson_distance):
            with self.assertRaises(ValueError):
                f([1, 2], [1, 2, 3])


class PipelineTests(unittest.TestCase):
    def run_cli(self, script, source, output, *options):
        subprocess.run([sys.executable, '-B', str(ROOT/script), '--input', str(source),
                        '--output', str(output), *options], check=True, capture_output=True, text=True)
        meta = json.loads((output/'run.json').read_text())
        with (output/'distance_matrix.csv').open() as f:
            rows = list(csv.reader(f))
        return meta, rows

    def test_same_upstream_and_downstream_as_dtw(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            data = root/'data'; (data/'processed').mkdir(parents=True)
            source = data/'processed'/'centered.csv'
            source.write_text('gene_id,t0,t10,t20,t30\n'
                              'a,0,1,0,-1\nb,0,2,0,-2\nc,1,0,-1,0\nd,2,0,-2,0\n')
            (data/'ground_truth.csv').write_text(
                'gene_id,true_cluster\na,early\nb,early\nc,late\nd,late\n')
            _, _, raw = read_matrix(source)
            x = prepare_expression(raw, preprocessing_config(.1))
            dtw_meta, dtw_rows = self.run_cli('dtw.py', source, root/'dtw_run')
            for script, method, f in [('euclidean.py', 'euclidean', euclidean_distance),
                                      ('pearson.py', 'pearson', pearson_distance)]:
                meta, rows = self.run_cli(script, source, root/method)
                # Same upstream: identical input record, preprocessing, gene order.
                for key in ('input', 'input_sha256', 'preprocessing', 'timepoints',
                            'clustering', 'clustering_seed', 'mode', 'labels_read'):
                    self.assertEqual(meta[key], dtw_meta[key], key)
                self.assertEqual(meta['distance_method'], method)
                self.assertEqual(rows[0], dtw_rows[0])
                for i, row in enumerate(rows[1:]):
                    for j, v in enumerate(row[1:]):
                        self.assertAlmostEqual(float(v), 0.0 if i == j else f(x[i], x[j]), places=12)
                # Same downstream: diagnose reads it and clusters with xqf by default.
                report = build_report(root/method, 2)
                payload = json.loads(re.search(
                    r'<script id="report-data" type="application/json">(.*?)</script>',
                    report.read_text(), re.S).group(1))
                self.assertEqual(payload['distance_method'], method)
                self.assertEqual(payload['clustering'], 'xqf')
                self.assertEqual(payload['runs'][0]['paths'], {})
                self.assertEqual(payload['runs'][0]['metrics']['ari'], 1)


if __name__ == '__main__':
    unittest.main()
