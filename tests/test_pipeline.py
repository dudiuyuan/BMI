import sys
import tempfile
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from dtw import dtw
from distance_matrix import distance_matrix
from clustering import cluster
from preprocessing import preprocess
from biological_validation import validate, bh


class AlgorithmTests(unittest.TestCase):
    def test_dtw_known_and_boundaries(self):
        self.assertEqual(dtw([1, 2, 3], [1, 2, 3]), 0)
        self.assertAlmostEqual(dtw([0, 1], [0, 2], 0), 1)
        self.assertEqual(dtw([0, 0, 1], [0, 1, 1], 1), 0)
        self.assertEqual(dtw([1, 1], [1, 1, 1], 1), 0)
        self.assertRaises(ValueError, dtw, [], [1])
        self.assertRaises(ValueError, dtw, [np.nan], [1])
        self.assertRaises(ValueError, dtw, [1], [1, 2, 3], 1)

    def test_zero_window_and_symmetry(self):
        X = np.random.default_rng(7).normal(size=(8, 7))
        euclidean = distance_matrix(X, 'euclidean')
        zero = distance_matrix(X, 'dtw', 0)
        np.testing.assert_allclose(euclidean, zero, rtol=1e-6)
        warped = distance_matrix(X, 'dtw', 2)
        np.testing.assert_allclose(warped, warped.T)
        self.assertTrue((warped <= zero+1e-6).all())

    def test_preprocessing_elapsed_time_and_duplicates(self):
        f = pd.DataFrame({'gene_id': ['a', 'a', 'b', 'c'], 't10': [10, 14, 5, np.nan],
                          't0': [0, 0, 5, np.nan], 't2': [np.nan, np.nan, 5, np.nan]})
        cleaned, log, audit = preprocess(f, 'none', .4)
        self.assertEqual(list(cleaned), ['gene_id', 't0', 't2', 't10'])
        np.testing.assert_allclose(cleaned.iloc[0, 1:].to_numpy(dtype=float), [0, 2.4, 12])
        self.assertEqual(log['dropped_genes'], 1)
        z, _, _ = preprocess(f, 'zscore', .4)
        np.testing.assert_allclose(z.iloc[1, 1:].to_numpy(dtype=float), 0)
        self.assertTrue((audit.missing_before >= 0).all())

    def test_cluster_clear_groups_and_ties(self):
        X = np.array([[0., 0.], [.1, .1], [8., 8.], [8.1, 8.1]])
        D = distance_matrix(X, 'euclidean')
        first = cluster(D, 2, seed=5)
        again = cluster(D, 2, seed=5)
        np.testing.assert_array_equal(first[0], again[0])
        self.assertEqual(first[0][0], first[0][1])
        self.assertNotEqual(first[0][0], first[0][2])
        labels, medoids, _, _ = cluster(np.zeros((5, 5)), 3)
        self.assertEqual(len(set(labels)), 3)
        self.assertEqual(len(set(medoids)), 3)

    def test_pearson_constants(self):
        D = distance_matrix([[1, 1], [2, 2], [0, 1], [1, 0]], 'pearson')
        self.assertEqual(D[0, 1], 0)
        self.assertEqual(D[0, 2], 1)
        self.assertAlmostEqual(D[2, 3], 2)

    def test_enrichment_background_and_fdr(self):
        a = pd.DataFrame({'gene_id': [f'g{i}' for i in range(10)], 'cluster_id': [0]*4+[1]*6})
        ref = pd.DataFrame({'gene_id': ['g0', 'g1', 'g2', 'g8', 'outside'], 'phase': ['G1']*3+['S']*2})
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parents[1]) as temp:
            summary = validate(a, ref, temp)
            self.assertEqual(summary['matched'], 4)
            result = pd.read_csv(Path(temp)/'cluster_enrichment.csv')
            self.assertAlmostEqual(result.p_value.iloc[0], hypergeom.sf(2, 10, 4, 4))
            self.assertEqual(int(result.iloc[0][['a', 'b', 'c', 'd']].sum()), 10)
            self.assertTrue((result.q_value_bh >= result.p_value-1e-14).all())
        np.testing.assert_allclose(bh([.01, .04, .03]), [.03, .04, .04])


if __name__ == '__main__':
    unittest.main()
