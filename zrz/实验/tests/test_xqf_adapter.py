"""Adapter tests use a stand-in module, so they check this side of the hand-off only."""
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import support.xqf_adapter as adapter

D = [[0, 1, 9, 9], [1, 0, 9, 9], [9, 9, 0, 1], [9, 9, 1, 0]]


def stub(body):
    folder = tempfile.mkdtemp()
    path = Path(folder) / 'k_medoids.py'
    path.write_text('def k_medoids(distances, k, seed=42, n_init=5, max_iter=100):\n' + body)
    return path


class XqfAdapterTest(unittest.TestCase):
    def setUp(self):
        self.saved = adapter.XQF_PATH

    def tearDown(self):
        adapter.XQF_PATH = self.saved

    def test_converts_labels_to_project_format(self):
        adapter.XQF_PATH = stub('    return [0, 0, 1, 1], [0, 2]\n')
        assignment, clusters = adapter.xqf_k_medoids(D, 2)
        self.assertEqual(assignment, [0, 0, 1, 1])
        self.assertEqual(clusters, [[0, 1], [2, 3]])

    def test_rejects_unassigned_genes(self):
        adapter.XQF_PATH = stub('    return [0, -1, 1, 1], [0, 2]\n')
        with self.assertRaises(ValueError):
            adapter.xqf_k_medoids(D, 2)

    def test_matrix_checked_before_her_code(self):
        adapter.XQF_PATH = stub('    raise AssertionError("should not be called")\n')
        with self.assertRaises(ValueError):
            adapter.xqf_k_medoids([[0, 1], [2, 0]], 1)


if __name__ == '__main__':
    unittest.main()
