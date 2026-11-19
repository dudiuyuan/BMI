"""Acceptance checks for the reproducible inputs and handoff contract."""
from __future__ import annotations

import csv
import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from preprocessing import clean, interpolate, preprocess_file, read_expression  # noqa: E402
from prepare_real_data import parse_matrix  # noqa: E402
from simulation import CLASSES, generate, write_dataset  # noqa: E402


class SimulationTests(unittest.TestCase):
    def test_reproducibility_and_separate_truth(self):
        with tempfile.TemporaryDirectory() as tmp:
            a, b, c = (Path(tmp) / x for x in "abc")
            write_dataset(a, genes_per_class=3, timepoints=9, seed=7)
            write_dataset(b, genes_per_class=3, timepoints=9, seed=7)
            write_dataset(c, genes_per_class=3, timepoints=9, seed=8)
            for filename in ("expression.csv", "ground_truth.csv", "gene_parameters.csv", "metadata.json"):
                self.assertEqual((a / filename).read_bytes(), (b / filename).read_bytes())
            self.assertNotEqual((a / "expression.csv").read_bytes(), (c / "expression.csv").read_bytes())
            with (a / "expression.csv").open() as f:
                expression = list(csv.DictReader(f))
            with (a / "ground_truth.csv").open() as f:
                truth = list(csv.DictReader(f))
            self.assertEqual(len(expression), 18)
            self.assertEqual(list(expression[0]), ["gene_id", "t0", "t15", "t30", "t45", "t60", "t75", "t90", "t105", "t120"])
            self.assertEqual({r["gene_id"] for r in expression}, {r["gene_id"] for r in truth})
            self.assertEqual({r["true_cluster"] for r in truth}, set(CLASSES))
            self.assertNotIn("true_cluster", expression[0])

    def test_paired_one_factor_design(self):
        common = dict(genes_per_class=2, timepoints=13, seed=71, noise_sd=0.05,
                      max_shift=0.03, amplitude_sd=0.15, stretch_sd=0.08)
        low = generate(**common, missing_rate=0.0)
        high = generate(**common, missing_rate=0.31)
        self.assertEqual(low[2], high[2])
        for left, right in zip(low[1], high[1]):
            self.assertEqual(left["gene_id"], right["gene_id"])
            for col in low[0]:
                if right[col]:
                    self.assertEqual(left[col], right[col])
        noisy = generate(**dict(common, noise_sd=0.25), missing_rate=0)
        self.assertEqual(low[2], noisy[2])
        self.assertTrue(any(a != b for a, b in zip(low[1], noisy[1])))
        self.assertEqual([r["width"] for r in low[3]], [r["width"] for r in noisy[3]])

    def test_each_factor_changes_observations(self):
        base = dict(genes_per_class=3, timepoints=13, seed=29,
                    noise_sd=0.05, max_shift=0.03, amplitude_sd=0.15,
                    missing_rate=0.0, stretch_sd=0.08)
        reference = generate(**base)[1]
        for factor, value in (("noise_sd", 0.25), ("max_shift", 0.16),
                              ("amplitude_sd", 0.60), ("missing_rate", 0.31),
                              ("stretch_sd", 0.40)):
            with self.subTest(factor=factor):
                changed = generate(**dict(base, **{factor: value}))[1]
                self.assertEqual([x["gene_id"] for x in reference], [x["gene_id"] for x in changed])
                self.assertNotEqual(reference, changed)


class PreprocessingTests(unittest.TestCase):
    def test_interpolation_filter_and_constant_policy(self):
        times = [0.0, 2.0, 10.0, 20.0]
        interpolated = interpolate(times, [None, 2.0, None, 8.0])
        self.assertEqual((interpolated[0], interpolated[1], interpolated[3]), (2.0, 2.0, 8.0))
        self.assertAlmostEqual(interpolated[2], 2.0 + 6.0 * 8.0 / 18.0)
        kept, log = clean([("g1", [1.0, None, 1.0, 1.0]), ("g2", [None] * 4)], times, max_missing=0.4)
        self.assertEqual(len(kept), 1)
        self.assertEqual(log[1]["decision"], "all_missing")
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "input.csv"
            source.write_text("gene_id,t0,t2,t10,t20\ng1,1,,1,1\ng2,,,,\n", encoding="utf-8")
            report = preprocess_file(source, Path(tmp) / "processed")
            self.assertEqual((report["rows_kept"], report["imputed_cells"], report["constant_rows_after_imputation"]), (1, 1, 1))
            _, _, rows = read_expression(Path(tmp) / "processed" / "zscore.csv")
            self.assertEqual(rows[0][1], [0.0] * 4)
            self.assertEqual(json.loads((Path(tmp) / "processed" / "preprocessing_log.json").read_text())["rows_all_missing"], 1)

    def test_reject_unsorted_time_and_duplicate_gene(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "bad.csv"
            path.write_text("gene_id,t0,t10,t5\na,1,2,3\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "increase strictly"):
                read_expression(path)
            path.write_text("gene_id,t0,t5,t10\na,1,2,3\na,4,5,6\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "duplicate"):
                read_expression(path)


class ArtifactTests(unittest.TestCase):
    def test_real_metadata_and_processed_contract(self):
        expected = {"GSE3635": (13, 1), "GSE4987": (50, 2)}
        for accession, (samples, runs) in expected.items():
            with self.subTest(accession=accession):
                source = ROOT / accession / f"{accession}_series_matrix"
                _, _, rows, _, _, _, groups = parse_matrix(source)
                summary = json.loads((ROOT / "outputs" / "real" / accession / "data_preview_summary.json").read_text())
                self.assertEqual((len(rows), len(groups), summary["samples"]), (6228, runs, samples))
                self.assertEqual(summary["source_sha256"], hashlib.sha256(source.read_bytes()).hexdigest())
                self.assertEqual(summary["source_features"], summary["orf_features"] + summary["excluded_features"])
                for group in summary["run_outputs"]:
                    dirname = ROOT / "outputs" / "real" / accession / group["run"]
                    _, _, raw = read_expression(dirname / "genes_unprocessed.csv")
                    _, _, processed = read_expression(dirname / "processed" / "centered.csv")
                    log = json.loads((dirname / "processed" / "preprocessing_log.json").read_text())
                    self.assertEqual(len(raw), summary["orf_features"])
                    self.assertEqual(len(processed), log["rows_kept"])
                    self.assertTrue(all(all(v is not None for v in values) for _, values in processed))

    def test_sweep_index_paths_and_checksums(self):
        sweep = ROOT / "outputs" / "simulated" / "sweep"
        with (sweep / "sweep_index.csv").open() as f:
            entries = list(csv.DictReader(f))
        self.assertEqual(len(entries), 100)
        self.assertEqual({e["factor"] for e in entries}, {"noise_sd", "max_shift", "amplitude_sd", "missing_rate", "stretch_sd"})
        for entry in entries:
            path = sweep / entry["relative_dir"]
            for filename, key in (("expression.csv", "expression_sha256"), ("processed/centered.csv", "centered_sha256")):
                self.assertEqual(hashlib.sha256((path / filename).read_bytes()).hexdigest(), entry[key])
            _, _, matrix = read_expression(path / "processed" / "centered.csv")
            self.assertEqual(len(matrix), int(entry["genes_kept"]))
            self.assertTrue(all(all(v is not None for v in values) for _, values in matrix))


if __name__ == "__main__":
    unittest.main()
