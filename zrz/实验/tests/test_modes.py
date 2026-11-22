"""dtw_modes.py: each fixed variant equals dtw.py with the same parameters written out."""
import csv
import json
import random
import subprocess
import sys
import tempfile
import unittest
from dataclasses import asdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import dtw
from dtw_modes import MODES, config, dtw_arguments


def run(script, *args):
    return subprocess.run([sys.executable, '-B', str(ROOT/script), *args], capture_output=True, text=True)


def matrix(folder):
    with (folder/'distance_matrix.csv').open() as f:
        return list(csv.reader(f))


class ModeTests(unittest.TestCase):
    def test_four_modes_with_expected_parameters(self):
        self.assertEqual(list(MODES), ['time', 'time-dynamics', 'diagonal-time-dynamics', 'diagonal-time'])
        self.assertEqual(asdict(config('time')), asdict(dtw.Config(time_penalty=2)))
        self.assertEqual(config('diagonal-time').window, 0.1667)
        self.assertEqual(config('time-dynamics').window, 0)
        self.assertEqual(config('diagonal-time-dynamics').dynamics, 0.3)
        for mode in MODES:
            self.assertEqual(dtw_arguments(mode)[-2:], ['--scale-floor', '0'])

    def test_each_mode_matches_dtw_py_and_records_variant(self):
        rng = random.Random(3)
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp); source = tmp/'in.csv'
            with source.open('w', newline='') as f:
                w = csv.writer(f); w.writerow(['gene_id'] + [f't{10*i}' for i in range(13)])
                for g in range(8):
                    w.writerow([f'g{g}'] + [round(rng.gauss(0, 1), 4) for _ in range(13)])
            for mode in MODES:
                a = run('dtw_modes.py', '--mode', mode, '--input', str(source), f'--output={tmp/mode}')
                self.assertEqual(a.returncode, 0, a.stderr)
                b = run('dtw.py', '--input', str(source), '--output', str(tmp/(mode+'_manual')), *dtw_arguments(mode))
                self.assertEqual(b.returncode, 0, b.stderr)
                self.assertEqual(matrix(tmp/mode), matrix(tmp/(mode+'_manual')))
                meta = json.loads((tmp/mode/'run.json').read_text())
                manual = json.loads((tmp/(mode+'_manual')/'run.json').read_text())
                self.assertEqual(meta['variant']['name'], mode)
                self.assertEqual(meta['config'], manual['config'])
                self.assertEqual(meta['preprocessing'], {'method': 'none', 'version': 1})
            # pass-through options still work
            c = run('dtw_modes.py', '--mode', 'time', '--input', str(source), '--output', str(tmp/'pt'),
                    '--limit', '5', '--clustering', 'kmedoids')
            self.assertEqual(c.returncode, 0, c.stderr)
            meta = json.loads((tmp/'pt'/'run.json').read_text())
            self.assertEqual((meta['selected_genes'], meta['clustering']), (5, 'kmedoids'))

    def test_fixed_parameters_cannot_be_overridden(self):
        for extra in (['--time-penalty', '5'], ['--scale-floor=0.1'], ['--config-file', 'x.json']):
            r = run('dtw_modes.py', '--mode', 'time', '--input', 'x.csv', '--output', 'y', *extra)
            self.assertNotEqual(r.returncode, 0)
            self.assertIn('fixed by --mode', r.stderr)
        self.assertNotEqual(run('dtw_modes.py', '--input', 'x.csv', '--output', 'y').returncode, 0)
        listing = run('dtw_modes.py', '--list')
        self.assertEqual(listing.returncode, 0)
        for mode in MODES:
            self.assertIn(mode, listing.stdout)


if __name__ == '__main__':
    unittest.main()
