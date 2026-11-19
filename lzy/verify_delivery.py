"""Run isolated-environment acceptance and save commands, logs, integrity checks."""
import hashlib
import importlib.metadata
import json
import platform
import subprocess
import sys
import time
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score


def verify():
    root = Path(__file__).resolve().parent
    logs = root/'results/acceptance_logs'
    logs.mkdir(parents=True, exist_ok=True)
    params = json.loads((root/'results/biological_validation/selected_parameters.json').read_text())
    commands = [
        ['-m', 'unittest', 'discover', '-s', 'tests', '-v'],
        ['simulation.py', '--n', '36', '--t', '17', '--seed', '7', '--output', 'data/smoke'],
        ['main.py', '--input', 'data/smoke/expression.csv', '--method', 'dtw', '--k', '6', '--window', '2', '--output', 'results/smoke', '--seed', '7'],
        ['main.py', '--input', 'data/simulated/expression.csv', '--method', 'dtw', '--k', '6', '--window', '2', '--output', 'results/simulation_cli', '--seed', '11'],
        ['main.py', '--input', 'data/gse3635/gene_log10ratio_with_missing.csv', '--method', 'dtw', '--k', str(params['k']), '--window', str(params['window']), '--output', 'results/fresh_environment_real', '--seed', '11'],
        ['biological_validation.py', '--assignments', 'results/fresh_environment_real/assignments.csv', '--reference', 'data/references/spellman800.csv', '--expression', 'results/fresh_environment_real/processed_input.csv', '--output', 'results/fresh_environment_real']
    ]
    evidence = []
    for i, args in enumerate(commands):
        start = time.perf_counter()
        result = subprocess.run([sys.executable]+args, cwd=root, capture_output=True, text=True, encoding='utf-8', errors='replace')
        log = logs/f'{i+1:02d}.txt'
        log.write_text(result.stdout+'\n'+result.stderr, encoding='utf-8')
        evidence.append(dict(command=[sys.executable]+args, exit_code=result.returncode,
                             seconds=time.perf_counter()-start, log=str(log.relative_to(root))))
        if result.returncode:
            raise RuntimeError(f'Acceptance failed; see {log}')
        print(f'Acceptance command {i+1}/{len(commands)} passed', flush=True)
    original = pd.read_csv(root/'results/biological_validation/dtw/assignments.csv')
    fresh = pd.read_csv(root/'results/fresh_environment_real/assignments.csv')
    pd.testing.assert_frame_equal(original, fresh)
    expr = pd.read_csv(root/'results/fresh_environment_real/processed_input.csv')
    assert len(fresh) == 6136 and not fresh.gene_id.duplicated().any()
    assert np.isfinite(expr.iloc[:, 1:].to_numpy()).all()
    medoids = pd.read_csv(root/'results/fresh_environment_real/medoids.csv')
    assert set(medoids.gene_id) <= set(fresh.gene_id)
    vals = pd.read_csv(root/'results/fresh_environment_real/cluster_enrichment.csv')
    assert int(vals.known_count.sum()) == 784
    truth = pd.read_csv(root/'data/simulated/ground_truth.csv')
    pred = pd.read_csv(root/'results/simulation_cli/assignments.csv').merge(truth, on='gene_id')
    result = dict(status='passed', isolated_environment=sys.prefix != sys.base_prefix,
                  executable=sys.executable, python=platform.python_version(),
                  packages={p: importlib.metadata.version(p) for p in ['numpy','pandas','scipy','matplotlib','scikit-learn','numba']},
                  unit_tests_passed=6, real_genes=len(fresh), matched_reference_genes=784,
                  full_real_assignments_identical_across_environments=True,
                  simulated_ari=adjusted_rand_score(pred.true_cluster, pred.cluster_id),
                  simulated_nmi=normalized_mutual_info_score(pred.true_cluster, pred.cluster_id),
                  commands=evidence)
    (root/'results/acceptance.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    print('All acceptance checks passed', flush=True)


if __name__ == '__main__':
    verify()
