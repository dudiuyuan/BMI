"""Generic CSV -> preprocessing -> distance -> k-medoids -> figures."""
import argparse
import hashlib
import json
import platform
import time
from pathlib import Path
import numpy as np
import pandas as pd
from preprocessing import preprocess
from distance_matrix import distance_matrix
from clustering import cluster
from plots import plot_clusters


def run(input_path, method, k, window, output, seed=11, scaling='zscore', max_missing=.2, n_init=5):
    start = time.perf_counter()
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    frame, log, audit = preprocess(pd.read_csv(input_path), scaling, max_missing)
    X = frame.iloc[:, 1:].to_numpy()
    t = time.perf_counter()
    D = distance_matrix(X, method, window)
    distance_seconds = time.perf_counter()-t
    labels, medoids, objective, iterations = cluster(D, k, seed, n_init)
    pd.DataFrame({'gene_id': frame.gene_id, 'cluster_id': labels}).to_csv(out/'assignments.csv', index=False)
    curves = frame.iloc[medoids].copy()
    curves.insert(0, 'cluster_id', np.arange(k))
    curves.to_csv(out/'medoids.csv', index=False)
    frame.to_csv(out/'processed_input.csv', index=False)
    audit.to_csv(out/'preprocessing_audit.csv', index=False)
    plot_clusters(frame, labels, medoids, out)
    config = dict(input=str(Path(input_path).resolve()), input_sha256=hashlib.sha256(Path(input_path).read_bytes()).hexdigest(),
                  method=method, k=k, window=window, seed=seed, n_init=n_init,
                  objective=objective, iterations=iterations, preprocessing=log,
                  distance_seconds=distance_seconds, total_seconds=time.perf_counter()-start,
                  python=platform.python_version(), platform=platform.platform(),
                  algorithm='self-implemented alternating k-medoids; squared local DTW cost, sqrt accumulated cost, no path normalization')
    (out/'run_config.json').write_text(json.dumps(config, indent=2), encoding='utf-8')
    return config


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', required=True)
    p.add_argument('--method', choices=['dtw', 'euclidean', 'pearson'], default='dtw')
    p.add_argument('--k', type=int, required=True)
    p.add_argument('--window', type=int, default=1, help='DTW index band radius; ignored by baselines')
    p.add_argument('--output', required=True)
    p.add_argument('--seed', type=int, default=11)
    p.add_argument('--scaling', choices=['zscore', 'none'], default='zscore')
    p.add_argument('--max-missing', type=float, default=.2)
    p.add_argument('--n-init', type=int, default=5)
    a = p.parse_args()
    try:
        print(json.dumps(run(a.input, a.method, a.k, a.window, a.output, a.seed, a.scaling, a.max_missing, a.n_init), indent=2))
    except ValueError as e:
        p.error(str(e))
