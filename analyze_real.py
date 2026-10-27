"""Unsupervised k/window selection followed by independent reference-set annotation."""
import json
import time
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import silhouette_score
from preprocessing import preprocess
from distance_matrix import distance_matrix
from clustering import cluster
from main import run
from biological_validation import validate


def analyze():
    out = Path('results/biological_validation')
    out.mkdir(parents=True, exist_ok=True)
    # Start from masked, pre-imputation gene matrix so missing-value audit is preserved.
    input_path = Path('data/gse3635/gene_log10ratio_with_missing.csv')
    clean, _, _ = preprocess(pd.read_csv(input_path))
    X = clean.iloc[:, 1:].to_numpy()
    rng = np.random.default_rng(11)
    sampled = np.sort(rng.choice(len(X), min(800, len(X)), replace=False))
    pd.DataFrame({'gene_id': clean.gene_id.iloc[sampled]}).to_csv(out/'selection_score_genes.csv', index=False)
    rows = []
    for window in [0, 1, 2]:
        start = time.perf_counter()
        D = distance_matrix(X, 'dtw', window)
        elapsed = time.perf_counter()-start
        for k in [3, 4, 5, 6, 7, 8]:
            labels, medoids, objective, iterations = cluster(D, k, seed=11, n_init=5)
            counts = np.bincount(labels, minlength=k)
            sample_labels = labels[sampled]
            score = float(silhouette_score(D[np.ix_(sampled, sampled)], sample_labels, metric='precomputed'))
            rows.append(dict(window=window, k=k, silhouette=score, minimum_cluster_size=int(counts.min()),
                             objective=objective, distance_seconds=elapsed, iterations=iterations))
        pd.DataFrame(rows).to_csv(out/'selection.csv', index=False)
        del D
        print(f'window {window} selection complete', flush=True)
    selection = pd.DataFrame(rows)
    # Rule fixed before reading Spellman labels; exclude clusters smaller than 1% of retained genes.
    eligible = selection[selection.minimum_cluster_size >= max(5, int(np.ceil(.01*len(X))))]
    if eligible.empty:
        raise ValueError('No candidate passed minimum cluster-size rule')
    best = eligible.sort_values(['silhouette', 'window', 'k'], ascending=[False, True, True]).iloc[0]
    k, window = int(best.k), int(best.window)
    selected = dict(k=k, window=window, silhouette=float(best.silhouette), seed=11,
                    rule='maximum sampled silhouette, min cluster size >=1%; ties smaller window then k',
                    candidate_k=[3,4,5,6,7,8], candidate_windows=[0,1,2], score_sample_size=len(sampled),
                    caveat='silhouettes use different distance geometries; descriptive heuristic, not independent biological model selection')
    (out/'selected_parameters.json').write_text(json.dumps(selected, indent=2))
    # Baselines use same k, data, seed, preprocessing and optimizer.
    for method in ['dtw', 'euclidean', 'pearson']:
        target = out/method
        run(input_path, method, k, window, target, seed=11)
    # Only now read reference labels. Secondary analyses must not be used to retune parameters.
    reference = pd.read_csv('data/references/spellman800.csv')
    for method in ['dtw', 'euclidean', 'pearson']:
        target = out/method
        print(method, validate(pd.read_csv(target/'assignments.csv'), reference, target,
                               pd.read_csv(target/'processed_input.csv')), flush=True)


if __name__ == '__main__':
    analyze()
