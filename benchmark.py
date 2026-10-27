"""Small supplementary integration benchmark: five one-factor challenges, three seeds."""
import argparse
import json
import time
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
import matplotlib.pyplot as plt
from simulation import generate
from preprocessing import preprocess
from distance_matrix import distance_matrix
from clustering import cluster


def benchmark(output):
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    conditions = dict(noise=[0., .15, .4], shift=[0., .08, .2], amplitude=[0., .4, .8],
                      missing=[0., .1, .25], stretch=[0., .3, .7])
    seeds = [101, 102, 103]
    rows = []
    # Warm compilation outside timing.
    distance_matrix(np.zeros((2, 25)), 'dtw', 2)
    for challenge, levels in conditions.items():
        for level in levels:
            for seed in seeds:
                frame, truth = generate(seed=seed, **{challenge: level})
                for scaling in ['zscore', 'none'] if challenge == 'amplitude' else ['zscore']:
                    clean, _, _ = preprocess(frame, scaling, max_missing=.4)
                    labels_true = truth.set_index('gene_id').loc[clean.gene_id, 'true_cluster']
                    for method, window in [('euclidean', 0), ('pearson', 0), ('dtw', 2), ('dtw', 24)]:
                        start = time.perf_counter()
                        D = distance_matrix(clean.iloc[:, 1:].to_numpy(), method, window)
                        labels, _, _, _ = cluster(D, 6, seed, n_init=3)
                        rows.append(dict(challenge=challenge, level=level, seed=seed, scaling=scaling,
                                         method=method, window=window, retained=len(clean),
                                         ari=adjusted_rand_score(labels_true, labels),
                                         nmi=normalized_mutual_info_score(labels_true, labels),
                                         seconds=time.perf_counter()-start))
    result = pd.DataFrame(rows)
    result.to_csv(out/'benchmark_runs.csv', index=False)
    result.groupby(['challenge', 'level', 'scaling', 'method', 'window'])[['ari', 'nmi', 'seconds']].agg(['mean', 'std']).to_csv(out/'benchmark_summary.csv')
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    for ax, (challenge, _) in zip(axes.flat, conditions.items()):
        subset = result[(result.challenge == challenge) & (result.scaling == 'zscore')]
        for (method, window), group in subset.groupby(['method', 'window']):
            stats = group.groupby('level').ari.agg(['mean', 'std'])
            ax.errorbar(stats.index, stats['mean'], yerr=stats['std'], marker='o', capsize=3,
                        label=f'{method} w={window}')
        ax.set_title(challenge)
        ax.set_xlabel('Perturbation level')
        ax.set_ylabel('ARI (mean +/- sample SD; 3 seeds)')
        ax.legend(fontsize=7)
    axes.flat[-1].set_visible(False)
    fig.tight_layout()
    fig.savefig(out/'benchmark_ari.png', dpi=160)
    plt.close(fig)
    (out/'experiment_config.json').write_text(json.dumps(dict(conditions=conditions, seeds=seeds, n=120, t=25,
        k=6, truth='phase-sensitive', max_missing=.4, n_init=3,
        limitation='early and late peaks share shape: DTW can merge them; zscore removes amplitude; supplementary integration benchmark, not Member D full benchmark'), indent=2))
    print('benchmark runs:', len(result))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', default='results/benchmarks')
    benchmark(p.parse_args().output)
