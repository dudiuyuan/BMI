import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def plot_clusters(frame, labels, medoids, output):
    cols = list(frame.columns[1:])
    times = np.array([float(c[1:]) for c in cols])
    X = frame[cols].to_numpy()
    k = len(medoids)
    fig, axes = plt.subplots((k+1)//2, 2, figsize=(12, 2.8*((k+1)//2)), squeeze=False)
    for c, ax in enumerate(axes.flat):
        if c >= k:
            ax.set_visible(False)
            continue
        idx = np.flatnonzero(labels == c)
        # All genes are in heatmap; at most 100 evenly selected trajectories per panel.
        shown = idx[np.linspace(0, len(idx)-1, min(100, len(idx))).astype(int)]
        ax.plot(times, X[shown].T, color='steelblue', alpha=.1, lw=.6)
        ax.plot(times, X[medoids[c]], color='crimson', lw=2, label=str(frame.gene_id.iloc[medoids[c]]))
        ax.set_title(f'Cluster {c}: n={len(idx)}')
        ax.set_xlabel('Time (input units)')
        ax.set_ylabel('Processed expression')
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(output/'medoid_curves.png', dpi=160)
    plt.close(fig)
    order = np.lexsort((np.argmax(X, axis=1), labels))
    fig, ax = plt.subplots(figsize=(9, 9))
    limit = max(float(np.quantile(np.abs(X), .99)), .01)
    im = ax.imshow(X[order], aspect='auto', cmap='RdBu_r', vmin=-limit, vmax=limit, interpolation='nearest')
    ax.set_xticks(np.arange(len(times)), [f'{t:g}' for t in times], rotation=45)
    for boundary in np.flatnonzero(np.diff(labels[order])):
        ax.axhline(boundary+.5, color='black', lw=.5)
    ax.set_xlabel('Time (input units)')
    ax.set_ylabel('All retained genes, ordered by cluster and observed peak')
    fig.colorbar(im, ax=ax, label='Processed expression; 1% tails clipped for display')
    fig.tight_layout()
    fig.savefig(output/'heatmap.png', dpi=180)
    plt.close(fig)
