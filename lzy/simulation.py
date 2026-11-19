"""Phase-sensitive six-class integration benchmark; not a replacement for Member A's design."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd


def generate(n=120, t=25, seed=11, noise=.1, shift=0., amplitude=0., missing=0., stretch=0.):
    if n < 6 or t < 3 or min(noise, shift, amplitude, missing, stretch) < 0 or missing >= 1:
        raise ValueError('Invalid simulation parameters')
    rng = np.random.default_rng(seed)
    labels = np.arange(n) % 6
    rng.shuffle(labels)
    times = np.linspace(0, 1, t)
    X = []
    for label in labels:
        u = np.clip(times + rng.uniform(-shift, shift), 0, 1)
        # Positive exponent gives monotonic nonlinear warping with fixed endpoints.
        u = u ** np.exp(rng.uniform(-stretch, stretch))
        templates = [np.exp(-((u-.2)/.12)**2), np.exp(-((u-.8)/.12)**2),
                     np.sin(4*np.pi*u), np.exp(-((u-.5)/.08)**2),
                     -np.exp(-((u-.5)/.08)**2), np.zeros(t)]
        row = np.exp(rng.normal(0, amplitude))*templates[label] + rng.normal(0, noise, t)
        row[rng.random(t) < missing] = np.nan
        X.append(row)
    frame = pd.DataFrame(X, columns=[f't{i}' for i in range(t)])
    frame.insert(0, 'gene_id', [f'sim_{i:04d}' for i in range(n)])
    truth = pd.DataFrame({'gene_id': frame.gene_id, 'true_cluster': labels})
    return frame, truth


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', default='data/simulated')
    p.add_argument('--n', type=int, default=120)
    p.add_argument('--t', type=int, default=25)
    p.add_argument('--seed', type=int, default=11)
    for key, default in [('noise', .1), ('shift', 0.), ('amplitude', 0.), ('missing', 0.), ('stretch', 0.)]:
        p.add_argument('--'+key, type=float, default=default)
    a = p.parse_args()
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=True)
    frame, truth = generate(a.n, a.t, a.seed, a.noise, a.shift, a.amplitude, a.missing, a.stretch)
    frame.to_csv(out/'expression.csv', index=False)
    truth.to_csv(out/'ground_truth.csv', index=False)
