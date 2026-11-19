"""Shared preprocessing for synthetic and real gene-by-time matrices."""
import re
import numpy as np
import pandas as pd


def preprocess(frame, scaling='zscore', max_missing=0.2):
    if not 0 <= max_missing < 1:
        raise ValueError('max_missing must be in [0, 1)')
    if 'gene_id' not in frame or frame.gene_id.isna().any():
        raise ValueError('A nonempty gene_id column is required')
    frame = frame.copy()
    frame['gene_id'] = frame.gene_id.astype(str).str.strip()
    if (frame.gene_id == '').any():
        raise ValueError('Empty gene_id')
    cols = [c for c in frame if c != 'gene_id']
    if len(cols) < 2 or any(not re.fullmatch(r't-?\d+(?:\.\d+)?', c) for c in cols):
        raise ValueError('Time columns must be numeric t0,t10,...; at least two are required')
    cols.sort(key=lambda c: float(c[1:]))
    times = np.array([float(c[1:]) for c in cols])
    if len(set(times)) != len(times):
        raise ValueError('Duplicate numeric times')
    # Invalid text is an error; blank/NaN and +/-infinity are missing.
    values = frame[cols].apply(pd.to_numeric, errors='raise').replace([np.inf, -np.inf], np.nan)
    values.index = frame.gene_id
    before = len(values)
    values = values.groupby(level=0, sort=True).median()
    missing = values.isna().sum(axis=1)
    keep = (missing / len(cols) <= max_missing) & (values.notna().sum(axis=1) >= 2)
    audit = pd.DataFrame({'gene_id': values.index, 'missing_before': missing.to_numpy(),
                          'retained': keep.to_numpy()})
    kept = values.loc[keep].copy()
    X = kept.to_numpy(dtype=float, copy=True)
    for row in X:
        good = np.isfinite(row)
        # Interpolate on actual elapsed time; nearest value at either endpoint.
        row[~good] = np.interp(times[~good], times[good], row[good])
    if not len(X):
        raise ValueError('No genes remain after missing-value filtering')
    std = X.std(axis=1)
    constant = std <= 1e-12
    if scaling == 'zscore':
        X = (X - X.mean(axis=1, keepdims=True)) / np.where(constant, 1, std)[:, None]
    elif scaling != 'none':
        raise ValueError('scaling must be zscore or none')
    out = pd.DataFrame(X, columns=cols)
    out.insert(0, 'gene_id', kept.index)
    log = dict(input_rows=before, unique_genes=len(values), retained_genes=len(out),
               dropped_genes=int((~keep).sum()), imputed_cells=int(missing.loc[keep].sum()),
               constant_genes=int(constant.sum()), scaling=scaling, max_missing=max_missing,
               duplicate_rule='per-time median', interpolation='linear elapsed time; nearest endpoints',
               time_columns=cols)
    return out, log, audit
