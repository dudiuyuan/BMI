"""Post hoc validation only. No reference information is passed to clustering."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


def bh(p):
    p = np.asarray(p, dtype=float)
    if not len(p):
        return p
    order = np.argsort(p)
    q = np.minimum.accumulate((p[order]*len(p)/np.arange(1, len(p)+1))[::-1])[::-1]
    result = np.empty_like(q)
    result[order] = np.minimum(q, 1)
    return result


def validate(assignments, reference, output, expression=None):
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    if assignments.gene_id.duplicated().any() or reference.gene_id.duplicated().any():
        raise ValueError('Duplicate gene IDs in assignments or reference')
    if assignments[['gene_id', 'cluster_id']].isna().any().any():
        raise ValueError('Missing assignment IDs/labels')
    universe = set(assignments.gene_id)
    known = set(reference.gene_id)
    hits = universe & known
    N, K = len(universe), len(hits)
    if not N:
        raise ValueError('Empty assignment universe')
    merged = assignments.merge(reference, on='gene_id', how='left', validate='one_to_one')
    merged['spellman_known'] = merged.gene_id.isin(known)
    merged.to_csv(out/'reference_matches.csv', index=False)
    pd.DataFrame({'gene_id': sorted(known-universe)}).to_csv(out/'reference_not_in_input.csv', index=False)
    rows = []
    phase_rows = []
    phases = sorted(reference.phase.dropna().unique()) if 'phase' in reference else []
    for c, group in assignments.groupby('cluster_id', sort=True):
        genes = set(group.gene_id)
        n, a = len(genes), len(genes & hits)
        table = [[a, n-a], [K-a, N-n-K+a]]
        odds, p = fisher_exact(table, alternative='greater')
        rows.append(dict(cluster_id=c, cluster_size=n, known_count=a, known_fraction=a/n,
                         background_size=N, background_known=K, background_fraction=K/N,
                         fold_enrichment=(a/n)/(K/N) if K else np.nan,
                         odds_ratio=odds, p_value=p, a=a, b=n-a, c=K-a, d=N-n-K+a))
        for phase in phases:
            phase_genes = set(reference.loc[reference.phase == phase, 'gene_id']) & universe
            P, h = len(phase_genes), len(genes & phase_genes)
            # Conditional phase validation: background = matched cell-cycle genes only.
            _, pp = fisher_exact([[h, a-h], [P-h, K-a-P+h]], alternative='greater')
            phase_rows.append(dict(cluster_id=c, phase=phase, phase_hits=h, matched_in_cluster=a,
                                   fraction_among_matched=h/a if a else np.nan,
                                   matched_background=K, phase_background=P, p_value=pp))
    result = pd.DataFrame(rows)
    result['q_value_bh'] = bh(result.p_value)
    result.to_csv(out/'cluster_enrichment.csv', index=False)
    if phase_rows:
        phase_table = pd.DataFrame(phase_rows)
        phase_table['q_value_bh'] = bh(phase_table.p_value)
        phase_table.to_csv(out/'phase_enrichment.csv', index=False)
        counts = phase_table.pivot(index='cluster_id', columns='phase', values='phase_hits')
        counts.to_csv(out/'phase_counts.csv')
        import matplotlib.pyplot as plt
        counts.plot.bar(stacked=True, figsize=(9, 5))
        plt.ylabel('Matched Spellman genes')
        plt.tight_layout()
        plt.savefig(out/'phase_distribution.png', dpi=160)
        plt.close()
    if expression is not None:
        peaks = expression.set_index('gene_id')
        peak = peaks.idxmax(axis=1).rename('observed_peak_time')
        merged.merge(peak, on='gene_id').to_csv(out/'peak_times.csv', index=False)
    summary = dict(input_genes=N, reference_genes=len(known), matched=K,
                   reference_coverage=K/len(known) if known else None,
                   background_fraction=K/N, significant_clusters=result.loc[result.q_value_bh < .05, 'cluster_id'].tolist(),
                   test='one-sided Fisher exact; BH across clusters; separate BH across all cluster-phase tests',
                   phase_background='only input genes with Spellman phase annotation',
                   interpretation='biological support, not clustering accuracy; not independent experimental replication')
    (out/'validation_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    return summary


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--assignments', required=True)
    p.add_argument('--reference', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--expression')
    args = p.parse_args()
    print(validate(pd.read_csv(args.assignments), pd.read_csv(args.reference), args.output,
                   pd.read_csv(args.expression) if args.expression else None))
