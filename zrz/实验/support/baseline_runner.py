"""Command-line runner shared by euclidean.py and pearson.py.

Upstream and downstream are the same as for dtw.py: the same CSV reader, the
same per-gene preprocessing (default centre + scale with floor 0.1), the same
distance_matrix.csv + run.json output, so diagnose.py and the clustering code
read a baseline experiment exactly like a DTW experiment. Only the pairwise
distance function differs.
"""
import argparse
import hashlib
import json
import time
from pathlib import Path

from dtw import read_matrix, preprocessing_config, prepare_expression, write_distance_matrix


def run_baseline(method, distance, definition, description, argv=None):
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--scale-floor', type=float, default=0.1,
                        help='Per-gene SD floor for centering/scaling (default 0.1, same as dtw.py); 0 uses input unchanged')
    parser.add_argument('--clustering', choices=['kmedoids','average','xqf'], default='xqf',
                        help='Preferred downstream clustering for diagnose (default xqf, same as dtw.py)')
    parser.add_argument('--clustering-seed', type=int, default=42,
                        help='Random seed saved for --clustering xqf (default 42)')
    parser.add_argument('--limit', type=int, help='Explicitly use first N genes (demonstration only)')
    args = parser.parse_args(argv)
    try:
        ids, times, values = read_matrix(args.input)
        original_count = len(ids)
        if args.limit is not None:
            if args.limit < 1:
                raise ValueError('limit must be positive')
            ids, values = ids[:args.limit], values[:args.limit]
        preprocessing = preprocessing_config(args.scale_floor)
        values = prepare_expression(values, preprocessing)
        args.output.mkdir(parents=True, exist_ok=True)
        started = time.perf_counter()
        metadata = {'distance_method': method, 'input': str(args.input.resolve()),
                    'input_sha256': hashlib.sha256(args.input.read_bytes()).hexdigest(),
                    'original_genes': original_count, 'selected_genes': len(ids), 'timepoints': times,
                    'distance_definition': definition,
                    'preprocessing': preprocessing, 'clustering': args.clustering,
                    **({'clustering_seed': args.clustering_seed} if args.clustering == 'xqf' else {}),
                    'expression_transform': ('centered per gene, divided by max(population SD, scale_floor)'
                                             if preprocessing['method'] == 'center_scale' else 'none; input values used unchanged'),
                    'labels_read': False, 'mode': 'matrix'}
        write_distance_matrix(args.output/'distance_matrix.csv', ids, values,
                              lambda i, j: distance(values[i], values[j]))
        metadata['pairs_computed'] = len(ids)*(len(ids)-1)//2
        metadata['elapsed_seconds'] = time.perf_counter()-started
        (args.output/'run.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
        print(f'Completed: {args.output.resolve()} ({metadata["elapsed_seconds"]:.3f}s)')
    except (ValueError, OSError, OverflowError, TypeError) as exc:
        parser.error(str(exc))
