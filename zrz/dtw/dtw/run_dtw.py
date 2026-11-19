#!/usr/bin/env python3
"""Run DTW with explicitly specified parameters; no method presets."""
import argparse
import csv
import hashlib
import json
import re
import time
from dataclasses import asdict
from html import escape
from pathlib import Path
from dtw import Config, align
from math import isfinite, fsum, sqrt

ROOT = Path(__file__).resolve().parent



def read_matrix(path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        reader = csv.reader(f)
        header = next(reader, [])
        if len(header) < 2 or header[0] != 'gene_id':
            raise ValueError('Expected gene_id,t0,t10,... header')
        if not all(re.fullmatch(r't[-+]?\d+(?:\.\d+)?', h) for h in header[1:]):
            raise ValueError('Invalid time-column names')
        times = [float(h[1:]) for h in header[1:]]
        if any(b <= a for a, b in zip(times, times[1:])):
            raise ValueError('Time columns must strictly increase')
        ids, values, seen = [], [], set()
        for line, row in enumerate(reader, 2):
            if len(row) != len(header):
                raise ValueError(f'Line {line}: wrong number of fields')
            if not row[0].strip() or row[0] in seen:
                raise ValueError(f'Line {line}: empty or duplicate gene_id')
            try:
                x = [float(v) for v in row[1:]]
            except ValueError as exc:
                raise ValueError(f'Line {line}: nonnumeric/missing data; use processed/centered.csv') from exc
            if not all(isfinite(v) for v in x):
                raise ValueError(f'Line {line}: nonfinite data')
            ids.append(row[0]); values.append(x); seen.add(row[0])
    if not ids:
        raise ValueError('No genes in input')
    return ids, times, values


def preprocessing_config(scale_floor):
    """A positive floor enables guarded, per-gene amplitude normalization."""
    if not isfinite(scale_floor) or scale_floor < 0:
        raise ValueError('scale_floor must be finite and nonnegative; 0 disables normalization')
    if scale_floor == 0:
        return {'method': 'none', 'version': 1}
    return {'method': 'center_scale', 'scale_floor': scale_floor, 'version': 1}


def prepare_expression(values, preprocessing):
    """Replay the saved transform before both DTW and diagnostic plotting.

    No labels or other genes are used. This makes pair, subset, and full-matrix
    runs agree. Missing preprocessing metadata means legacy unchanged input.
    """
    spec = preprocessing or {'method': 'none', 'version': 1}
    if spec.get('version', 1) != 1:
        raise ValueError('Unsupported expression preprocessing version')
    if spec.get('method') == 'none':
        return [list(row) for row in values]
    if spec.get('method') != 'center_scale':
        raise ValueError('Unknown expression preprocessing method')
    floor = spec['scale_floor']
    if not isfinite(floor) or floor <= 0:
        raise ValueError('Saved center_scale preprocessing requires a positive scale_floor')
    result = []
    for row in values:
        mean = fsum(row) / len(row)
        centered = [v - mean for v in row]
        sd = sqrt(fsum(v*v for v in centered) / len(centered))
        result.append([v / max(sd, floor) for v in centered])
    return result


def plot_pair(path, ids, times, x, y, result):
    """Standalone SVG: original time positions and the exact warping links."""
    lo, hi = min(x+y), max(x+y)
    span = hi-lo or 1
    tspan = times[-1]-times[0] or 1
    px = lambda i: 65+850*(times[i]-times[0])/tspan
    py = lambda v, base: base-105*(v-lo)/span
    parts = ['<svg xmlns="http://www.w3.org/2000/svg" width="980" height="470" viewBox="0 0 980 470">',
             '<rect width="980" height="470" fill="#faf9f6"/>',
             '<g font-family="sans-serif" font-size="14" fill="#202c35">',
             '<text x="45" y="32" font-size="22">DTW alignment</text>']
    for i, j in result.path:
        parts.append(f'<line x1="{px(i)}" y1="{py(x[i],190)}" x2="{px(j)}" y2="{py(y[j],375)}" stroke="#c5cbd0"/>')
    for series, base, color, label in [(x,190,'#007f82',ids[0]),(y,375,'#c35a35',ids[1])]:
        points = ' '.join(f'{px(i)},{py(v,base)}' for i,v in enumerate(series))
        parts.append(f'<polyline points="{points}" fill="none" stroke="{color}" stroke-width="3"/>')
        parts.append(f'<text x="45" y="{base-120}" fill="{color}">{escape(label)}</text>')
        for i,v in enumerate(series):
            parts.append(f'<circle cx="{px(i)}" cy="{py(v,base)}" r="3" fill="{color}"/>')
    for i,t in enumerate(times):
        parts.append(f'<text x="{px(i)}" y="410" text-anchor="middle">{t:g}</text>')
    parts.append(f'<text x="45" y="448">Time (input coordinates) · distance = {result.distance:.6g} · max normalized lag = {result.diagnostics["max_time_displacement"]:.3g}</text></g></svg>')
    path.write_text('\n'.join(parts), encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--scale-floor', type=float, default=0.1,
                        help='Per-gene SD floor for centering/scaling (default 0.1); 0 uses input unchanged')
    parser.add_argument('--clustering', choices=['kmedoids','average'], default='kmedoids',
                        help='Preferred downstream clustering for diagnose; does not change DTW distances')
    parser.add_argument('--window', type=float, default=0.0, help='Normalized time window; 0 disables constraint')
    parser.add_argument('--step-penalty', type=float, default=0.0)
    parser.add_argument('--time-penalty', type=float, default=0.0)
    parser.add_argument('--power', type=float, default=2.0)
    parser.add_argument('--dynamics', type=float, default=0.0)
    parser.add_argument('--open-penalty', type=float, default=0.0)
    parser.add_argument('--extend-penalty', type=float, default=0.0)
    parser.add_argument('--pair', nargs=2, metavar=('GENE_A','GENE_B'), help='Only align this pair of gene IDs')
    parser.add_argument('--limit', type=int, help='Explicitly use first N genes (demonstration only)')
    parser.add_argument('--allow-large', action='store_true', help='Allow full pairwise matrix above 1000 genes')
    parser.add_argument('--config-file', type=Path, help='Config JSON or run.json; run.json also reuses saved preprocessing and clustering')
    args = parser.parse_args()
    try:
        ids, times, values = read_matrix(args.input)
        original_count = len(ids)
        if args.limit is not None:
            if args.limit < 1:
                raise ValueError('limit must be positive')
            ids, values = ids[:args.limit], values[:args.limit]
        saved = json.loads(args.config_file.read_text(encoding='utf-8')) if args.config_file else None
        if saved is not None and not isinstance(saved, dict):
            raise ValueError('Config file must contain a JSON object')
        cfg = Config(**saved.get('config', saved)) if saved is not None else Config(
            window=args.window, step_penalty=args.step_penalty,
            time_penalty=args.time_penalty, power=args.power, dynamics=args.dynamics,
            open_penalty=args.open_penalty, extend_penalty=args.extend_penalty)
        # Old run.json files explicitly reproduce their old, unscaled inputs.
        saved_run = saved is not None and 'config' in saved
        preprocessing = (saved.get('preprocessing', {'method': 'none', 'version': 1})
                         if saved_run else preprocessing_config(args.scale_floor))
        clustering = saved.get('clustering', 'average') if saved_run else args.clustering
        if clustering not in ('average', 'kmedoids'):
            raise ValueError('Unknown saved clustering algorithm')
        values = prepare_expression(values, preprocessing)
        if not args.pair and len(ids) > 1000 and not args.allow_large:
            raise ValueError('Over 1000 genes: use --pair, --limit N, or --allow-large. Full matrix may take a long time.')
        if args.pair and any(g not in ids for g in args.pair):
            raise ValueError('Pair gene ID not present in selected input')
        args.output.mkdir(parents=True, exist_ok=True)
        started = time.perf_counter()
        metadata = {'config':asdict(cfg), 'input':str(args.input.resolve()),
                    'input_sha256':hashlib.sha256(args.input.read_bytes()).hexdigest(),
                    'original_genes':original_count, 'selected_genes':len(ids), 'timepoints':times,
                    'distance_definition':'sqrt(minimum cumulative local and transition cost)',
                    'kernel_sha256':hashlib.sha256((ROOT/'dtw.py').read_bytes()).hexdigest(),
                    'preprocessing':preprocessing, 'clustering':clustering,
                    'expression_transform':('centered per gene, divided by max(population SD, scale_floor)'
                                            if preprocessing['method']=='center_scale' else 'none; input values used unchanged'),
                    'labels_read':False, 'mode':'pair' if args.pair else 'matrix'}
        if args.pair:
            a, b = [ids.index(g) for g in args.pair]
            result = align(values[a],values[b],cfg,times_x=times,times_y=times)
            with (args.output/'alignment.csv').open('w',newline='') as f:
                w=csv.writer(f); w.writerow(['index_a','index_b','time_a','time_b','value_a','value_b'])
                for i,j in result.path:
                    w.writerow([i,j,times[i],times[j],values[a][i],values[b][j]])
            plot_pair(args.output/'alignment.svg',args.pair,times,values[a],values[b],result)
            metadata.update(pair=args.pair, cost=result.cost, distance=result.distance, diagnostics=result.diagnostics)
        else:
            # Stream rows to disk: O(N) matrix-output memory, compute each unordered pair once.
            # Lower triangle is kept until its row is written, using a temporary disk store.
            import sqlite3
            import tempfile
            with tempfile.TemporaryDirectory(prefix='dtw-') as temp:
                db=sqlite3.connect(str(Path(temp)/'distances.db'))
                db.execute('CREATE TABLE d (i INTEGER,j INTEGER,v REAL,PRIMARY KEY(i,j))')
                try:
                    with (args.output/'distance_matrix.csv').open('w',newline='') as f:
                        w=csv.writer(f); w.writerow(['gene_id']+ids)
                        for i,x in enumerate(values):
                            lower=dict(db.execute('SELECT j,v FROM d WHERE i=?',(i,)))
                            row=[lower[j] for j in range(i)]+[0.0]
                            db.execute('DELETE FROM d WHERE i=?',(i,))
                            for j in range(i+1,len(ids)):
                                d=align(x,values[j],cfg,times_x=times,times_y=times,return_path=False).distance
                                row.append(d); db.execute('INSERT INTO d VALUES (?,?,?)',(j,i,d))
                            w.writerow([ids[i]]+row)
                            db.commit()
                            if i % 100 == 0:
                                print(f'Rows: {i+1}/{len(ids)}',flush=True)
                finally:
                    db.close()
            metadata['pairs_computed']=len(ids)*(len(ids)-1)//2
        metadata['elapsed_seconds']=time.perf_counter()-started
        (args.output/'run.json').write_text(json.dumps(metadata,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
        print(f'Completed: {args.output.resolve()} ({metadata["elapsed_seconds"]:.3f}s)')
    except (ValueError,OSError,OverflowError,TypeError) as exc:
        parser.error(str(exc))


if __name__ == '__main__':
    main()
