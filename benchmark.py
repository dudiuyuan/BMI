import argparse
import csv
import json
import sys
import time
from collections import Counter
from dataclasses import asdict
from math import isfinite
from pathlib import Path

from metrics import clustering_scores

METHODS = ['standard', 'time', 'time-dynamics', 'diagonal-time-dynamics',
           'diagonal-time', 'euclidean', 'pearson']
TASKS = [('phase6', 6, 'true_cluster'), ('shape5', 5, 'shape_cluster')]
SEED = 42
N_INIT = 5
MAX_ITER = 100


def read_csv(path):
    with open(path, encoding='utf-8-sig', newline='') as file:
        return list(csv.DictReader(file))


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', encoding='utf-8', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def get_datasets(data_root, groups, data_folder=None):
    if data_folder is not None:
        return [{'dataset_id': 'custom', 'dataset_group': 'custom', 'factor': '',
                 'level_index': '', 'level': '', 'data_seed': '', 'folder': data_folder}]

    datasets = []
    for group in groups:
        if group == 'sweep':
            entries = read_csv(data_root / 'sweep' / 'sweep_index.csv')
            for entry in entries:
                relative = entry['relative_dir'].replace('\\', '/')
                datasets.append({
                    'dataset_id': f"sweep_{entry['factor']}_L{entry['level_index']}_S{entry['seed']}",
                    'dataset_group': 'sweep', 'factor': entry['factor'],
                    'level_index': entry['level_index'], 'level': entry['level'],
                    'data_seed': entry['seed'], 'folder': data_root / 'sweep' / relative,
                })
        else:
            folder = data_root / group
            with open(folder / 'metadata.json', encoding='utf-8') as file:
                metadata = json.load(file)
            datasets.append({'dataset_id': group, 'dataset_group': group, 'factor': '',
                             'level_index': '', 'level': '',
                             'data_seed': metadata['config']['seed'], 'folder': folder})
    return datasets


def get_config(dtw, mode):
    parameters = dtw.MODES[mode]
    if isinstance(parameters, tuple):
        parameters = parameters[1]
    return dtw.Config(**parameters)


def calculate_distances(values, times, method, dtw, euclidean, pearson):
    parameters = {}
    if method in dtw.MODES:
        config = get_config(dtw, method)
        parameters = asdict(config)

    start = time.perf_counter()
    n = len(values)
    distances = []
    for i in range(n):
        distances.append([0.0] * n)

    for i in range(n):
        for j in range(i + 1, n):
            if method == 'euclidean':
                value = euclidean.euclidean_distance(values[i], values[j], times)
            elif method == 'pearson':
                value = pearson.pearson_distance(values[i], values[j], times)
            else:
                value = dtw.dtw_distance(values[i], values[j], times, config)
            if not isfinite(value) or value < 0:
                raise ValueError('Distances must be finite and nonnegative')
            distances[i][j] = value
            distances[j][i] = value

    return distances, time.perf_counter() - start, parameters


def get_true_labels(folder, gene_ids, column):
    rows = read_csv(folder / 'ground_truth.csv')
    labels = {}
    for row in rows:
        gene = row['gene_id']
        if not gene.strip() or gene in labels:
            raise ValueError('Empty or duplicate gene ID in truth labels')
        labels[gene] = row[column]

    ordered = []
    for gene in gene_ids:
        if gene not in labels or not labels[gene]:
            raise ValueError(f'Missing true label: {gene}')
        ordered.append(labels[gene])
    return ordered, len(rows)


def run(args):
    sys.path.insert(0, str(args.team_root / 'zrz' / 'dtw'))
    sys.path.insert(0, str(args.team_root / 'xqf' / 'xqf_clustering'))
    import dtw
    import euclidean
    import pearson
    import k_medoids

    datasets = get_datasets(args.data_root, args.groups, args.data_folder)
    args.output.mkdir(parents=True, exist_ok=False)
    settings = {'team_root': str(args.team_root), 'data_root': str(args.data_root),
                'data_folder': str(args.data_folder) if args.data_folder else None,
                'groups': args.groups, 'methods': METHODS, 'tasks': TASKS,
                'seed': SEED, 'n_init': N_INIT, 'max_iter': MAX_ITER,
                'preprocessing': 'centered input; no extra scaling', 'python': sys.version,
                'parameters': {mode: asdict(get_config(dtw, mode)) for mode in METHODS[:5]}}
    with open(args.output / 'settings.json', 'w', encoding='utf-8') as file:
        json.dump(settings, file, indent=2, ensure_ascii=True)

    results = []
    for index, dataset in enumerate(datasets, 1):
        folder = dataset['folder']
        gene_ids, times, values = dtw.read_matrix(folder / 'processed' / 'centered.csv')
        for method in METHODS:
            distances, distance_seconds, parameters = calculate_distances(values, times, method, dtw, euclidean, pearson)
            output = args.output / 'predictions' / dataset['dataset_id'] / method
            output.mkdir(parents=True)
            with open(output / 'distance_matrix.csv', 'w', encoding='utf-8', newline='') as file:
                writer = csv.writer(file)
                writer.writerow(['gene_id'] + gene_ids)
                for gene, row in zip(gene_ids, distances):
                    writer.writerow([gene] + row)

            for task, k, column in TASKS:
                start = time.perf_counter()
                predicted, medoids = k_medoids.k_medoids(distances, k, seed=SEED, n_init=N_INIT, max_iter=MAX_ITER)
                clustering_seconds = time.perf_counter() - start
                if len(predicted) != len(gene_ids) or set(predicted) != set(range(k)):
                    raise ValueError('Invalid clustering output')

                # Read truth after clustering; match by ID because preprocessing can remove rows.
                truth, n_input = get_true_labels(folder, gene_ids, column)
                scores = clustering_scores(truth, predicted)
                prediction_file = output / f'{task}.csv'
                rows = []
                for i in range(len(gene_ids)):
                    rows.append({'gene_id': gene_ids[i], 'predicted_cluster': predicted[i],
                                 'true_label': truth[i], 'is_medoid': int(i in medoids)})
                write_csv(prediction_file, rows)

                composition = []
                for (label, cluster), count in sorted(Counter(zip(truth, predicted)).items()):
                    composition.append({'true_label': label, 'predicted_cluster': cluster, 'count': count})
                write_csv(output / f'{task}_composition.csv', composition)
                results.append({
                    'dataset_id': dataset['dataset_id'], 'dataset_group': dataset['dataset_group'],
                    'factor': dataset['factor'], 'level_index': dataset['level_index'], 'level': dataset['level'],
                    'data_seed': dataset['data_seed'], 'method': method, 'task': task, 'k': k,
                    'n_input_genes': n_input, 'n_genes': len(gene_ids), 'retained_fraction': len(gene_ids) / n_input,
                    'n_timepoints': len(times), 'ari': scores['ari'], 'nmi': scores['nmi'],
                    'distance_seconds': distance_seconds, 'clustering_seconds': clustering_seconds,
                    'clustering_seed': SEED, 'n_init': N_INIT, 'max_iter': MAX_ITER,
                    'parameters_json': json.dumps(parameters, sort_keys=True),
                    'prediction_file': prediction_file.relative_to(args.output).as_posix(),
                })
            write_csv(args.output / 'results.csv', results)
        print(f"dataset {index}/{len(datasets)}: {dataset['dataset_id']}; {len(results)} score records", flush=True)
    print('done:', args.output)


def get_args():
    parser = argparse.ArgumentParser(description='Benchmark simulated expression clusters')
    parser.add_argument('--team-root', type=Path, required=True, help='BMI directory containing zrz and xqf')
    parser.add_argument('--data-root', type=Path, help='prepared simulated directory containing baseline and sweep')
    parser.add_argument('--groups', nargs='+', choices=['smoke', 'baseline', 'challenge', 'sweep'], default=['baseline'])
    parser.add_argument('--data-folder', type=Path, help='run one dataset instead of --groups')
    parser.add_argument('--output', type=Path, required=True, help='new output folder')
    args = parser.parse_args()
    if args.data_root is None and args.data_folder is None:
        parser.error('provide --data-root or --data-folder')
    if args.output.exists():
        parser.error('output folder already exists; use a new name')
    if len(args.groups) != len(set(args.groups)):
        parser.error('do not repeat dataset groups')
    return args


def main():
    args = get_args()
    run(args)


if __name__ == '__main__':
    main()
