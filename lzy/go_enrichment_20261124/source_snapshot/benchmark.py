"""Run the team's frozen distance and clustering code on existing datasets.

Python 3.10+, standard library only. No algorithm or input data is rewritten.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
import math
import os
import platform
import shutil
import sys
import time
from collections import Counter
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path

from metrics import clustering_scores

HERE = Path(__file__).resolve().parent
DEFAULT_TEAM = HERE.parent / "小组成员压缩包" / "BMI"
METHODS = ["standard", "time", "time-dynamics", "diagonal-time-dynamics", "diagonal-time", "euclidean", "pearson"]
TASKS = {"phase6": (6, "true_cluster"), "shape5": (5, "shape_cluster")}
CODE_PATHS = ["zrz/dtw/dtw.py", "zrz/dtw/euclidean.py", "zrz/dtw/pearson.py", "xqf/xqf_clustering/k_medoids.py"]
FIELDS = ["dataset_id", "dataset_group", "factor", "level_index", "level", "data_seed", "method", "task", "k",
          "n_input_genes", "n_genes", "retained_fraction", "n_timepoints", "ari", "nmi", "distance_seconds",
          "clustering_seconds", "clustering_seed", "n_init", "max_iter", "input_sha256", "matrix_sha256",
          "distance_cache_hit", "clustering_cache_hit", "distance_measurement_id", "clustering_measurement_id",
          "parameters_json", "prediction_file"]


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def signature(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, ensure_ascii=True).encode()).hexdigest()


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def source_lock(team, output):
    originals = {name: digest(team / name) for name in CODE_PATHS}
    target = output / "source_snapshot"
    manifest = target / "source_manifest.json"
    if manifest.exists():
        saved = json.loads(manifest.read_text(encoding="utf-8"))
        if saved["sha256"] != originals:
            raise ValueError("Team code changed after the run was frozen. Use a new output directory.")
        for name in CODE_PATHS:
            if digest(target / name) != originals[name]:
                raise ValueError(f"Frozen source changed: {name}")
    else:
        for name in CODE_PATHS:
            destination = target / name
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(team / name, destination)
        write_json(manifest, {"team_root": str(team), "sha256": originals})
    dtw = load_module("benchmark_frozen_dtw", target / CODE_PATHS[0])
    # Baselines import run from the adjacent final dtw.py, not lzy or the old code.
    sys.modules["dtw"] = dtw
    euclidean = load_module("benchmark_frozen_euclidean", target / CODE_PATHS[1])
    pearson = load_module("benchmark_frozen_pearson", target / CODE_PATHS[2])
    cluster = load_module("benchmark_frozen_kmedoids", target / CODE_PATHS[3])
    return dtw, euclidean, pearson, cluster, originals


def collect_datasets(team, groups, dtw):
    simulated = team / "wpj" / "predata" / "数据准备" / "outputs" / "simulated"
    datasets = []
    for group in groups:
        if group == "sweep":
            entries = read_csv(simulated / "sweep" / "sweep_index.csv")
            if len(entries) != 100:
                raise ValueError("Expected 100 formal sweep cases")
            for entry in entries:
                folder = simulated / "sweep" / Path(entry["relative_dir"])
                datasets.append(dict(dataset_id=f"sweep_{entry['factor']}_L{entry['level_index']}_S{entry['seed']}",
                                     dataset_group="sweep", factor=entry["factor"], level_index=int(entry["level_index"]),
                                     level=float(entry["level"]), data_seed=int(entry["seed"]),
                                     folder=folder, expected_input_sha256=entry["centered_sha256"],
                                     expected_n_genes=int(entry["genes_kept"])))
        else:
            folder = simulated / group
            meta = json.loads((folder / "metadata.json").read_text(encoding="utf-8"))
            datasets.append(dict(dataset_id=group, dataset_group=group, factor="", level_index="", level="",
                                 data_seed=meta["config"]["seed"], folder=folder))
    # Input checks do not read truth; truth is joined only after predictions.
    manifests = []
    for dataset in datasets:
        source = dataset["folder"] / "processed" / "centered.csv"
        actual_sha = digest(source)
        if dataset.get("expected_input_sha256", actual_sha) != actual_sha:
            raise ValueError(f"Input checksum mismatch: {source}")
        ids, times, values = dtw.read_matrix(source)
        if len(ids) != dataset.get("expected_n_genes", len(ids)):
            raise ValueError(f"Retained gene count mismatch: {source}")
        dataset.update(input=source, input_sha256=actual_sha,
                       truth_sha256=digest(dataset["folder"] / "ground_truth.csv"),
                       n_genes=len(ids), n_timepoints=len(times))
        manifests.append({key: str(value) if isinstance(value, Path) else value for key, value in dataset.items()})
    return datasets, manifests


def validate_matrix(matrix, ids):
    n = len(ids)
    if len(matrix) != n or any(len(row) != n for row in matrix):
        raise ValueError("Wrong matrix shape")
    if any(not math.isfinite(value) or value < 0 for row in matrix for value in row):
        raise ValueError("Invalid matrix values")
    if any(matrix[i][i] != 0 for i in range(n)):
        raise ValueError("Nonzero diagonal")
    if any(matrix[i][j] != matrix[j][i] for i in range(n) for j in range(i)):
        raise ValueError("Asymmetric distance matrix")


def distance_matrix(ids, values, times, distance):
    start = time.perf_counter()
    matrix = [[0.0] * len(ids) for _ in ids]
    for i in range(len(ids)):
        for j in range(i + 1, len(ids)):
            matrix[i][j] = matrix[j][i] = distance(values[i], values[j], times)
    elapsed = time.perf_counter() - start
    validate_matrix(matrix, ids)
    return matrix, elapsed


def get_matrix(output, dataset, method, params, dtw, euclidean, pearson, source_hashes):
    # Byte-identical shared baseline cases safely share work, but not replication counts.
    key = signature({"input_sha256": dataset["input_sha256"], "method": method,
                     "params": params, "source_sha256": source_hashes, "format_version": 1})
    cache = output / "distance_cache" / key
    matrix_file, metadata_file = cache / "distance_matrix.csv.gz", cache / "run.json"
    ids, times, values = dtw.read_matrix(dataset["input"])
    if metadata_file.exists():
        meta = json.loads(metadata_file.read_text(encoding="utf-8"))
        if digest(matrix_file) != meta["compressed_file_sha256"]:
            raise ValueError(f"Corrupt cached matrix: {cache}")
        with gzip.open(matrix_file, "rt", encoding="utf-8", newline="") as handle:
            rows = list(csv.reader(handle))
        if rows[0] != ["gene_id", *ids] or [row[0] for row in rows[1:]] != ids:
            raise ValueError("Cached matrix IDs do not match input")
        matrix = [[float(value) for value in row[1:]] for row in rows[1:]]
        validate_matrix(matrix, ids)
        return ids, matrix, meta, True
    if method in dtw.MODES:
        config = dtw.Config(**params)
        distance = lambda x, y, t: dtw.dtw_distance(x, y, t, config)
    else:
        distance = euclidean.euclidean_distance if method == "euclidean" else pearson.pearson_distance
    matrix, seconds = distance_matrix(ids, values, times, distance)
    cache.mkdir(parents=True, exist_ok=True)
    with gzip.open(matrix_file, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["gene_id", *ids])
        writer.writerows([gene, *row] for gene, row in zip(ids, matrix))
    with gzip.open(matrix_file, "rb") as handle:
        matrix_sha256 = hashlib.sha256(handle.read()).hexdigest()
    meta = dict(measurement_id=key, method=method, params=params, input_sha256=dataset["input_sha256"],
                first_dataset_id=dataset["dataset_id"], n_genes=len(ids), n_timepoints=len(times),
                distance_seconds=seconds, matrix_sha256=matrix_sha256,
                compressed_file_sha256=digest(matrix_file),
                timing="matrix allocation + all unordered pair calls; excludes input, validation, truth, CSV/gzip writes and plotting",
                preprocessing="centered input unchanged", source_sha256=source_hashes)
    write_json(metadata_file, meta)
    return ids, matrix, meta, False


def get_prediction(output, matrix_meta, matrix, k, cluster, cluster_options, cluster_sha256):
    key = signature({"distance_measurement_id": matrix_meta["measurement_id"], "k": k,
                     "options": cluster_options, "code_sha256": cluster_sha256})
    filename = output / "clustering_cache" / f"{key}.json"
    if filename.exists():
        result = json.loads(filename.read_text(encoding="utf-8"))
        return result, True
    start = time.perf_counter()
    labels, medoids = cluster.k_medoids(matrix, k, **cluster_options)
    elapsed = time.perf_counter() - start
    if len(labels) != len(matrix) or set(labels) != set(range(k)) or len(set(medoids)) != k:
        raise ValueError("Clustering returned incomplete labels or medoids")
    result = dict(measurement_id=key, labels=labels, medoids=medoids, k=k, clustering_seconds=elapsed,
                  cluster_options=cluster_options,
                  objective=sum(matrix[i][medoids[label]] for i, label in enumerate(labels)))
    write_json(filename, result)
    return result, False


def align_truth(folder, ids, label_column):
    truth = read_csv(folder / "ground_truth.csv")
    lookup = {row["gene_id"]: row for row in truth}
    if len(lookup) != len(truth):
        raise ValueError("Duplicate ground truth gene IDs")
    if any(gene not in lookup or not lookup[gene].get(label_column) for gene in ids):
        raise ValueError("Missing ground truth for a retained gene")
    return [lookup[gene][label_column] for gene in ids], len(truth)


def run(args):
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    dtw, euclidean, pearson, cluster, source_hashes = source_lock(args.team_root.resolve(), output)
    datasets, manifest = collect_datasets(args.team_root.resolve(), args.groups, dtw)
    parameter_table = {method: asdict(dtw.Config(**dtw.MODES[method][1])) if method in dtw.MODES else {}
                       for method in METHODS}
    cluster_options = dict(seed=args.clustering_seed, n_init=args.n_init, max_iter=args.max_iter)
    locked = dict(version=1, team_root=str(args.team_root.resolve()), preprocessing="centered input unchanged",
                  methods=parameter_table, clustering="xqf alternating random multi-start k-medoids",
                  cluster_options=cluster_options, tasks=TASKS, nmi_average_method="arithmetic", source_sha256=source_hashes,
                  runner_sha256=digest(Path(__file__)), metrics_sha256=digest(HERE / "metrics.py"))
    config_file = output / "experiment_config.json"
    if config_file.exists():
        if json.loads(config_file.read_text(encoding="utf-8")) != json.loads(json.dumps(locked)):
            raise ValueError("Locked settings/code differ; use a new output directory")
    else:
        write_json(config_file, locked)
        write_json(output / "environment.json", dict(python=sys.version, platform=platform.platform(),
                   processor=platform.processor(), cpu_count=os.cpu_count(),
                   implementation=platform.python_implementation(), parallelism="single process, sequential methods",
                   created_utc=datetime.now(timezone.utc).isoformat(), argv=sys.argv))
    old_manifest_file = output / "dataset_manifest.json"
    merged_manifest = {item["dataset_id"]: item for item in json.loads(old_manifest_file.read_text(encoding="utf-8"))} if old_manifest_file.exists() else {}
    for item in manifest:
        if item["dataset_id"] in merged_manifest and merged_manifest[item["dataset_id"]] != item:
            raise ValueError(f"Dataset changed: {item['dataset_id']}")
        merged_manifest[item["dataset_id"]] = item
    write_json(old_manifest_file, list(merged_manifest.values()))
    results_file = output / "results.csv"
    results = read_csv(results_file) if results_file.exists() else []
    existing = {(row["dataset_id"], row["method"], row["task"]) for row in results}
    started = time.perf_counter()
    total = len(datasets) * len(METHODS) * len(TASKS)
    processed = 0
    for index, dataset in enumerate(datasets, 1):
        for method in METHODS:
            wanted = [(task, k, label) for task, (k, label) in TASKS.items()
                      if (dataset["dataset_id"], method, task) not in existing]
            if not wanted:
                processed += len(TASKS)
                continue
            ids, matrix, matrix_meta, distance_reused = get_matrix(output, dataset, method, parameter_table[method],
                                                                 dtw, euclidean, pearson, source_hashes)
            for task, k, label_column in wanted:
                prediction, cluster_reused = get_prediction(output, matrix_meta, matrix, k, cluster, cluster_options,
                                                            source_hashes[CODE_PATHS[3]])
                truth, n_input = align_truth(dataset["folder"], ids, label_column)
                metric_values = clustering_scores(truth, prediction["labels"])
                if set(truth) != set(row[label_column] for row in read_csv(dataset["folder"] / "ground_truth.csv")):
                    raise ValueError("A true class disappeared after filtering; revisit K before scoring")
                prediction_relative = Path("predictions") / dataset["dataset_id"] / method / f"{task}.csv"
                prediction_rows = [dict(gene_id=gene, predicted_cluster=predicted, true_label=actual,
                                       is_medoid=int(i in prediction["medoids"]))
                                   for i, (gene, predicted, actual) in enumerate(zip(ids, prediction["labels"], truth))]
                write_csv(output / prediction_relative, ["gene_id", "predicted_cluster", "true_label", "is_medoid"], prediction_rows)
                contingency = Counter(zip(truth, prediction["labels"]))
                write_csv(output / prediction_relative.with_name(task + "_composition.csv"),
                          ["true_label", "predicted_cluster", "count"],
                          [dict(true_label=actual, predicted_cluster=predicted, count=count)
                           for (actual, predicted), count in sorted(contingency.items())])
                row = {key: dataset[key] for key in ["dataset_id", "dataset_group", "factor", "level_index", "level", "data_seed"]}
                row.update(method=method, task=task, k=k, n_input_genes=n_input, n_genes=len(ids),
                           retained_fraction=len(ids)/n_input, n_timepoints=dataset["n_timepoints"], **metric_values,
                           distance_seconds=matrix_meta["distance_seconds"], clustering_seconds=prediction["clustering_seconds"],
                           clustering_seed=args.clustering_seed, n_init=args.n_init, max_iter=args.max_iter,
                           input_sha256=dataset["input_sha256"], matrix_sha256=matrix_meta["matrix_sha256"],
                           distance_cache_hit=int(distance_reused), clustering_cache_hit=int(cluster_reused),
                           distance_measurement_id=matrix_meta["measurement_id"], clustering_measurement_id=prediction["measurement_id"],
                           parameters_json=json.dumps(parameter_table[method], sort_keys=True), prediction_file=prediction_relative.as_posix())
                results.append(row)
                existing.add((dataset["dataset_id"], method, task))
            processed += len(TASKS)
            write_csv(results_file, FIELDS, results)
        elapsed = time.perf_counter() - started
        print(f"[{index}/{len(datasets)}] {dataset['dataset_id']}: {processed}/{total} task records; this invocation {elapsed:.1f}s", flush=True)
        write_json(output / "progress.json", dict(dataset=dataset["dataset_id"], datasets_completed=index,
                   datasets_requested=len(datasets), invocation_elapsed_seconds=elapsed, total_result_rows=len(results)))
    selected_ids = {dataset["dataset_id"] for dataset in datasets}
    selected_results = [row for row in results if row["dataset_id"] in selected_ids]
    assert len(selected_results) == total
    verify(output, results)
    print(f"COMPLETE: {total} requested task records; results={results_file}", flush=True)


def verify(output, results):
    seen, groups, condition_seeds = set(), Counter(), {}
    for row in results:
        key = (row["dataset_id"], row["method"], row["task"])
        if key in seen:
            raise ValueError(f"Duplicate result: {key}")
        seen.add(key)
        predictions = read_csv(output / row["prediction_file"])
        if len(predictions) != int(row["n_genes"]) or len({item['gene_id'] for item in predictions}) != len(predictions):
            raise ValueError(f"Invalid saved predictions: {key}")
        actual = clustering_scores([item["true_label"] for item in predictions], [item["predicted_cluster"] for item in predictions])
        if any(not math.isclose(actual[name], float(row[name]), abs_tol=1e-12) for name in ("ari", "nmi")):
            raise ValueError(f"Saved score mismatch: {key}")
        if not -1-1e-10 <= float(row["ari"]) <= 1+1e-10 or not -1e-10 <= float(row["nmi"]) <= 1+1e-10:
            raise ValueError("Score outside expected bounds")
        groups[row["dataset_group"]] += 1
        if row["dataset_group"] == "sweep":
            condition = (row["method"], row["task"], row["factor"], str(row["level"]))
            condition_seeds.setdefault(condition, set()).add(int(row["data_seed"]))
    if groups.get("sweep") == 1400:
        assert len(condition_seeds) == 280
        assert all(seeds == set(range(3100, 3105)) for seeds in condition_seeds.values())
    write_json(output / "verification.json", dict(result_rows=len(results), groups=dict(groups),
               unique_result_keys=len(seen), scores_recomputed_from_saved_predictions=True,
               expected_complete_sweep_rows=1400, complete_sweep=groups.get("sweep") == 1400,
               timing_note="Cache hits reference original measured compute times, never fake zero timings. Repeated inputs are not independent repeats."))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--team-root", type=Path, default=DEFAULT_TEAM)
    parser.add_argument("--output", type=Path, default=HERE / "results")
    parser.add_argument("--groups", nargs="+", choices=["smoke", "baseline", "challenge", "sweep"], default=["baseline", "sweep"])
    parser.add_argument("--clustering-seed", type=int, default=42)
    parser.add_argument("--n-init", type=int, default=5)
    parser.add_argument("--max-iter", type=int, default=100)
    args = parser.parse_args()
    if args.n_init < 1 or args.max_iter < 1:
        parser.error("n-init and max-iter must be positive")
    run(args)


if __name__ == "__main__":
    main()
