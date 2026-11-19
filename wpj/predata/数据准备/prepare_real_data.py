"""Audit GEO series matrices and export separate gene-by-time runs."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import statistics
from datetime import date
from pathlib import Path

from preprocessing import preprocess_file

ORF = re.compile(r"Y[A-P][LR]\d{3}[CW](?:-[A-Z])?$")
TIME = re.compile(r"time point\s+(\d+)\s+min", re.IGNORECASE)


def sha256(path: Path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def parse_matrix(path: Path):
    meta, table = {}, []
    in_table = False
    with path.open(newline="", encoding="utf-8") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            cells = next(csv.reader([line], delimiter="\t"))
            if in_table:
                table.append(cells)
            elif cells and cells[0].startswith("!"):
                meta.setdefault(cells[0], []).append(cells[1:])
    if not table or table[0][0] != "ID_REF":
        raise ValueError(f"No GEO ID_REF table in {path}")
    header, rows = table[0], table[1:]
    if any(len(row) != len(header) for row in rows):
        raise ValueError("Inconsistent number of sample values")
    sample_ids = meta["!Sample_geo_accession"][0]
    titles = meta["!Sample_title"][0]
    if header[1:] != sample_ids or len(titles) != len(sample_ids):
        raise ValueError("Sample metadata and table columns disagree")
    times = []
    for title in titles:
        match = TIME.search(title)
        if not match:
            raise ValueError(f"Cannot parse sample time: {title}")
        times.append(int(match.group(1)))
    runs, current = [], []
    for j, time in enumerate(times):
        if current and time <= times[current[-1]]:
            runs.append(current)
            current = []
        current.append(j)
    runs.append(current)
    for run in runs:
        if times[run[0]] != 0 or any(times[b] <= times[a] for a, b in zip(run, run[1:])):
            raise ValueError("Runs must begin at 0 and have increasing times")
    return meta, header, rows, sample_ids, titles, times, runs


def _number(token: str):
    if token.strip().lower() in ("", "na", "nan", "null"):
        return None
    val = float(token)
    if not math.isfinite(val):
        return None
    return val


def _csv(path: Path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(fields)
        writer.writerows(rows)


def prepare(accession: str, source: Path, output_root: Path, *, max_missing: float = 0.4):
    meta, header, rows, sample_ids, titles, times, runs = parse_matrix(source)
    observed_accession = meta["!Series_geo_accession"][0][0]
    if observed_accession != accession:
        raise ValueError(f"Expected {accession}, found {observed_accession}")
    output_root.mkdir(parents=True, exist_ok=True)
    ids = [row[0] for row in rows]
    if len(set(ids)) != len(ids):
        raise ValueError("Duplicate source ID_REF: aggregation rule would be required")
    decisions, kept = [], []
    for row in rows:
        original = row[0]
        canonical = original.upper()
        reason = "kept_orf" if ORF.fullmatch(canonical) else "control_or_noncanonical_id"
        decisions.append((original, canonical if reason == "kept_orf" else "", reason))
        if reason == "kept_orf":
            kept.append((canonical, row[1:]))
    canonical_ids = [x[0] for x in kept]
    if len(set(canonical_ids)) != len(canonical_ids):
        raise ValueError("Canonical ORF mapping would create duplicates")
    _csv(output_root / "feature_decisions.csv", ["source_id_ref", "gene_id", "decision"], decisions)
    summary = {"accession": accession, "source_file": str(source), "source_sha256": sha256(source),
               "processed_date": date.today().isoformat(),
               "geo_url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}",
               "platform": meta["!Series_platform_id"][0][0],
               "value_description": meta["!Sample_data_processing"][0][0],
               "overall_design": meta["!Series_overall_design"][0][0],
               "source_features": len(rows), "orf_features": len(kept),
               "excluded_features": len(rows) - len(kept),
               "samples": len(sample_ids), "runs": len(runs),
               "missing_cells_all_features": sum(_number(x) is None for row in rows for x in row[1:]),
               "missing_cells_orf": sum(_number(x) is None for _, vals in kept for x in vals),
               "run_outputs": []}
    sample_manifest = []
    for run_number, indices in enumerate(runs, 1):
        run_name = f"run{run_number:02d}"
        run_dir = output_root / run_name
        run_dir.mkdir(exist_ok=True)
        time_names = [f"t{times[j]}" for j in indices]
        _csv(run_dir / "genes_unprocessed.csv", ["gene_id", *time_names],
             ([gene_id, *[values[j] for j in indices]] for gene_id, values in kept))
        report = preprocess_file(run_dir / "genes_unprocessed.csv", run_dir / "processed", max_missing=max_missing)
        summary["run_outputs"].append({"run": run_name, "times_min": [times[j] for j in indices],
                                        "sample_ids": [sample_ids[j] for j in indices],
                                        "rows_kept_after_missing_filter": report["rows_kept"],
                                        "imputed_cells": report["imputed_cells"]})
        for j in indices:
            col = [_number(row[j + 1]) for row in rows]
            valid = [x for x in col if x is not None]
            sample_manifest.append((run_name, sample_ids[j], times[j], titles[j], len(valid),
                                    len(col) - len(valid), statistics.mean(valid),
                                    statistics.median(valid), statistics.pstdev(valid), min(valid), max(valid)))
    _csv(output_root / "sample_manifest.csv",
         ["run", "gsm", "time_min", "title", "valid_features", "missing_features", "mean", "median", "sd", "min", "max"],
         sample_manifest)
    # This is a source preview, not a separate analysis dataset.
    _csv(output_root / "data_preview.csv", header, rows[:5])
    (output_root / "data_preview_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    cleaning_log = {
        "accession": accession,
        "feature_rule": "upper-case ID_REF; retain canonical yeast nuclear ORF Y[A-P][LR][0-9]{3}[CW] with optional -letter suffix",
        "source_features": len(rows), "kept_orf": len(kept),
        "excluded_control_or_noncanonical": len(rows) - len(kept),
        "case_normalized_orf": sum(original != canonical and reason == "kept_orf"
                                   for original, canonical, reason in decisions),
        "missingness_threshold": max_missing,
        "imputation": "linear in actual time for internal gaps; nearest observed endpoint for edge gaps",
        "runs": summary["run_outputs"],
        "details": "feature_decisions.csv plus each run's processed/row_log.csv and preprocessing_log.json",
    }
    (output_root / "real_data_cleaning_log.json").write_text(json.dumps(cleaning_log, indent=2) + "\n", encoding="utf-8")
    return summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--accession", required=True, choices=("GSE3635", "GSE4987"))
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--max-missing", type=float, default=0.4)
    a = p.parse_args()
    print(json.dumps(prepare(a.accession, a.source, a.output, max_missing=a.max_missing), indent=2))


if __name__ == "__main__":
    main()
