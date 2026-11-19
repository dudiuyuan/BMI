"""Label-blind cleaning of gene-by-time expression CSV files (Member A)."""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path

TIME_RE = re.compile(r"^t(\d+(?:\.\d+)?)$")
MISSING = {"", "na", "nan", "null"}


def read_expression(path: Path):
    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fields = reader.fieldnames or []
        if len(fields) < 4 or fields[0] != "gene_id":
            raise ValueError("Expected gene_id followed by at least three time columns")
        names = fields[1:]
        if len(set(names)) != len(names) or any(not TIME_RE.fullmatch(x) for x in names):
            raise ValueError("Time columns must be unique names such as t0,t10,t20")
        times = [float(TIME_RE.fullmatch(x).group(1)) for x in names]
        if any(b <= a for a, b in zip(times, times[1:])):
            raise ValueError("Time columns must increase strictly")
        seen, rows = set(), []
        for line_number, row in enumerate(reader, 2):
            gene_id = row.get("gene_id", "").strip()
            if not gene_id or gene_id in seen:
                raise ValueError(f"Blank or duplicate gene_id at CSV line {line_number}: {gene_id!r}")
            seen.add(gene_id)
            vals = []
            for name in names:
                token = row.get(name)
                if token is None:
                    raise ValueError(f"Missing column {name} at line {line_number}")
                if token.strip().lower() in MISSING:
                    vals.append(None)
                else:
                    try:
                        val = float(token)
                    except ValueError as e:
                        raise ValueError(f"Invalid number at line {line_number}, {name}: {token}") from e
                    if not math.isfinite(val):
                        raise ValueError(f"Nonfinite number at line {line_number}, {name}")
                    vals.append(val)
            rows.append((gene_id, vals))
    return names, times, rows


def interpolate(times: list[float], values: list[float | None]):
    observed = [i for i, value in enumerate(values) if value is not None]
    if not observed:
        raise ValueError("Cannot interpolate an entirely missing trajectory")
    result = [float(v) if v is not None else None for v in values]
    for j in range(len(values)):
        if result[j] is not None:
            continue
        left = max((i for i in observed if i < j), default=None)
        right = min((i for i in observed if i > j), default=None)
        if left is None:
            result[j] = result[right]  # nearest observed endpoint
        elif right is None:
            result[j] = result[left]
        else:
            weight = (times[j] - times[left]) / (times[right] - times[left])
            result[j] = result[left] + weight * (result[right] - result[left])
    return result


def clean(rows, times, *, max_missing: float = 0.4):
    if not 0 <= max_missing < 1:
        raise ValueError("max_missing must be in [0,1)")
    processed, row_log = [], []
    for gene_id, vals in rows:
        missing = sum(v is None for v in vals)
        fraction = missing / len(vals)
        reason = "all_missing" if missing == len(vals) else ("too_missing" if fraction > max_missing else "kept")
        row_log.append({"gene_id": gene_id, "missing_before": missing,
                        "missing_fraction": format(fraction, ".6g"), "decision": reason,
                        "imputed_cells": missing if reason == "kept" else 0})
        if reason == "kept":
            processed.append((gene_id, interpolate(times, vals)))
    return processed, row_log


def _write_csv(path: Path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(fields)
        w.writerows(rows)


def preprocess_file(input_path: Path, output_dir: Path, *, max_missing: float = 0.4):
    names, times, rows = read_expression(input_path)
    processed, log = clean(rows, times, max_missing=max_missing)
    output_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(output_dir / "row_log.csv", ["gene_id", "missing_before", "missing_fraction", "decision", "imputed_cells"],
               ([r[k] for k in ("gene_id", "missing_before", "missing_fraction", "decision", "imputed_cells")] for r in log))
    if not processed:
        raise ValueError("No genes remain after missingness filtering")
    _write_csv(output_dir / "raw_imputed.csv", ["gene_id", *names],
               ([gene_id, *vals] for gene_id, vals in processed))
    centered, zscored = [], []
    constant_count = 0
    for gene_id, vals in processed:
        mean = sum(vals) / len(vals)
        centered_values = [v - mean for v in vals]
        sd = math.sqrt(sum(v * v for v in centered_values) / len(vals))
        if sd < 1e-12:
            zscores = [0.0] * len(vals)
            constant_count += 1
        else:
            zscores = [v / sd for v in centered_values]
        centered.append([gene_id, *centered_values])
        zscored.append([gene_id, *zscores])
    _write_csv(output_dir / "centered.csv", ["gene_id", *names], centered)
    _write_csv(output_dir / "zscore.csv", ["gene_id", *names], zscored)
    report = {"input": str(input_path), "rows_input": len(rows), "rows_kept": len(processed),
              "rows_too_missing": sum(x["decision"] == "too_missing" for x in log),
              "rows_all_missing": sum(x["decision"] == "all_missing" for x in log),
              "missing_cells_input": sum(x["missing_before"] for x in log),
              "imputed_cells": sum(x["imputed_cells"] for x in log),
              "constant_rows_after_imputation": constant_count, "max_missing": max_missing,
              "time_columns": names, "time_values": times,
              "method": "filter > threshold; linear interpolation on actual time; nearest endpoint; per-gene population SD"}
    (output_dir / "preprocessing_log.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--max-missing", type=float, default=0.4)
    a = p.parse_args()
    print(json.dumps(preprocess_file(a.input, a.output, max_missing=a.max_missing), indent=2))


if __name__ == "__main__":
    main()
