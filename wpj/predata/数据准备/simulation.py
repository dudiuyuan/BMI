"""Reproducible, label-separated gene trajectory simulator (Member A)."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
from pathlib import Path

CLASSES = (
    "early_peak", "late_peak", "oscillatory", "transient_induction",
    "transient_repression", "flat_noise",
)


def _rng(seed: int, gene_id: str, stream: str) -> random.Random:
    digest = hashlib.sha256(f"{seed}|{gene_id}|{stream}".encode()).digest()
    return random.Random(int.from_bytes(digest[:16], "big"))


def template(kind: str, u: float, width: float = 0.085) -> float:
    """Dimensionless family curve evaluated at normalized time u in [0, 1]."""
    if kind in ("early_peak", "late_peak"):
        center = 0.20 if kind == "early_peak" else 0.78
        return math.exp(-0.5 * ((u - center) / width) ** 2)
    if kind == "oscillatory":
        return math.sin(4 * math.pi * u)
    if kind in ("transient_induction", "transient_repression"):
        # A finite central plateau distinguishes this pulse from a Gaussian peak.
        rise = 1 / (1 + math.exp(-min(60.0, max(-60.0, (u - 0.39) / 0.025))))
        fall = 1 / (1 + math.exp(-min(60.0, max(-60.0, (0.64 - u) / 0.025))))
        pulse = rise * fall
        return pulse if kind == "transient_induction" else -pulse
    if kind == "flat_noise":
        return 0.0
    raise ValueError(f"Unknown class: {kind}")


def generate(*, genes_per_class: int = 25, timepoints: int = 13, seed: int = 2026,
             noise_sd: float = 0.05, max_shift: float = 0.03,
             amplitude_sd: float = 0.15, missing_rate: float = 0.0,
             stretch_sd: float = 0.08, duration: float = 120.0):
    if genes_per_class < 1 or timepoints < 3 or seed < 0 or duration <= 0:
        raise ValueError("Require genes_per_class>=1, timepoints>=3, seed>=0, duration>0")
    if min(noise_sd, max_shift, amplitude_sd, missing_rate, stretch_sd) < 0 or missing_rate >= 1:
        raise ValueError("Perturbations must be nonnegative and missing_rate < 1")
    if max_shift > 1:
        raise ValueError("max_shift is a fraction of total duration and must be <= 1")
    times = [duration * j / (timepoints - 1) for j in range(timepoints)]
    columns = [f"t{t:g}" for t in times]
    expression, truth, parameters = [], [], []
    for kind in CLASSES:
        for index in range(genes_per_class):
            gene_id = f"SIM_{kind}_{index + 1:04d}"
            morphology = _rng(seed, gene_id, "morphology")
            perturb = _rng(seed, gene_id, "perturbation")
            errors = _rng(seed, gene_id, "measurement")
            masks = _rng(seed, gene_id, "missingness")
            width = 0.075 + 0.020 * morphology.random()
            baseline = morphology.uniform(-0.1, 0.1)
            amplitude_z = perturb.gauss(0, 1)
            shift_z = perturb.uniform(-1, 1)
            stretch_z = perturb.gauss(0, 1)
            amplitude = math.exp(amplitude_sd * amplitude_z)
            shift = max_shift * shift_z
            gamma = math.exp(stretch_sd * stretch_z)
            values = []
            for j in range(timepoints):
                u = j / (timepoints - 1)
                warped = min(1.0, max(0.0, u ** gamma - shift))
                value = baseline + amplitude * template(kind, warped, width)
                value += noise_sd * errors.gauss(0, 1)
                values.append("" if masks.random() < missing_rate else format(value, ".12g"))
            expression.append(dict(gene_id=gene_id, **dict(zip(columns, values))))
            truth.append({"gene_id": gene_id, "true_cluster": kind,
                          "shape_cluster": "broad_peak" if kind in ("early_peak", "late_peak") else kind})
            parameters.append({"gene_id": gene_id, "true_cluster": kind,
                               "width": width, "baseline": baseline,
                               "amplitude_z": amplitude_z, "amplitude": amplitude,
                               "shift_z": shift_z, "shift": shift,
                               "stretch_z": stretch_z, "gamma": gamma,
                               "missing_count": values.count("")})
    config = dict(genes_per_class=genes_per_class, timepoints=timepoints, seed=seed,
                  noise_sd=noise_sd, max_shift=max_shift, amplitude_sd=amplitude_sd,
                  missing_rate=missing_rate, stretch_sd=stretch_sd, duration=duration)
    return columns, expression, truth, parameters, config


def _write_csv(path: Path, fields: list[str], rows: list[dict]):
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_dataset(output: Path, *, purpose: str = "benchmark input", **kwargs):
    columns, expression, truth, parameters, config = generate(**kwargs)
    output.mkdir(parents=True, exist_ok=True)
    _write_csv(output / "expression.csv", ["gene_id", *columns], expression)
    _write_csv(output / "ground_truth.csv", ["gene_id", "true_cluster", "shape_cluster"], truth)
    _write_csv(output / "gene_parameters.csv", list(parameters[0]), parameters)
    metadata = {"generator": "simulation.py", "model_version": 1, "purpose": purpose,
                "config": config, "classes": list(CLASSES), "rows": len(expression),
                "time_columns": columns, "missing_cells": sum(p["missing_count"] for p in parameters),
                "primary_truth": "true_cluster: six phase-sensitive classes",
                "secondary_truth": "shape_cluster: early/late peaks merged (five classes)"}
    (output / "metadata.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return metadata


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--purpose", default="benchmark input")
    p.add_argument("--genes-per-class", type=int, default=25)
    p.add_argument("--timepoints", type=int, default=13)
    p.add_argument("--seed", type=int, default=2026)
    p.add_argument("--noise-sd", type=float, default=0.05)
    p.add_argument("--max-shift", type=float, default=0.03)
    p.add_argument("--amplitude-sd", type=float, default=0.15)
    p.add_argument("--missing-rate", type=float, default=0.0)
    p.add_argument("--stretch-sd", type=float, default=0.08)
    p.add_argument("--duration", type=float, default=120.0)
    a = vars(p.parse_args())
    output = a.pop("output")
    print(json.dumps(write_dataset(output, **a), indent=2))


if __name__ == "__main__":
    main()
