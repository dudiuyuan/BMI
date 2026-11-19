"""Create paired, one-factor-at-a-time benchmark inputs (no clustering)."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

from preprocessing import preprocess_file
from simulation import write_dataset

LEVELS = {
    "noise_sd": (0.0, 0.05, 0.12, 0.25),
    "max_shift": (0.0, 0.03, 0.08, 0.16),
    "amplitude_sd": (0.0, 0.15, 0.35, 0.60),
    "missing_rate": (0.0, 0.08, 0.16, 0.31),
    "stretch_sd": (0.0, 0.08, 0.20, 0.40),
}
BASE = {"noise_sd": 0.05, "max_shift": 0.03, "amplitude_sd": 0.15,
        "missing_rate": 0.0, "stretch_sd": 0.08}


def digest(path: Path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build(output: Path, *, genes_per_class: int = 25, timepoints: int = 13,
          first_seed: int = 3100, seeds: int = 5, max_missing: float = 0.4):
    if seeds < 1:
        raise ValueError("seeds must be positive")
    output.mkdir(parents=True, exist_ok=True)
    entries = []
    for factor, levels in LEVELS.items():
        for level_index, level in enumerate(levels):
            for seed in range(first_seed, first_seed + seeds):
                directory = output / factor / f"level{level_index}_{level:g}" / f"seed{seed}"
                config = dict(BASE, genes_per_class=genes_per_class, timepoints=timepoints, seed=seed)
                config[factor] = level
                meta = write_dataset(directory, purpose=f"paired sweep: {factor} level {level_index}", **config)
                prep = preprocess_file(directory / "expression.csv", directory / "processed", max_missing=max_missing)
                entries.append({"factor": factor, "level_index": level_index, "level": level,
                                "seed": seed, "relative_dir": str(directory.relative_to(output)),
                                "genes_input": meta["rows"], "genes_kept": prep["rows_kept"],
                                "missing_cells": meta["missing_cells"],
                                "expression_sha256": digest(directory / "expression.csv"),
                                "centered_sha256": digest(directory / "processed" / "centered.csv")})
    fields = list(entries[0])
    with (output / "sweep_index.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(entries)
    design = {"factors": {key: list(value) for key, value in LEVELS.items()},
              "fixed_others": BASE, "genes_per_class": genes_per_class,
              "timepoints": timepoints, "seeds": list(range(first_seed, first_seed + seeds)),
              "max_missing": max_missing, "datasets": len(entries),
              "paired_rule": "For a fixed seed, gene IDs and independent random streams are stable across levels. Non-target parameters and latent draws stay identical."}
    (output / "sweep_design.json").write_text(json.dumps(design, indent=2) + "\n", encoding="utf-8")
    return design


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--genes-per-class", type=int, default=25)
    p.add_argument("--timepoints", type=int, default=13)
    p.add_argument("--first-seed", type=int, default=3100)
    p.add_argument("--seeds", type=int, default=5)
    p.add_argument("--max-missing", type=float, default=0.4)
    a = p.parse_args()
    print(json.dumps(build(**vars(a)), indent=2))


if __name__ == "__main__":
    main()
