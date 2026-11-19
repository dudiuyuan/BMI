"""Regenerate every Member A data product and SVG figure from source files."""
from __future__ import annotations

import json
from pathlib import Path

from audit_gpr import audit
from generate_sweep import build
from make_figures import make
from prepare_real_data import prepare
from preprocessing import preprocess_file
from simulation import write_dataset

ROOT = Path(__file__).resolve().parent


def main():
    real = ROOT / "outputs" / "real"
    simulated = ROOT / "outputs" / "simulated"
    report = {}
    for accession in ("GSE3635", "GSE4987"):
        source = ROOT / accession / f"{accession}_series_matrix"
        report[accession] = prepare(accession, source, real / accession)
    gpr_dir = ROOT.parent / "A" / "data" / "GSE3635_RAW"
    if gpr_dir.exists():
        report["gpr_audit"] = audit(gpr_dir, ROOT / "GSE3635" / "GSE3635_series_matrix",
                                    real / "GSE3635" / "gpr_audit.json")
    examples = (
        ("smoke", dict(genes_per_class=2, timepoints=7, seed=101), "Small interface check"),
        ("baseline", dict(genes_per_class=25, timepoints=13, seed=2026), "Baseline 150-gene example"),
        ("challenge", dict(genes_per_class=25, timepoints=13, seed=2027, noise_sd=0.25,
                           max_shift=0.16, amplitude_sd=0.60, missing_rate=0.31,
                           stretch_sd=0.40), "Difficult integration example, not benchmark result"),
    )
    for name, kwargs, purpose in examples:
        target = simulated / name
        write_dataset(target, purpose=purpose, **kwargs)
        report[name] = preprocess_file(target / "expression.csv", target / "processed")
    report["sweep"] = build(simulated / "sweep")
    report["figures"] = make(ROOT, ROOT / "figures")
    print(json.dumps({key: (len(value) if isinstance(value, list) else value.get("rows_kept", value.get("datasets", "done")))
                      for key, value in report.items()}, indent=2))


if __name__ == "__main__":
    main()
