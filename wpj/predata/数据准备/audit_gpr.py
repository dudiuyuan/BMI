"""Read-only spot-level audit of the 13 archived GSE3635 GenePix files."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
from pathlib import Path

from prepare_real_data import parse_matrix, sha256


def audit(gpr_dir: Path, matrix: Path, output: Path):
    _, _, _, sample_ids, _, times, _ = parse_matrix(matrix)
    sample_time = dict(zip(sample_ids, times))
    records = []
    for path in sorted(gpr_dir.glob("*.gpr.gz")):
        gsm = path.name.split(".")[0]
        if gsm not in sample_time:
            raise ValueError(f"GPR sample absent from series matrix: {gsm}")
        with gzip.open(path, "rt", encoding="latin-1", newline="") as f:
            reader = csv.reader(f, delimiter="\t")
            for header in reader:
                if header and header[0] == "Block" and "Log Ratio" in header:
                    break
            else:
                raise ValueError(f"No GenePix spot table: {path}")
            columns = {name: j for j, name in enumerate(header)}
            count = blank_id = non_orf_id = missing_log_ratio = flagged = 0
            unique_ids = set()
            for row in reader:
                if not row or len(row) < len(header):
                    continue
                count += 1
                gene_id = row[columns["ID"]].strip()
                unique_ids.add(gene_id)
                blank_id += not bool(gene_id)
                non_orf_id += not (gene_id.upper().startswith("Y") and len(gene_id) >= 7)
                try:
                    ratio = float(row[columns["Log Ratio"]])
                    missing_log_ratio += not math.isfinite(ratio)
                except ValueError:
                    missing_log_ratio += 1
                flag = row[columns["Flags"]].strip()
                flagged += flag not in ("", "0")
        records.append({"file": str(path), "gsm": gsm, "time_min": sample_time[gsm],
                        "sha256": sha256(path), "spot_rows": count, "unique_spot_ids": len(unique_ids),
                        "blank_id_spots": blank_id, "non_orf_like_spots": non_orf_id,
                        "missing_or_nonfinite_log_ratio": missing_log_ratio, "nonzero_flags": flagged,
                        "fields": list(columns)})
    if set(x["gsm"] for x in records) != set(sample_ids):
        raise ValueError("GPR directory does not contain exactly the matrix samples")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps({"files": records, "warning": "These are repeated spots, not independent genes. GEO processed VALUE is the analysis input."},
                                 indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return {"files": len(records), "spot_rows_range": [min(x["spot_rows"] for x in records), max(x["spot_rows"] for x in records)]}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gpr-dir", type=Path, required=True)
    p.add_argument("--matrix", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    print(json.dumps(audit(a.gpr_dir, a.matrix, a.output), indent=2))


if __name__ == "__main__":
    main()
