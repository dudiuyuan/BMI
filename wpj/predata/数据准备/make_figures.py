"""Make data-driven SVG QC and simulation figures without plotting dependencies."""
from __future__ import annotations

import argparse
import csv
import html
import json
from pathlib import Path

from simulation import CLASSES, template

COLORS = ("#146c94", "#c7522a", "#227a53", "#824b9e", "#a07812", "#353a42")


def _read(path: Path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _num(value):
    return None if value is None or value == "" else float(value)


def plot(path: Path, title: str, panels: list[dict], *, y_label: str):
    cols = 2
    cell_w, cell_h = 520, 320
    width, height = 1100, 90 + ((len(panels) + cols - 1) // cols) * cell_h
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
           '<rect width="100%" height="100%" fill="white"/>',
           f'<text x="35" y="40" font-family="Arial,sans-serif" font-size="23" font-weight="bold">{html.escape(title)}</text>']
    for index, panel in enumerate(panels):
        x0 = 65 + (index % cols) * 540
        y0 = 90 + (index // cols) * cell_h
        left, top, w, h = x0 + 50, y0 + 30, 410, 195
        all_points = [(x, y) for line in panel["lines"] for x, y in line[1] if y is not None]
        all_x = [x for line in panel["lines"] for x, _ in line[1]]
        if not all_points:
            continue
        min_x, max_x = min(all_x), max(all_x)
        min_y, max_y = min(y for _, y in all_points), max(y for _, y in all_points)
        pad = max(0.05, (max_y - min_y) * 0.12)
        min_y -= pad
        max_y += pad
        if max_x == min_x:
            max_x += 1
        X = lambda x: left + w * (x - min_x) / (max_x - min_x)
        Y = lambda y: top + h * (max_y - y) / (max_y - min_y)
        out.append(f'<text x="{left}" y="{y0 + 18}" font-family="Arial,sans-serif" font-size="17" font-weight="bold">{html.escape(panel["title"])}</text>')
        for fraction in (0, 0.5, 1):
            gy = top + h * fraction
            label = max_y - (max_y - min_y) * fraction
            out.append(f'<line x1="{left}" y1="{gy:.1f}" x2="{left + w}" y2="{gy:.1f}" stroke="#e1e5ea"/>')
            out.append(f'<text x="{left - 8}" y="{gy + 4:.1f}" text-anchor="end" font-family="Arial,sans-serif" font-size="11">{label:.2f}</text>')
        out.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + h}" stroke="#48515a"/>')
        out.append(f'<line x1="{left}" y1="{top + h}" x2="{left + w}" y2="{top + h}" stroke="#48515a"/>')
        for tick in (min_x, (min_x + max_x) / 2, max_x):
            out.append(f'<text x="{X(tick):.1f}" y="{top + h + 18}" text-anchor="middle" font-family="Arial,sans-serif" font-size="11">{tick:g}</text>')
        for j, (label, points) in enumerate(panel["lines"]):
            color = COLORS[j % len(COLORS)]
            segment = []
            for x, y in points + [(max_x, None)]:
                if y is None:
                    if len(segment) > 1:
                        coords = " ".join(f"{X(px):.1f},{Y(py):.1f}" for px, py in segment)
                        out.append(f'<polyline points="{coords}" fill="none" stroke="{color}" stroke-width="2.2"/>')
                    segment = []
                else:
                    segment.append((x, y))
                    out.append(f'<circle cx="{X(x):.1f}" cy="{Y(y):.1f}" r="2.4" fill="{color}"/>')
            legend_x = left + (j % 2) * 205
            legend_y = top + h + 40 + (j // 2) * 17
            out.append(f'<line x1="{legend_x}" y1="{legend_y - 4}" x2="{legend_x + 18}" y2="{legend_y - 4}" stroke="{color}" stroke-width="3"/>')
            out.append(f'<text x="{legend_x + 24}" y="{legend_y}" font-family="Arial,sans-serif" font-size="11">{html.escape(label)}</text>')
        missing_x = sorted({x for _, points in panel["lines"] for x, y in points if y is None})
        for x in missing_x:
            px, py = X(x), top + h - 8
            out.append(f'<path d="M {px - 5:.1f},{py - 5:.1f} L {px + 5:.1f},{py + 5:.1f} M {px - 5:.1f},{py + 5:.1f} L {px + 5:.1f},{py - 5:.1f}" stroke="#c72b35" stroke-width="2"/>')
        out.append(f'<text x="{left + w / 2}" y="{top + h + 75}" text-anchor="middle" font-family="Arial,sans-serif" font-size="12">Time after release (min)</text>')
        out.append(f'<text x="{left - 50}" y="{top + h / 2}" transform="rotate(-90 {left - 50} {top + h / 2})" text-anchor="middle" font-family="Arial,sans-serif" font-size="12">{html.escape(y_label)}</text>')
    out.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(out) + "\n", encoding="utf-8")


def make(root: Path, output: Path, *, seed: int = 3100):
    output.mkdir(parents=True, exist_ok=True)
    times = [i * 2 for i in range(61)]
    panels = [{"title": kind.replace("_", " "), "lines": [(kind, [(t, template(kind, t / 120)) for t in times])]}
              for kind in CLASSES]
    plot(output / "six_templates.svg", "Six unperturbed trajectory families", panels,
         y_label="Dimensionless template value")
    sweep = root / "outputs" / "simulated" / "sweep"
    design = json.loads((sweep / "sweep_design.json").read_text(encoding="utf-8"))
    for factor, levels in design["factors"].items():
        lines = []
        for level_index, level in enumerate(levels):
            path = sweep / factor / f"level{level_index}_{level:g}" / f"seed{seed}" / "expression.csv"
            rows = _read(path)
            row = next(r for r in rows if r["gene_id"] == "SIM_early_peak_0001")
            points = [(float(name[1:]), _num(row[name])) for name in list(row)[1:]]
            lines.append((f"{factor}={level:g}", points))
        if factor == "missing_rate":
            panels = [{"title": f"p={level:g}; red X = missing", "lines": [lines[i]]}
                      for i, level in enumerate(levels)]
        else:
            panels = [{"title": "SIM_early_peak_0001", "lines": lines}]
        plot(output / f"effect_{factor}.svg", f"Effect of {factor} on one fixed gene (seed {seed})",
             panels, y_label="Simulated expression (a.u.)")
    for accession in ("GSE3635", "GSE4987"):
        real = root / "outputs" / "real" / accession
        manifest = _read(real / "sample_manifest.csv")
        panels = []
        for run in sorted({r["run"] for r in manifest}):
            records = [r for r in manifest if r["run"] == run]
            panels.append({"title": run,
                           "lines": [("median", [(float(r["time_min"]), float(r["median"])) for r in records]),
                                     ("SD", [(float(r["time_min"]), float(r["sd"])) for r in records])]})
        plot(output / f"qc_{accession}.svg", f"{accession}: sample distributions", panels,
             y_label="Normalized log-ratio")
    rows = _read(root / "outputs" / "real" / "GSE3635" / "run01" / "genes_unprocessed.csv")
    gene_rows = {r["gene_id"]: r for r in rows}
    lines = []
    for symbol, gene_id in (("MCM3", "YEL032W"), ("CLN2", "YPL256C"),
                            ("CLB2", "YPR119W"), ("HTA1", "YDR225W")):
        row = gene_rows[gene_id]
        lines.append((f"{symbol} ({gene_id})", [(float(name[1:]), _num(row[name])) for name in list(row)[1:]]))
    plot(output / "gse3635_example_genes.svg", "GSE3635: four traceable genes", [{"title": "Processed GEO values", "lines": lines}],
         y_label="Normalized log-ratio")
    return sorted(p.name for p in output.glob("*.svg"))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=Path(__file__).resolve().parent)
    p.add_argument("--output", type=Path)
    p.add_argument("--seed", type=int, default=3100)
    a = p.parse_args()
    print("\n".join(make(a.root, a.output or a.root / "figures", seed=a.seed)))


if __name__ == "__main__":
    main()
