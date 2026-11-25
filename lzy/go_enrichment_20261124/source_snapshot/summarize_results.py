"""Summarize fixed benchmark scores and draw exportable scientific charts.

Uses Python standard library and Pillow; SVGs are editable vector figures.
"""
from __future__ import annotations

import argparse
import csv
import html
import json
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

METHODS = ["standard", "time", "time-dynamics", "diagonal-time-dynamics", "diagonal-time", "euclidean", "pearson"]
LABELS = {"standard": "Standard DTW", "time": "Time DTW", "time-dynamics": "Time + dynamics DTW",
          "diagonal-time-dynamics": "Window + time + dynamics", "diagonal-time": "Window + time DTW",
          "euclidean": "Euclidean", "pearson": "Pearson (1-r)"}
COLORS = {method: color for method, color in zip(METHODS, ["#50545a", "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#332288"])}
FACTORS = ["noise_sd", "max_shift", "amplitude_sd", "missing_rate", "stretch_sd"]
FACTOR_LABELS = {"noise_sd": "Noise SD", "max_shift": "Maximum shift / duration", "amplitude_sd": "Amplitude log SD",
                 "missing_rate": "Missing probability", "stretch_sd": "Stretch log SD"}
METRICS = ["ari", "nmi", "distance_seconds", "clustering_seconds"]


def read_csv(filename):
    with filename.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(filename, rows):
    if not rows:
        raise ValueError(f"No rows to write: {filename}")
    filename.parent.mkdir(parents=True, exist_ok=True)
    with filename.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def mean_sd(values):
    return statistics.mean(values), statistics.stdev(values) if len(values) > 1 else 0.0


def summarize(rows):
    groups = defaultdict(list)
    seen = set()
    for row in rows:
        identity = (row["dataset_id"], row["method"], row["task"])
        if identity in seen:
            raise ValueError(f"Duplicate result {identity}")
        seen.add(identity)
        if row["dataset_group"] == "sweep":
            key = (row["task"], row["method"], row["factor"], float(row["level"]))
            groups[key].append(row)
    summary = []
    for (task, method, factor, level), values in sorted(groups.items()):
        seeds = {int(row["data_seed"]) for row in values}
        if len(values) != 5 or seeds != set(range(3100, 3105)):
            raise ValueError(f"Incomplete five-seed condition: {(task, method, factor, level)}")
        record = dict(task=task, method=method, factor=factor, level_index=int(values[0]["level_index"]),
                      level=level, n_data_seeds=len(values), data_seeds=";".join(map(str, sorted(seeds))))
        for metric in METRICS:
            record[metric + "_mean"], record[metric + "_sd"] = mean_sd([float(row[metric]) for row in values])
        record.update(n_genes_min=min(int(row["n_genes"]) for row in values),
                      n_genes_max=max(int(row["n_genes"]) for row in values))
        summary.append(record)
    if len(summary) != 280:
        raise ValueError(f"Expected 280 complete condition summaries, found {len(summary)}")
    return summary


def paired(rows):
    sweep = [row for row in rows if row["dataset_group"] == "sweep"]
    lookup = {(row["dataset_id"], row["task"], row["method"]): row for row in sweep}
    groups, raw = defaultdict(list), []
    for row in sweep:
        if row["method"] == "standard":
            continue
        control = lookup[row["dataset_id"], row["task"], "standard"]
        record = {key: row[key] for key in ["dataset_id", "task", "method", "factor", "level", "data_seed"]}
        record.update(ari_delta=float(row["ari"])-float(control["ari"]),
                      nmi_delta=float(row["nmi"])-float(control["nmi"]))
        raw.append(record)
        groups[row["task"], row["method"], row["factor"], float(row["level"])].append(record)
    aggregate = []
    for (task, method, factor, level), values in sorted(groups.items()):
        record = dict(task=task, method=method, factor=factor, level=level, n_pairs=len(values))
        for metric in ("ari", "nmi"):
            deltas = [row[metric + "_delta"] for row in values]
            record[metric + "_delta_mean"], record[metric + "_delta_sd"] = mean_sd(deltas)
            record[metric + "_wins"] = sum(value > 1e-12 for value in deltas)
            record[metric + "_ties"] = sum(abs(value) <= 1e-12 for value in deltas)
            record[metric + "_losses"] = sum(value < -1e-12 for value in deltas)
        aggregate.append(record)
    return raw, aggregate


class Chart:
    """Draw identical geometry to a high-resolution raster and a vector SVG."""
    def __init__(self, width, height):
        self.width, self.height = width, height
        self.image = Image.new("RGB", (width, height), "white")
        self.draw = ImageDraw.Draw(self.image)
        self.svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
                    '<rect width="100%" height="100%" fill="white"/>']
        self.fonts = {}

    def font(self, size):
        if size not in self.fonts:
            for filename in ["C:/Windows/Fonts/arial.ttf", "C:/Windows/Fonts/calibri.ttf", "DejaVuSans.ttf"]:
                try:
                    self.fonts[size] = ImageFont.truetype(filename, size)
                    break
                except OSError:
                    continue
            else:
                self.fonts[size] = ImageFont.load_default()
        return self.fonts[size]

    def text(self, x, y, text, size=19, color="#222222", anchor="left"):
        label = str(text)
        font = self.font(size)
        box = self.draw.textbbox((0, 0), label, font=font)
        width = box[2] - box[0]
        xx = x - width/2 if anchor == "middle" else x-width if anchor == "right" else x
        self.draw.text((xx, y), label, font=font, fill=color)
        svg_anchor = {"middle": "middle", "right": "end", "left": "start"}[anchor]
        self.svg.append(f'<text x="{x}" y="{y+size}" font-family="Arial, sans-serif" font-size="{size}" text-anchor="{svg_anchor}" fill="{color}">{html.escape(label)}</text>')

    def line(self, points, color="#222222", width=2, dashed=False):
        points = [(float(x), float(y)) for x, y in points]
        if dashed:
            for (x1, y1), (x2, y2) in zip(points, points[1:]):
                length = math.hypot(x2-x1, y2-y1)
                for offset in range(0, int(length)+1, 15):
                    low, high = offset/max(length, 1), min(offset+8, length)/max(length, 1)
                    self.draw.line((x1+(x2-x1)*low, y1+(y2-y1)*low,
                                    x1+(x2-x1)*high, y1+(y2-y1)*high), fill=color, width=width)
        else:
            self.draw.line(points, fill=color, width=width)
        dash = ' stroke-dasharray="8,7"' if dashed else ''
        self.svg.append(f'<polyline points="{" ".join(f"{x:.2f},{y:.2f}" for x,y in points)}" fill="none" stroke="{color}" stroke-width="{width}"{dash}/>')

    def marker(self, x, y, color, shape):
        radius = 5
        if shape % 3 == 0:
            self.draw.ellipse((x-radius, y-radius, x+radius, y+radius), fill="white", outline=color, width=3)
            self.svg.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{radius}" fill="white" stroke="{color}" stroke-width="3"/>')
        elif shape % 3 == 1:
            self.draw.rectangle((x-radius, y-radius, x+radius, y+radius), fill="white", outline=color, width=3)
            self.svg.append(f'<rect x="{x-radius:.2f}" y="{y-radius:.2f}" width="10" height="10" fill="white" stroke="{color}" stroke-width="3"/>')
        else:
            points = [(x, y-radius-1), (x+radius+1, y+radius), (x-radius-1, y+radius), (x, y-radius-1)]
            self.draw.polygon(points, fill="white", outline=color, width=3)
            self.svg.append(f'<polygon points="{" ".join(f"{a:.2f},{b:.2f}" for a,b in points)}" fill="white" stroke="{color}" stroke-width="3"/>')

    def save(self, path):
        self.image.save(path.with_suffix(".png"), dpi=(220, 220))
        path.with_suffix(".svg").write_text("\n".join(self.svg + ["</svg>"]), encoding="utf-8")


def plot(summary, output, task, metric):
    selected = [row for row in summary if row["task"] == task]
    timing = metric.endswith("seconds")
    title_metric = {"ari": "Adjusted Rand index (ARI)", "nmi": "Normalized mutual information (NMI)",
                    "distance_seconds": "Distance computation time", "clustering_seconds": "Clustering computation time"}[metric]
    chart = Chart(1980, 1400)
    chart.text(55, 24, f"{title_metric} - {'phase-sensitive, K=6' if task=='phase6' else 'shape-based, K=5'}", 32)
    chart.text(55, 73, "Mean +/- sample SD across five simulated datasets per condition; fixed clustering seed=42, n_init=5", 20, "#555555")
    all_lows = [row[metric+"_mean"]-row[metric+"_sd"] for row in selected]
    all_highs = [row[metric+"_mean"]+row[metric+"_sd"] for row in selected]
    if timing:
        positive = [value for value in all_lows if value > 0]
        ymin = 10 ** math.floor(math.log10(min(positive)))
        ymax = 10 ** math.ceil(math.log10(max(all_highs)))
        if ymax <= ymin:
            ymax = ymin * 10
        y_ticks = [10**power for power in range(int(round(math.log10(ymin))), int(round(math.log10(ymax)))+1)]
    else:
        ymin, ymax = min(-0.05, min(all_lows)-0.03), max(1.05, max(all_highs)+0.03)
        ymin = math.floor(ymin*10)/10
        y_ticks = [value/5 for value in range(math.ceil(ymin*5), 6)]
    panel_w, panel_h = 940, 380
    positions = [(45, 140), (1025, 140), (45, 540), (1025, 540), (45, 940)]
    for factor, (left, top) in zip(FACTORS, positions):
        chart.text(left+65, top, FACTOR_LABELS[factor], 24)
        x0, x1, y0, y1 = left+80, left+panel_w-35, top+48, top+panel_h-62
        rows = [row for row in selected if row["factor"] == factor]
        levels = sorted({row["level"] for row in rows})
        pad = (levels[-1]-levels[0])*.045
        xmin, xmax = levels[0]-pad, levels[-1]+pad
        xx = lambda value: x0+(value-xmin)/(xmax-xmin)*(x1-x0)
        if timing:
            yy = lambda value: y1-(math.log10(max(value,ymin))-math.log10(ymin))/(math.log10(ymax)-math.log10(ymin))*(y1-y0)
        else:
            yy = lambda value: y1-(value-ymin)/(ymax-ymin)*(y1-y0)
        for tick in y_ticks:
            chart.line([(x0, yy(tick)), (x1, yy(tick))], "#e2e4e7", 1)
            chart.text(x0-12, yy(tick)-12, f"{tick:g}", 17, anchor="right")
        chart.line([(x0, y0), (x0, y1), (x1, y1)], width=2)
        for tick in levels:
            chart.line([(xx(tick), y1), (xx(tick), y1+5)], width=2)
            chart.text(xx(tick), y1+10, f"{tick:g}", 18, anchor="middle")
        chart.text(x0, y1+37, "Perturbation strength", 17, "#555555")
        for index, method in enumerate(METHODS):
            series = sorted([row for row in rows if row["method"] == method], key=lambda row: row["level"])
            color = COLORS[method]
            points = [(xx(row["level"]), yy(row[metric+"_mean"])) for row in series]
            chart.line(points, color, 3, dashed=method in ["euclidean", "pearson"])
            for row, (x,y) in zip(series, points):
                low = max(ymin, row[metric+"_mean"]-row[metric+"_sd"])
                high = min(ymax, row[metric+"_mean"]+row[metric+"_sd"])
                chart.line([(x, yy(low)), (x, yy(high))], color, 2)
                chart.line([(x-4, yy(low)), (x+4, yy(low))], color, 2)
                chart.line([(x-4, yy(high)), (x+4, yy(high))], color, 2)
                chart.marker(x, y, color, index)
    legend_x, legend_y = 1100, 963
    chart.text(legend_x, legend_y-20, "Methods", 24)
    for index, method in enumerate(METHODS):
        y = legend_y+25+index*34
        chart.line([(legend_x, y+12), (legend_x+55, y+12)], COLORS[method], 3, method in ["euclidean", "pearson"])
        chart.marker(legend_x+27, y+12, COLORS[method], index)
        chart.text(legend_x+72, y, LABELS[method], 20)
    chart.text(legend_x, 1260, "Y: seconds (log scale)" if timing else "Y: score; higher is better", 18, "#555555")
    chart.text(55, 1360, "Error bars show dataset variation, not confidence intervals. Timing excludes file I/O; identical inputs reuse measured times.", 17, "#555555")
    chart.save(output / "figures" / f"{task}_{metric}")


def report(rows, summary, paired_summary, output):
    lookup = {(row["task"], row["method"], row["factor"], row["level"]): row for row in summary}
    max_levels = {factor: max(row["level"] for row in summary if row["factor"] == factor) for factor in FACTORS}
    baseline = [row for row in rows if row["dataset_group"] == "baseline"]
    lines = ["# Project 11：第一版正式模拟 benchmark 结果", "",
             f"共 {len(rows)} 行评分：baseline {len(baseline)} 行、正式 sweep {sum(row['dataset_group']=='sweep' for row in rows)} 行。",
             "使用 zrz 最终版七种距离、xqf 交替更新 K-medoids；直接使用 centered，聚类 seed=42、n_init=5、max_iter=100。",
             "ARI/NMI按保留基因ID关联后计算。以下 ± 为五个数据种子的样本标准差，不是置信区间。", "",
             "## 主要可观察结果", ""]
    for factor in FACTORS:
        high = max_levels[factor]
        standard = lookup['phase6', 'standard', factor, high]
        candidates = [lookup['phase6', method, factor, high] for method in METHODS]
        best = max(candidates, key=lambda row: row['ari_mean'])
        lines.append(f"- 六类任务，{factor}={high:g}：标准DTW ARI={standard['ari_mean']:.3f}±{standard['ari_sd']:.3f}；"
                     f"本条件七种方法中观察到的最高平均值为 {best['method']}，{best['ari_mean']:.3f}±{best['ari_sd']:.3f}。"
                     "这仅描述该条件，不等于通用最优或统计显著。")
    shape_standard = lookup['shape5', 'standard', 'stretch_sd', max_levels['stretch_sd']]
    shape_modifications = [lookup['shape5', method, 'stretch_sd', max_levels['stretch_sd']] for method in METHODS[1:5]]
    shape_best_modification = max(shape_modifications, key=lambda row: row['ari_mean'])
    lines += ["", "两种标签目标呈现不同取舍：例如五类形状任务的最高拉伸条件，"
              f"标准DTW ARI={shape_standard['ari_mean']:.3f}±{shape_standard['ari_sd']:.3f}，"
              f"四个修改配置中最高平均值为{shape_best_modification['method']}，"
              f"{shape_best_modification['ari_mean']:.3f}±{shape_best_modification['ari_sd']:.3f}。"
              "这与惩罚保留时间相位、但可能限制形状对齐的假设相容；单凭分数不能证明具体错误原因，需结合composition表核查。",
              "高幅度差异时Pearson更稳健，与其相关距离弱化幅度尺度的定义相容；它也因此可能忽略具有生物学意义的表达强度差异。"]
    for task, description in [('phase6', '六类相位任务（K=6）'), ('shape5', '五类形状任务（K=5）')]:
        lines += ["", f"## {description}：各因素最高强度下的 ARI", "",
                  "| 方法 | 噪声0.25 | 偏移0.16 | 幅度logSD0.6 | 缺失0.31 | 伸缩logSD0.4 |",
                  "| --- | ---: | ---: | ---: | ---: | ---: |"]
        for method in METHODS:
            cells = [f"{lookup[task,method,factor,max_levels[factor]]['ari_mean']:.3f} ± {lookup[task,method,factor,max_levels[factor]]['ari_sd']:.3f}" for factor in FACTORS]
            lines.append("| " + " | ".join([method, *cells]) + " |")
    lines += ["", "## 与标准 DTW 的配对比较", "",
              "每个差值来自相同数据seed、同一任务：该方法分数－standard分数。完整逐次差值见paired_runs.csv，汇总见paired_vs_standard.csv。", "",
              "| 六类任务：最高噪声0.25 | ARI差值均值 ± SD | 五个seed胜／平／负 |", "| --- | ---: | ---: |"]
    for method in METHODS[1:]:
        row = next(row for row in paired_summary if row['task']=='phase6' and row['method']==method and row['factor']=='noise_sd' and row['level']==max_levels['noise_sd'])
        lines.append(f"| {method} | {row['ari_delta_mean']:.3f} ± {row['ari_delta_sd']:.3f} | {row['ari_wins']}/{row['ari_ties']}/{row['ari_losses']} |")
    lines += ["", "## baseline：开发数据结果，单独报告", "", "baseline参与过选参，不是独立最终测试；固定一个聚类seed的结果与此前20聚类seed平均值不必相同。", "",
              "| 方法 | 六类ARI | 六类NMI | 五类ARI | 五类NMI |", "| --- | ---: | ---: | ---: | ---: |"]
    for method in METHODS:
        pair = {row['task']: row for row in baseline if row['method']==method}
        if len(pair) == 2:
            lines.append(f"| {method} | {float(pair['phase6']['ari']):.3f} | {float(pair['phase6']['nmi']):.3f} | {float(pair['shape5']['ari']):.3f} | {float(pair['shape5']['nmi']):.3f} |")
    missing = [row for row in rows if row['dataset_group']=='sweep' and row['factor']=='missing_rate' and float(row['level'])==max_levels['missing_rate']]
    lines += ["", "## 运行时间、过滤与解释范围", "",
              f"- 最高缺失条件保留 {min(int(row['n_genes']) for row in missing)}–{max(int(row['n_genes']) for row in missing)} / 150 个基因。评分是在各自保留集上完成，不能忽略样本变化。",
              "- 图表中的距离时间和聚类时间分别记录，使用对数纵轴。两个任务复用同一距离，不能把它们相加。重复基础输入沿用原测量时间，有measurement_id和cache标记；复用不等于算法耗时0。",
              "- 每个条件有五套独立模拟，只有一个固定外部聚类seed；内部五起点不是五套数据。这一版没有独立测量聚类初始化敏感性。",
              "- 100个目录含重复基础条件，只有80份不同矩阵。不能将全部成绩平均成总体准确率，也不能当作100个独立重复。",
              "- 这批参数按六类标签选择，五类任务是单列补充。它合并早晚峰，broad_peak清洗前占50/150，其他类别各25。",
              "- 某模式同时改变多个参数；成绩说明整个配置的效果，不能仅凭模式间分差断言单个动态项或窗口的贡献。",
              "- 平坦/噪声类是零模板加噪声，没有共同动态形状。高噪声会造成真实类别重叠。",
              "- 耗时是本机一次计算测量；不作细微速度差的稳定性或显著性断言。",
              "- 未做N/T规模扩展、真实GEO生物学验证或新参数搜索。真实数据与模拟ARI/NMI不能混为同一评价。", "",
              "## 怎样检查错误", "",
              "到 ../predictions/<dataset_id>/<method>/ 查看phase6或shape5的预测CSV和composition交叉计数。关注早晚峰合并、同一类型被拆分、平坦类混入其他类；相同ARI可能对应不同错误。", "",
              "## 主要文件", "", "- summary.csv：280个条件的均值、样本SD及保留数。",
              "- paired_runs.csv：相同数据seed与标准DTW的逐次分差。", "- paired_vs_standard.csv：240个条件的配对分差汇总。",
              "- figures/：六张PNG与对应SVG，包含两任务ARI/NMI和六类任务距离/聚类耗时。", "",
              "结果解读由程序从实测评分生成，用于组内分析；正式提交前请自行审阅并按课程要求披露AI使用。"]
    (output / "结果解读.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = read_csv(args.results)
    summary = summarize(rows)
    pairs, paired_summary = paired(rows)
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "figures").mkdir(exist_ok=True)
    write_csv(args.output / "summary.csv", summary)
    write_csv(args.output / "paired_runs.csv", pairs)
    write_csv(args.output / "paired_vs_standard.csv", paired_summary)
    write_csv(args.output / "baseline_scores.csv", [row for row in rows if row["dataset_group"] == "baseline"])
    for task in ['phase6', 'shape5']:
        for metric in ['ari', 'nmi']:
            plot(summary, args.output, task, metric)
    for metric in ['distance_seconds', 'clustering_seconds']:
        plot(summary, args.output, 'phase6', metric)
    report(rows, summary, paired_summary, args.output)
    (args.output / "summary_verification.json").write_text(json.dumps(dict(
        input_rows=len(rows), condition_summaries=len(summary), paired_runs=len(pairs),
        paired_condition_summaries=len(paired_summary), independent_data_seeds_per_condition=5,
        uncertainty="sample SD (ddof=1), not confidence interval", plots=6), indent=2), encoding="utf-8")
    print(json.dumps(dict(output=str(args.output), rows=len(rows), summaries=len(summary), plots=6), ensure_ascii=False))


if __name__ == '__main__':
    main()
