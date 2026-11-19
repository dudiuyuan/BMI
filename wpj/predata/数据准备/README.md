# BMI3 ICA Project 11：A 同学的数据准备

这个目录保存 A 同学的数据准备工作：下载原件、生成和清洗代码、给组员使用的输入，以及分享材料。**这里还没有 DTW、聚类或 ARI/NMI 结果。**

## 从哪里开始

| 你想做什么 | 先看这个文件 |
| --- | --- |
| 给组员介绍做了什么 | [6 页 PPT](reports/A_data_preparation_share.pptx) 和 [逐页讲解稿](分享讲解稿.md) |
| 快速了解数据和实测结果 | [数据准备报告](reports/data_preparation_report.md) |
| 理解六类曲线、五种扰动和清洗理由 | [方法说明](simulation_design.md) |
| 接手代码或重跑数据 | 从下方的“快速开始”和“交接给 B/C/D/E”阅读 |

目录分为四块：`GSE3635/`、`GSE4987/` 是下载原件；根目录的 `.py` 是生成代码；`outputs/`、`figures/` 是代码产物；`reports/` 和本目录的 `.md` 是交接材料。`outputs/` 可重建，两个 GEO 原件请保留。

## 快速开始

要求 Python 3.10+。核心数据处理、图表和测试只用 Python 标准库；不需要安装第三方包。在本目录执行：

```bash
cd 'ICA/数据准备'
python3 build_all.py
python3 -m unittest discover -s tests -v
```

`build_all.py` 会读取现有 `GSE3635/GSE3635_series_matrix` 与 `GSE4987/GSE4987_series_matrix`，重新生成 `outputs/` 和 `figures/`。它不会改动 GEO 下载原件。若保留了 `../A/data/GSE3635_RAW/` 中的 13 个 GenePix 文件，也会生成只读的 GPR 审查记录。生成 100 套模拟数据约需数秒至数十秒，取决于磁盘速度。

## 目录与文件契约

| 路径 | 用途及为什么需要 |
| --- | --- |
| `simulation.py` | 本组自写模拟生成算法；先确定类别，再加入噪声、相位、幅度、缺失和伸缩。CLI 可改变每类基因数、时间点数和 seed。 |
| `preprocessing.py` | 对模拟和真实数据执行相同的标签盲清洗；B/C 使用这里导出的无缺失矩阵。 |
| `prepare_real_data.py` | 解析 GEO series matrix 的元数据、按时间顺序拆分实验、筛掉控制或非标准 ORF，保留审计记录。 |
| `audit_gpr.py` | 审查既有 GSE3635 原始 spot 文件是否可读；spot 不是独立基因，**不作为聚类输入**。 |
| `generate_sweep.py` | 一次只改一个扰动强度，生成 100 个供 D 运行 benchmark 的输入条件。 |
| `make_figures.py` | 从实际输出绘制六类模板、五因素影响和 GEO 样本 QC 的 SVG；不包含虚构算法成绩。 |
| `build_all.py` | 一条命令重建全部 A 产物。 |
| `tests/test_data_preparation.py` | 验收配对设计、固定 seed、标签隔离、插值、GEO 结构和 100 份文件的校验和。 |
| `outputs/simulated/{smoke,baseline,challenge}/` | 小型接口检查、150 基因基础样例、较难联调样例；**它们不是已完成的 benchmark**。 |
| `outputs/simulated/sweep/` | 100 套正式受控输入及 `sweep_index.csv`；每套附配置、真标签、逐基因参数和预处理结果。 |
| `outputs/real/{GSE3635,GSE4987}/` | GEO 来源、样本、特征筛选、缺失处理、独立时间序列和校验和。 |
| `figures/` | 可复现的 SVG 图；图中模拟值为任意单位，GEO 值为已归一化 log ratio。 |
| `reports/data_preparation_report.md` | 本次运行的实测数据摘要、图表、验收和未完成事项；适合直接转交组员。 |

每套模拟数据的 `expression.csv` 仅包含 `gene_id,t0,t10,...` 及数值。**聚类程序只读取这个文件或 `processed/centered.csv`，不读取真标签**。`ground_truth.csv` 单独保存 `true_cluster`（六类，相位敏感）与 `shape_cluster`（早峰/晚峰合并为五类）。`gene_parameters.csv` 用于审计各基因实际抽到的宽度、基线、幅度、位移和伸缩；`metadata.json` 记录 seed、时间列、扰动级别和用途。缺失测量在 CSV 中为空白。

每个 `processed/` 有三张相同行顺序的矩阵：`raw_imputed.csv`（保留原数值尺度，但已插补，**不是原始下载文件**）、`centered.csv`（逐基因减均值，建议作为主比较输入）和 `zscore.csv`（再除以逐基因总体标准差，单列为形状敏感性分析）。另有 `row_log.csv` 和 `preprocessing_log.json` 记录剔除、插补和常值行。以 `gene_id` 连接清洗后的行与真标签，不能按原行号盲目比较。

## 已生成数据的实际规模

| 数据 | 来源特征 / 样本 | 清洗前缺失 | ORF 行 | 预处理后行 | 说明 |
| --- | ---: | ---: | ---: | ---: | --- |
| GSE3635 | 6,228 / 13 | 172 | 6,215 | 6,214 | `run01`，0–120 min，每 10 min；1 行超出 40% 缺失阈值。 |
| GSE4987 | 6,228 / 50 | 440 | 6,215 | 每组 6,215 | `run01`、`run02` 各 25 点，0–120 min，每 5 min；保留两组身份。 |
| smoke | 12 / 7 | 0 | — | 12 | 每类 2 条；只检查接口。 |
| baseline | 150 / 13 | 0 | — | 150 | 每类 25 条；基础样例。 |
| challenge | 150 / 13 | 629 | — | 122 | 5 种扰动同时较强；448 个单元格插补，28 行剔除。 |
| sweep | 100 套 × 150 / 13 | 随因素而变 | — | 见索引 | 五因素 × 四级 × 五 seed；每套独立日志。 |

两个 GEO 文件均为 GPL1914 酵母微阵列，GEO 描述数值为“normalized log ratio using Rosetta Resolver”，**不是 RNA-seq read counts**。GSE3635 官方说明无重复；GSE4987 说明包含 dye-swap 技术重复。两者的时间网格与实验批次不同，所以没有直接拼接或平均。GSE4987 两组也没有在此阶段平均。官方元数据可在 [GSE3635](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE3635) 与 [GSE4987](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE4987) 核对。

## 单独运行与修改参数

以下是各模块可独立重跑的命令；完整默认流程已通过 `build_all.py` 实跑。更改输出目录可保留已有结果。

```bash
python3 simulation.py --output outputs/simulated/smoke --genes-per-class 2 --timepoints 7 --seed 101 --purpose 'Small interface check'
python3 preprocessing.py --input outputs/simulated/smoke/expression.csv --output outputs/simulated/smoke/processed
python3 simulation.py --output outputs/simulated/baseline --genes-per-class 25 --timepoints 13 --seed 2026 --noise-sd 0.05 --max-shift 0.03 --amplitude-sd 0.15 --missing-rate 0 --stretch-sd 0.08 --purpose 'Baseline 150-gene example'
python3 generate_sweep.py --output outputs/simulated/sweep --genes-per-class 25 --timepoints 13 --first-seed 3100 --seeds 5
python3 prepare_real_data.py --accession GSE3635 --source GSE3635/GSE3635_series_matrix --output outputs/real/GSE3635
python3 prepare_real_data.py --accession GSE4987 --source GSE4987/GSE4987_series_matrix --output outputs/real/GSE4987
python3 make_figures.py
```

`--max-shift` 是总时长的比例；`--amplitude-sd` 和 `--stretch-sd` 是自然对数尺度的标准差；`--missing-rate` 是每个时间点独立遮蔽概率。若改动 sweep 的级别或生成算法，需重跑 `build_all.py` 与测试，并保留新版 `sweep_design.json`。

## 交接给 B/C/D/E

- **B（DTW）与 C（基线/聚类）**：先用 `outputs/simulated/smoke/processed/centered.csv` 调通；之后使用 `baseline` 或任一 sweep 的同名文件。所有方法用同一矩阵、相同 gene_id 顺序。不可把 NaN 直接送入 DTW；常值行的 Pearson 距离需共同约定规则，见设计文档。
- **D（benchmark）**：读取 `sweep_index.csv` 中的 `relative_dir`，每套用 `processed/centered.csv` 聚类，结束后按 `gene_id` 连接 `ground_truth.csv`。主任务用 `true_cluster`，次要形状任务用 `shape_cluster`，结果分开汇报。调窗口的 seed 应与最终评估 seed 分开。
- **E（真实数据）**：先用 `outputs/real/GSE3635/run01/processed/centered.csv`；GSE4987 两组各自分析。`real_data_cleaning_log.json` 汇总规则和保留数量，`feature_decisions.csv`、`sample_manifest.csv`、`data_preview_summary.json` 和逐行日志可追溯筛选。真实数据没有六类模拟真标签；生物学注释属于后续独立验证。

提交的 proposal 计划比较标准 DTW、三种单独的惩罚改法及 Euclidean/Pearson，并以 k-medoids 为主、平均连接层次聚类为补充。A 的输入对这些方法共用；B/D 应依据老师反馈明确最终测试的具体方法数、改法与文献版本的关系，并定义“clustering accuracy”的计算方式。A 的文件不预先绑定某一种距离或聚类实现。

## 限制与学术说明

模拟曲线是有意简化的解析模型，独立高斯噪声和随机缺失不能完全代表微阵列误差。强相位变化可能让早峰/晚峰变得难分，强缺失会改变保留基因数；这正是要在 benchmark 中报告的现象。所有图仅展示数据输入和质量，不代表算法性能。课程要求说明生成式 AI 的使用：本目录的代码和中文说明由 AI 协助起草、运行和检查；小组提交前应由成员亲自复核模型假设、代码、统计结果与最终表述，并在正式提交材料中如实披露。
