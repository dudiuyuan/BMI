# Project 11 · Member E 可运行交付

本目录实现 GSE3635 来源核查、数据清洗、无监督聚类、生物学富集验证和通用 CSV 命令行。完整结果及方法解释见 [Member_E_summary.md](Member_E_summary.md)。

## 快速运行（Windows PowerShell）

建议使用 Python 3.13（本次实测版本），在此目录执行。数据已经随目录保存，重跑无需联网下载数据；首次安装依赖需要联网。

```powershell
cd D:\BMI\project11
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements-lock.txt
.\.venv\Scripts\python.exe -m unittest discover -s tests -v
.\.venv\Scripts\python.exe main.py --input data/gse3635/gene_log10ratio_with_missing.csv --method dtw --k 6 --window 1 --output results/example_cli --seed 11
```

这里的 `k=6, window=1` 是演示命令；正式分析使用 `results/biological_validation/selected_parameters.json` 中由预先指定规则选出的参数。所有 CLI 参数都可修改，无需编辑源码。基线切换为 `--method euclidean` 或 `--method pearson`；保留幅度信息用 `--scaling none`。`--window` 是采样点索引半径，在本数据上每单位对应 10 分钟，基线忽略该参数。

## 从原始文件完整重建

```powershell
# 如来源文件已在 data 中，跳过此下载步骤。
powershell -ExecutionPolicy Bypass -File download_sources.ps1
.\.venv\Scripts\python.exe prepare_gse3635.py --archive ../GSE3635_RAW.tar
.\.venv\Scripts\python.exe analyze_real.py
.\.venv\Scripts\python.exe simulation.py --n 120 --t 25 --seed 11 --output data/simulated
.\.venv\Scripts\python.exe main.py --input data/simulated/expression.csv --method dtw --k 6 --window 2 --output results/simulation_cli --seed 11
.\.venv\Scripts\python.exe benchmark.py
.\.venv\Scripts\python.exe build_report.py
```

已有原始 tar 也保存在 `data/gse3635/raw/GSE3635_RAW.tar`，若原目录未保留，可用这个路径替换 `--archive`。正式分析是全部保留基因分析；距离矩阵为 N×N、float32，单矩阵约 144 MiB，运行时另有中间数组，建议留出数 GB 内存。首次 Numba 编译会增加启动耗时。

## 输入输出协议

CSV 第一列 `gene_id`，其余列为数值时间名 `t0,t10,t20,...`；行是基因，列是时间点。时间列按数值排序，允许不等间隔；至少两个不同时间点。不要把 `true_cluster` 或其他注释列混入表达矩阵。缺失值留空/NaN；非数字文本报错；无限值按缺失处理。gene ID 需先与参考统一；通用 CLI 保留传入 ID，不擅自推断别名。

`main.py` 输出 `assignments.csv`、`medoids.csv`、`medoid_curves.png`、`heatmap.png`、`processed_input.csv`、`preprocessing_audit.csv`、`run_config.json`。重复 ID 逐时间点取中位数；默认剔除缺失比例大于 20% 或有效点少于 2 的基因；按实际时间线性插值，端点采用最近值；z-score 使用总体标准差（ddof=0），常值序列变为零。z-score 去掉幅度及基因均值，不应据此声称幅度恢复准确。

单独验证其他成员产生的分簇：

```powershell
.\.venv\Scripts\python.exe biological_validation.py --assignments results/biological_validation/dtw/assignments.csv --reference data/references/spellman800.csv --expression results/biological_validation/dtw/processed_input.csv --output results/revalidate
```

## 代码与组员对接

| 文件 | 接口 / 对接对象 |
|---|---|
| `preprocessing.py` | `preprocess(frame, scaling, max_missing) -> clean_frame, log, audit`；与 A 共享 |
| `dtw.py` | `dtw(x, y, window)`；与 B 对接；核心 DP 自实现，Numba 仅加速 |
| `distance_matrix.py` | `distance_matrix(X, method, window) -> D`；统一三种距离 |
| `clustering.py` | `cluster(D, k, seed, n_init) -> labels, medoid_indices, objective, iterations`；与 C 对接 |
| `main.py` | 通用入口；不导入生物学参考或模拟真值 |
| `analyze_real.py` | 18 组 k/window 无监督搜索、三种距离完整运行、最后才读参考 |
| `biological_validation.py` | Fisher、BH、phase 分布、参考匹配审计；E 主体 |
| `simulation.py` / `benchmark.py` | 五类扰动的集成补充验证；待与 A/D 正式版本合并 |
| `build_report.py` | 由结果表自动生成中文总结，不手填统计成绩 |

距离为 DTW 累加平方局部代价后开方，不按路径长度归一化；窗口 `|i-j|<=w`。不同长度的单对序列要求 `w>=abs(n-m)`。无窗口可用单对 `window=None`；CSV 矩阵中所有序列共享列。单对复杂度 O(Tw)，rolling rows 空间 O(T)，完整矩阵空间 O(N²)。窗口 0 与 Euclidean 一致。Pearson 距离 `1-r`，常值与非常值约定距离 1、两个常值距离 0。k-medoids 是多起点交替更新版本，不保证全局最优，不等同于完整 PAM swap。

当前工作区未提供 A–D 成员源码、真实版本号或成员姓名，因此本目录中的相关模块是供集成使用的独立实现。参见 `docs/integration_handover.md`，其待填字段须由组员提供，不能冒充已收集完成。

## 数据使用和结论边界

原始 GEO 与平台信息保留其来源；平台描述为非商业分发。yeastCC 包许可为 Artistic-2.0；其中历史 `CellCycle98.txt` 第一行另明确标注非商业用途，已原样保留。本交付用于课程研究；对外重新分发应保留原始声明并遵守其适用条款。

真实数据富集结果仅代表参考基因集合的生物学支持，不是聚类准确率。800 基因不是所有真实簇的完整 Ground Truth。模拟 ARI/NMI 单独报告；这套小规模集成 benchmark 不能替代 D 的正式评测、额外生成算法设计及规模实验。
