# DTW

调参过程和全部实验在同级的 `实验/` 文件夹里。

只需要 Python 3.10 或更新版本，不依赖任何第三方库。

| 文件 | 作用 |
| --- | --- |
| `dtw.py` | DTW，五种固定模式 |
| `euclidean.py` | 对照：欧氏距离 |
| `pearson.py` | 对照：皮尔森相关距离（1 − r） |

三个程序的输入、输出格式完全相同，都直接使用输入表里的数值（不做逐基因缩放）。

## 用法

```bash
python3 dtw.py --mode diagonal-time --input <centered.csv> --output outputs/<名字>
python3 dtw.py --help                                   # 查看五种模式

python3 euclidean.py --input <centered.csv> --output outputs/<名字>
python3 pearson.py   --input <centered.csv> --output outputs/<名字>
```

每次实验请用新的输出文件夹，否则会覆盖。

## 五种模式

| 模式 | 含义 | 参数 |
| --- | --- | --- |
| `standard` | 标准 DTW，作对照 | 全部关闭 |
| `time` | 时间惩罚 | 时间惩罚 2 |
| `time-dynamics` | 时间惩罚 + 表达变化速率 | 时间惩罚 0.5、速率 0.3、延伸惩罚 0.05 |
| `diagonal-time-dynamics` | 对角线窗口 + 时间惩罚 + 表达变化速率 | 窗口 0.25、时间惩罚 0.5、速率 0.3、延伸惩罚 0.05 |
| `diagonal-time` | 对角线窗口 + 时间惩罚 | 窗口 0.1667、时间惩罚 2、延伸惩罚 0.05 |

- **时间惩罚**：两个时间点相距越远，配在一起的代价越高，防止把相隔很远的峰硬对齐。
- **表达变化速率**：表达变化快的地方，加重时间惩罚。
- **对角线窗口**：最多只允许错开总时长的一定比例（0.25 = 25%，0.1667 ≈ 2 个采样点），也更快。
- **延伸惩罚**：同一段拉伸从第二步起每步额外收费，防止连续拉伸很长一段。

## 输入

CSV，第一列 `gene_id`，其余列为按时间排序的表达值，列名 `t0,t10,t20,...`（数字为分钟）；每行一个基因，不能有缺失值。直接用 wpj 数据准备输出的 `processed/centered.csv`（已补缺失、逐基因减均值）。不要把真标签放进输入表。

## 输出

- `distance_matrix.csv`：所有基因两两之间的距离（对称矩阵，第一行、第一列是 gene_id，对角线为 0），交给下游聚类。
- `run.json`：方法、模式和参数、输入文件路径、基因数、时间点。

## 下游聚类

`distance_matrix.csv` 交给xqf的 k-medoids（`xqf/xqf_clustering/k_medoids.py`）。

## 运行时间

150 个基因 × 13 个时间点实测：`diagonal-time` 1.3 秒，`diagonal-time-dynamics` 1.6 秒，`time`、`time-dynamics`、`standard` 约 3 秒，欧氏、皮尔森 0.3 秒。耗时随基因数的平方增长：GSE3635 全部 6000 多个基因约 1800 万对，带窗口的模式约半小时，不带窗口的一个多小时。

## 参数是怎么选的

1. 在 wpj 的 baseline 模拟数据（150 个基因）上做网格搜索：7 个参数的全部组合，去掉等价的后 2736 组；每组用徐秋凡的 k-medoids 跑 20 个种子，按 6 类真标签（`true_cluster`，早峰、晚峰算两类）算 ARI 平均。
2. 挑出 baseline 上 ARI > 0.95 的 23 组，看它们在一份没参与调参的中等难度新数据（`zrz/simulation_new`）上的成绩，选出四组。

ARI 平均（6 类标签，20 个种子）：

| | baseline | 新数据 |
| --- | ---: | ---: |
| `time` | 0.962 | 0.791 |
| `time-dynamics` | 0.953 | 0.761 |
| `diagonal-time-dynamics` | 0.953 | 0.766 |
| `diagonal-time` | 0.951 | 0.796 |
| `standard` | 0.642 | 0.515 |
| 欧氏距离 | 0.918 | 0.534 |
| 皮尔森 | 0.732 | 0.497 |

四组在新数据上都排在全部 2736 组的前 1%。参数是按 6 类标签选的。详细过程见 `实验/outputs/final_candidates/`。
