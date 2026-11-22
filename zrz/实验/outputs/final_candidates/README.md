# 最终候选 DTW 参数组合

## 怎么挑出来的

1. 网格搜索（`tune_dtw.ipynb`，2738 组），**不缩放**（`--scale-floor 0`），聚类用徐秋凡的 k-medoids，20 个种子取平均。
2. 在 wpj 的 baseline 数据上，挑出 `true_cluster`（6 类）ARI 平均 > 0.95 的组合，共 23 组（`final_selection.ipynb` 第 2 步）。
3. 同一网格在新的中等难度数据（`zrz/simulation_new`，种子 2028）上的成绩里查这 23 组，选出下面 10 组。

结果来源：`outputs/tuning_noscale/results.csv`（baseline）、`outputs/tuning_new_noscale/results.csv`（新数据）。

## 10 组候选（6 类标签 `true_cluster`）

| # | 配置 | 开启参数 | baseline ARI | baseline ARI 最差 | 新数据 ARI | 新数据 ARI 最差 | 新数据 NMI |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | window=0.1667, extend_penalty=0.05, time_penalty=2 | 3 | 0.9509 | 0.7414 | 0.7961 | 0.6345 | 0.8396 |
| 2 | window=0.25, time_penalty=2 | 2 | 0.9623 | 0.7414 | 0.7909 | 0.6428 | 0.8408 |
| 3 | time_penalty=2 | 1 | 0.9618 | 0.7320 | 0.7908 | 0.6428 | 0.8398 |
| 4 | window=0.1667, time_penalty=2 | 2 | 0.9623 | 0.7349 | 0.7821 | 0.6330 | 0.8345 |
| 5 | window=0.1667, extend_penalty=0.05, time_penalty=0.5, dynamics=0.3 | 4 | 0.9521 | 0.7320 | 0.7728 | 0.6338 | 0.8283 |
| 6 | window=0.25, extend_penalty=0.05, time_penalty=0.5, dynamics=0.3 | 4 | 0.9532 | 0.7471 | 0.7656 | 0.6304 | 0.8247 |
| 7 | window=0.1667, extend_penalty=0.5, time_penalty=2 | 3 | 0.9622 | 0.7364 | 0.7654 | 0.5995 | 0.8154 |
| 8 | extend_penalty=0.05, time_penalty=0.5, dynamics=0.3 | 3 | 0.9532 | 0.7471 | 0.7609 | 0.6363 | 0.8240 |
| 9 | window=0.25, step_penalty=0.1, extend_penalty=0.5 | 3 | 0.9537 | 0.7541 | 0.7606 | 0.6025 | 0.8165 |
| 10 | window=0.1667, extend_penalty=0.5, time_penalty=0.5, dynamics=0.3 | 4 | 0.9535 | 0.7541 | 0.7580 | 0.5828 | 0.8123 |

参照（ARI 平均）：欧氏 baseline 0.9175 / 新数据 0.5341；皮尔森 baseline 0.7324 / 新数据 0.4973；标准 DTW baseline 0.6420 / 新数据 0.5153。

## dtw.py 命令

```bash
# 1
python3 dtw.py --input <centered.csv> --output outputs/<name> --window 0.1667 --extend-penalty 0.05 --time-penalty 2 --scale-floor 0
# 2
python3 dtw.py --input <centered.csv> --output outputs/<name> --window 0.25 --time-penalty 2 --scale-floor 0
# 3
python3 dtw.py --input <centered.csv> --output outputs/<name> --time-penalty 2 --scale-floor 0
# 4
python3 dtw.py --input <centered.csv> --output outputs/<name> --window 0.1667 --time-penalty 2 --scale-floor 0
# 5
python3 dtw.py --input <centered.csv> --output outputs/<name> --window 0.1667 --extend-penalty 0.05 --time-penalty 0.5 --dynamics 0.3 --scale-floor 0
# 6
python3 dtw.py --input <centered.csv> --output outputs/<name> --window 0.25 --extend-penalty 0.05 --time-penalty 0.5 --dynamics 0.3 --scale-floor 0
# 7
python3 dtw.py --input <centered.csv> --output outputs/<name> --window 0.1667 --extend-penalty 0.5 --time-penalty 2 --scale-floor 0
# 8
python3 dtw.py --input <centered.csv> --output outputs/<name> --extend-penalty 0.05 --time-penalty 0.5 --dynamics 0.3 --scale-floor 0
# 9
python3 dtw.py --input <centered.csv> --output outputs/<name> --window 0.25 --step-penalty 0.1 --extend-penalty 0.5 --scale-floor 0
# 10
python3 dtw.py --input <centered.csv> --output outputs/<name> --window 0.1667 --extend-penalty 0.5 --time-penalty 0.5 --dynamics 0.3 --scale-floor 0
```

全部指标（两份数据 × 两种标签）见同目录 `candidates.csv`。
