# DTW 网格搜索结果

数据：`../../wpj/predata/数据准备/outputs/simulated/baseline`；scale_floor=0；聚类：徐秋凡的 k-medoids，20 个种子取平均。

网格：`{"window": [0, 0.0834, 0.1667, 0.25], "step_penalty": [0, 0.1, 0.5, 2], "open_penalty": [0, 0.2, 2], "extend_penalty": [0, 0.05, 0.5], "time_penalty": [0, 0.5, 2, 10], "power": [1, 2], "dynamics": [0, 0.3, 1]}`

排序：ARI 平均 → ARI 最差种子 → NMI → 真标签轮廓系数 → 开启的参数越少越好 → 越快越好。前几名之间差距很小时，可能只是聚类的随机波动。

### 标签 true_cluster（k=6），共 2738 组

| 排名 | 配置 | ARI 平均 | ARI 最差 | NMI | 真标签轮廓 | 早晚峰同簇率 | 秒 |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | window=0.25, extend_penalty=0.5, time_penalty=0.5, dynamics=0.3 | 0.9651 | 0.7284 | 0.9836 | 0.6301 | 0.00 | 0.0 |
| 2 | window=0.1667, time_penalty=2 | 0.9623 | 0.7349 | 0.9799 | 0.6396 | 0.00 | 0.0 |
| 3 | window=0.25, time_penalty=2 | 0.9623 | 0.7414 | 0.9808 | 0.6257 | 0.00 | 0.0 |
| 4 | window=0.1667, extend_penalty=0.5, time_penalty=2 | 0.9622 | 0.7364 | 0.9804 | 0.6418 | 0.00 | 0.0 |
| 5 | window=0.1667, extend_penalty=0.05, time_penalty=2, dynamics=0.3 | 0.9618 | 0.7233 | 0.9805 | 0.6342 | 0.00 | 0.0 |
| 6 | time_penalty=2 | 0.9618 | 0.7320 | 0.9804 | 0.6246 | 0.00 | 0.0 |
| 7 | window=0.1667, step_penalty=0.1, extend_penalty=0.5, time_penalty=0.5 | 0.9615 | 0.7284 | 0.9796 | 0.6250 | 0.00 | 0.0 |
| 8 | window=0.25, step_penalty=0.1, open_penalty=0.2 | 0.9613 | 0.7313 | 0.9792 | 0.5899 | 0.00 | 0.0 |
| 9 | window=0.1667, extend_penalty=0.5, time_penalty=2, dynamics=0.3 | 0.9613 | 0.7123 | 0.9801 | 0.6345 | 0.00 | 0.0 |
| 10 | window=0.1667, step_penalty=0.1, time_penalty=0.5, dynamics=0.3 | 0.9588 | 0.7219 | 0.9774 | 0.6228 | 0.00 | 0.0 |
| 11 | window=0.1667, time_penalty=0.5, dynamics=0.3 | 0.9539 | 0.7681 | 0.9792 | 0.6401 | 0.00 | 0.0 |
| 12 | window=0.25, step_penalty=0.1, extend_penalty=0.5 | 0.9537 | 0.7541 | 0.9793 | 0.6153 | 0.00 | 0.0 |
| 13 | extend_penalty=0.5, time_penalty=0.5, dynamics=0.3 | 0.9536 | 0.7284 | 0.9787 | 0.6294 | 0.00 | 0.0 |
| 14 | window=0.1667, extend_penalty=0.5, time_penalty=0.5, dynamics=0.3 | 0.9535 | 0.7541 | 0.9771 | 0.6446 | 0.00 | 0.0 |
| 15 | window=0.25, extend_penalty=0.05, time_penalty=0.5, dynamics=0.3 | 0.9532 | 0.7471 | 0.9787 | 0.6117 | 0.00 | 0.0 |
| 16 | extend_penalty=0.05, time_penalty=0.5, dynamics=0.3 | 0.9532 | 0.7471 | 0.9787 | 0.6069 | 0.00 | 0.0 |
| 17 | window=0.1667, extend_penalty=0.05, time_penalty=0.5, dynamics=0.3 | 0.9521 | 0.7320 | 0.9759 | 0.6433 | 0.00 | 0.0 |
| 18 | window=0.25, step_penalty=0.1, extend_penalty=0.5, time_penalty=0.5 | 0.9514 | 0.7364 | 0.9761 | 0.6198 | 0.00 | 0.0 |
| 19 | step_penalty=0.1, extend_penalty=0.5, time_penalty=0.5 | 0.9514 | 0.7364 | 0.9761 | 0.6196 | 0.00 | 0.0 |
| 20 | window=0.1667, extend_penalty=0.05, time_penalty=2 | 0.9509 | 0.7414 | 0.9751 | 0.6410 | 0.00 | 0.0 |
| 1587 | euclidean | 0.9175 | 0.7430 | 0.9633 | 0.5900 | 0.00 | 0.0 |
| 2737 | pearson | 0.7324 | 0.7053 | 0.7879 | 0.6556 | 0.00 | 0.1 |
| 2738 | standard DTW | 0.6420 | 0.4699 | 0.7754 | 0.4198 | 0.32 | 0.1 |

推荐：**window=0.25, extend_penalty=0.5, time_penalty=0.5, dynamics=0.3**  
`python3 dtw.py --input <centered.csv> --output outputs/<name> --window 0.25 --extend-penalty 0.5 --time-penalty 0.5 --dynamics 0.3 --scale-floor 0`

每个参数取某个值时，所有含这个值的组合的平均 ARI 和最好 ARI（看哪些值普遍有帮助）：

| 参数 | 值 | 组合数 | 平均 ARI | 最好 ARI |
| --- | ---: | ---: | ---: | ---: |
| window | 0 | 684 | 0.9199 | 0.9618 |
| window | 0.0834 | 684 | 0.9204 | 0.9432 |
| window | 0.1667 | 684 | 0.9202 | 0.9623 |
| window | 0.25 | 684 | 0.9203 | 0.9651 |
| step_penalty | 0 | 684 | 0.9227 | 0.9651 |
| step_penalty | 0.1 | 684 | 0.9215 | 0.9615 |
| step_penalty | 0.5 | 684 | 0.9190 | 0.9422 |
| step_penalty | 2 | 684 | 0.9175 | 0.9175 |
| open_penalty | 0 | 912 | 0.9222 | 0.9651 |
| open_penalty | 0.2 | 912 | 0.9207 | 0.9613 |
| open_penalty | 2 | 912 | 0.9176 | 0.9290 |
| extend_penalty | 0 | 912 | 0.9201 | 0.9623 |
| extend_penalty | 0.05 | 912 | 0.9201 | 0.9618 |
| extend_penalty | 0.5 | 912 | 0.9204 | 0.9651 |
| time_penalty | 0 | 144 | 0.9180 | 0.9613 |
| time_penalty | 0.5 | 864 | 0.9222 | 0.9651 |
| time_penalty | 2 | 864 | 0.9202 | 0.9623 |
| time_penalty | 10 | 864 | 0.9184 | 0.9317 |
| power | 1 | 1296 | 0.9185 | 0.9477 |
| power | 2 | 1296 | 0.9221 | 0.9651 |
| dynamics | 0 | 864 | 0.9219 | 0.9623 |
| dynamics | 0.3 | 864 | 0.9201 | 0.9651 |
| dynamics | 1 | 864 | 0.9189 | 0.9496 |

### 标签 shape_cluster（k=5），共 2738 组

| 排名 | 配置 | ARI 平均 | ARI 最差 | NMI | 真标签轮廓 | 早晚峰同簇率 | 秒 |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | open_penalty=0.2 | 0.7072 | 0.6346 | 0.8411 | 0.4492 | 0.31 | 0.0 |
| 2 | window=0.0834, extend_penalty=0.5, time_penalty=0.5, power=1, dynamics=0.3 | 0.6694 | 0.6275 | 0.8492 | 0.4076 | 0.05 | 0.0 |
| 3 | window=0.0834, extend_penalty=0.05, time_penalty=0.5, power=1, dynamics=0.3 | 0.6694 | 0.6275 | 0.8492 | 0.4076 | 0.05 | 0.0 |
| 4 | window=0.0834, time_penalty=0.5, power=1, dynamics=0.3 | 0.6694 | 0.6275 | 0.8492 | 0.4076 | 0.05 | 0.0 |
| 5 | window=0.0834, open_penalty=0.2, extend_penalty=0.05, time_penalty=0.5, power=1 | 0.6694 | 0.6275 | 0.8492 | 0.4041 | 0.05 | 0.0 |
| 6 | window=0.0834, open_penalty=0.2, extend_penalty=0.5, time_penalty=0.5, power=1 | 0.6694 | 0.6275 | 0.8492 | 0.4041 | 0.05 | 0.0 |
| 7 | window=0.0834, open_penalty=0.2, time_penalty=0.5, power=1 | 0.6694 | 0.6275 | 0.8492 | 0.4041 | 0.05 | 0.0 |
| 8 | window=0.0834, time_penalty=2, power=1 | 0.6694 | 0.6275 | 0.8492 | 0.3982 | 0.05 | 0.0 |
| 9 | window=0.0834, extend_penalty=0.05, time_penalty=2, power=1 | 0.6694 | 0.6275 | 0.8492 | 0.3982 | 0.05 | 0.0 |
| 10 | window=0.0834, extend_penalty=0.5, time_penalty=2, power=1 | 0.6694 | 0.6275 | 0.8492 | 0.3982 | 0.05 | 0.0 |
| 11 | extend_penalty=0.5, time_penalty=2, power=1 | 0.6694 | 0.6275 | 0.8492 | 0.3981 | 0.05 | 0.0 |
| 12 | window=0.1667, extend_penalty=0.5, time_penalty=2, power=1 | 0.6694 | 0.6275 | 0.8492 | 0.3981 | 0.05 | 0.0 |
| 13 | window=0.25, extend_penalty=0.5, time_penalty=2, power=1 | 0.6694 | 0.6275 | 0.8492 | 0.3981 | 0.05 | 0.0 |
| 14 | extend_penalty=0.05, time_penalty=2, power=1 | 0.6694 | 0.6275 | 0.8492 | 0.3980 | 0.05 | 0.0 |
| 15 | window=0.1667, extend_penalty=0.05, time_penalty=2, power=1 | 0.6694 | 0.6275 | 0.8492 | 0.3980 | 0.05 | 0.0 |
| 16 | window=0.25, extend_penalty=0.05, time_penalty=2, power=1 | 0.6694 | 0.6275 | 0.8492 | 0.3980 | 0.05 | 0.0 |
| 17 | window=0.0834, extend_penalty=0.05, time_penalty=0.5, dynamics=1 | 0.6687 | 0.6275 | 0.8491 | 0.4357 | 0.00 | 0.0 |
| 18 | window=0.0834, extend_penalty=0.5, time_penalty=0.5, dynamics=1 | 0.6687 | 0.6275 | 0.8491 | 0.4357 | 0.00 | 0.0 |
| 19 | window=0.0834, extend_penalty=0.5, time_penalty=2, dynamics=0.3 | 0.6687 | 0.6275 | 0.8491 | 0.4282 | 0.00 | 0.0 |
| 20 | window=0.0834, time_penalty=2, dynamics=0.3 | 0.6687 | 0.6275 | 0.8491 | 0.4282 | 0.00 | 0.0 |
| 1218 | euclidean | 0.6579 | 0.3987 | 0.8441 | 0.3905 | 0.05 | 0.0 |
| 2613 | standard DTW | 0.6487 | 0.4700 | 0.8117 | 0.3827 | 0.63 | 0.1 |
| 2738 | pearson | 0.5825 | 0.4761 | 0.7025 | 0.3758 | 0.00 | 0.1 |

推荐：**open_penalty=0.2**  
`python3 dtw.py --input <centered.csv> --output outputs/<name> --open-penalty 0.2 --scale-floor 0`

每个参数取某个值时，所有含这个值的组合的平均 ARI 和最好 ARI（看哪些值普遍有帮助）：

| 参数 | 值 | 组合数 | 平均 ARI | 最好 ARI |
| --- | ---: | ---: | ---: | ---: |
| window | 0 | 684 | 0.6563 | 0.7072 |
| window | 0.0834 | 684 | 0.6582 | 0.6694 |
| window | 0.1667 | 684 | 0.6567 | 0.6694 |
| window | 0.25 | 684 | 0.6562 | 0.6694 |
| step_penalty | 0 | 684 | 0.6560 | 0.7072 |
| step_penalty | 0.1 | 684 | 0.6558 | 0.6686 |
| step_penalty | 0.5 | 684 | 0.6577 | 0.6580 |
| step_penalty | 2 | 684 | 0.6579 | 0.6579 |
| open_penalty | 0 | 912 | 0.6559 | 0.6694 |
| open_penalty | 0.2 | 912 | 0.6567 | 0.7072 |
| open_penalty | 2 | 912 | 0.6579 | 0.6579 |
| extend_penalty | 0 | 912 | 0.6567 | 0.7072 |
| extend_penalty | 0.05 | 912 | 0.6567 | 0.6694 |
| extend_penalty | 0.5 | 912 | 0.6572 | 0.6694 |
| time_penalty | 0 | 144 | 0.6517 | 0.7072 |
| time_penalty | 0.5 | 864 | 0.6555 | 0.6694 |
| time_penalty | 2 | 864 | 0.6577 | 0.6694 |
| time_penalty | 10 | 864 | 0.6583 | 0.6672 |
| power | 1 | 1296 | 0.6578 | 0.6694 |
| power | 2 | 1296 | 0.6564 | 0.6687 |
| dynamics | 0 | 864 | 0.6568 | 0.6694 |
| dynamics | 0.3 | 864 | 0.6569 | 0.6694 |
| dynamics | 1 | 864 | 0.6577 | 0.6687 |
