# DTW 网格搜索结果

数据：`../simulation_new`；scale_floor=0；聚类：徐秋凡的 k-medoids，20 个种子取平均。

网格：`{"window": [0, 0.0834, 0.1667, 0.25], "step_penalty": [0, 0.1, 0.5, 2], "open_penalty": [0, 0.2, 2], "extend_penalty": [0, 0.05, 0.5], "time_penalty": [0, 0.5, 2, 10], "power": [1, 2], "dynamics": [0, 0.3, 1]}`

排序：ARI 平均 → ARI 最差种子 → NMI → 真标签轮廓系数 → 开启的参数越少越好 → 越快越好。前几名之间差距很小时，可能只是聚类的随机波动。

### 标签 true_cluster（k=6），共 2738 组

| 排名 | 配置 | ARI 平均 | ARI 最差 | NMI | 真标签轮廓 | 早晚峰同簇率 | 秒 |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | window=0.1667, extend_penalty=0.05, time_penalty=2 | 0.7961 | 0.6345 | 0.8396 | 0.2660 | 0.00 | 0.0 |
| 2 | window=0.25, time_penalty=2 | 0.7909 | 0.6428 | 0.8408 | 0.2591 | 0.00 | 0.0 |
| 3 | time_penalty=2 | 0.7908 | 0.6428 | 0.8398 | 0.2563 | 0.00 | 0.0 |
| 4 | window=0.1667, time_penalty=2 | 0.7821 | 0.6330 | 0.8345 | 0.2643 | 0.00 | 0.0 |
| 5 | window=0.25, extend_penalty=0.05, time_penalty=2 | 0.7803 | 0.6425 | 0.8352 | 0.2618 | 0.00 | 0.0 |
| 6 | extend_penalty=0.05, time_penalty=2 | 0.7803 | 0.6425 | 0.8352 | 0.2598 | 0.00 | 0.0 |
| 7 | window=0.1667, extend_penalty=0.05, time_penalty=0.5, dynamics=1 | 0.7796 | 0.6320 | 0.8271 | 0.2616 | 0.00 | 0.0 |
| 8 | window=0.1667, extend_penalty=0.05, time_penalty=0.5, dynamics=0.3 | 0.7728 | 0.6338 | 0.8283 | 0.2753 | 0.00 | 0.0 |
| 9 | window=0.1667, extend_penalty=0.5, time_penalty=0.5, dynamics=1 | 0.7684 | 0.6505 | 0.8215 | 0.2616 | 0.00 | 0.0 |
| 10 | window=0.25, extend_penalty=0.5, time_penalty=2 | 0.7667 | 0.5995 | 0.8171 | 0.2637 | 0.00 | 0.0 |
| 11 | extend_penalty=0.5, time_penalty=2 | 0.7667 | 0.5995 | 0.8171 | 0.2628 | 0.00 | 0.0 |
| 12 | window=0.25, extend_penalty=0.05, time_penalty=0.5, dynamics=1 | 0.7660 | 0.6180 | 0.8191 | 0.2583 | 0.00 | 0.0 |
| 13 | extend_penalty=0.05, time_penalty=0.5, dynamics=1 | 0.7660 | 0.6180 | 0.8191 | 0.2571 | 0.00 | 0.0 |
| 14 | window=0.25, extend_penalty=0.05, time_penalty=0.5, dynamics=0.3 | 0.7656 | 0.6304 | 0.8247 | 0.2647 | 0.00 | 0.0 |
| 15 | window=0.1667, extend_penalty=0.5, time_penalty=2 | 0.7654 | 0.5995 | 0.8154 | 0.2657 | 0.00 | 0.0 |
| 16 | extend_penalty=0.05, time_penalty=0.5, dynamics=0.3 | 0.7609 | 0.6363 | 0.8240 | 0.2591 | 0.01 | 0.0 |
| 17 | window=0.1667, step_penalty=0.1 | 0.7607 | 0.6259 | 0.8209 | 0.2520 | 0.00 | 0.0 |
| 18 | window=0.25, step_penalty=0.1, extend_penalty=0.5 | 0.7606 | 0.6025 | 0.8165 | 0.2472 | 0.00 | 0.0 |
| 19 | step_penalty=0.1, extend_penalty=0.5 | 0.7606 | 0.6025 | 0.8165 | 0.2447 | 0.00 | 0.0 |
| 20 | window=0.25, step_penalty=0.1 | 0.7601 | 0.6146 | 0.8186 | 0.2387 | 0.00 | 0.0 |
| 1377 | euclidean | 0.5341 | 0.3999 | 0.6668 | 0.1293 | 0.04 | 0.0 |
| 2737 | standard DTW | 0.5153 | 0.4097 | 0.6416 | 0.1972 | 0.26 | 0.1 |
| 2738 | pearson | 0.4973 | 0.4649 | 0.5664 | 0.2865 | 0.00 | 0.1 |

推荐：**window=0.1667, extend_penalty=0.05, time_penalty=2**  
`python3 dtw.py --input <centered.csv> --output outputs/<name> --window 0.1667 --extend-penalty 0.05 --time-penalty 2 --scale-floor 0`

每个参数取某个值时，所有含这个值的组合的平均 ARI 和最好 ARI（看哪些值普遍有帮助）：

| 参数 | 值 | 组合数 | 平均 ARI | 最好 ARI |
| --- | ---: | ---: | ---: | ---: |
| window | 0 | 684 | 0.5574 | 0.7908 |
| window | 0.0834 | 684 | 0.5575 | 0.7561 |
| window | 0.1667 | 684 | 0.5589 | 0.7961 |
| window | 0.25 | 684 | 0.5578 | 0.7909 |
| step_penalty | 0 | 684 | 0.5847 | 0.7961 |
| step_penalty | 0.1 | 684 | 0.5708 | 0.7607 |
| step_penalty | 0.5 | 684 | 0.5450 | 0.6058 |
| step_penalty | 2 | 684 | 0.5311 | 0.5353 |
| open_penalty | 0 | 912 | 0.5841 | 0.7961 |
| open_penalty | 0.2 | 912 | 0.5584 | 0.7263 |
| open_penalty | 2 | 912 | 0.5311 | 0.5353 |
| extend_penalty | 0 | 912 | 0.5580 | 0.7909 |
| extend_penalty | 0.05 | 912 | 0.5581 | 0.7961 |
| extend_penalty | 0.5 | 912 | 0.5576 | 0.7684 |
| time_penalty | 0 | 144 | 0.5891 | 0.7607 |
| time_penalty | 0.5 | 864 | 0.5721 | 0.7796 |
| time_penalty | 2 | 864 | 0.5569 | 0.7961 |
| time_penalty | 10 | 864 | 0.5395 | 0.6742 |
| power | 1 | 1296 | 0.5410 | 0.7400 |
| power | 2 | 1296 | 0.5713 | 0.7961 |
| dynamics | 0 | 864 | 0.5634 | 0.7961 |
| dynamics | 0.3 | 864 | 0.5563 | 0.7728 |
| dynamics | 1 | 864 | 0.5488 | 0.7796 |

### 标签 shape_cluster（k=5），共 2738 组

| 排名 | 配置 | ARI 平均 | ARI 最差 | NMI | 真标签轮廓 | 早晚峰同簇率 | 秒 |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | window=0.25, time_penalty=2 | 0.5476 | 0.4793 | 0.7010 | 0.1485 | 0.01 | 0.0 |
| 2 | time_penalty=2 | 0.5453 | 0.4793 | 0.6982 | 0.1459 | 0.01 | 0.0 |
| 3 | extend_penalty=0.05, time_penalty=0.5 | 0.5447 | 0.4507 | 0.6994 | 0.1445 | 0.03 | 0.0 |
| 4 | extend_penalty=0.05, time_penalty=2 | 0.5439 | 0.4661 | 0.6972 | 0.1477 | 0.00 | 0.0 |
| 5 | window=0.25, extend_penalty=0.05, time_penalty=2 | 0.5432 | 0.4661 | 0.6971 | 0.1497 | 0.01 | 0.0 |
| 6 | extend_penalty=0.05, time_penalty=0.5, dynamics=0.3 | 0.5422 | 0.4353 | 0.6969 | 0.1462 | 0.02 | 0.0 |
| 7 | extend_penalty=0.5, time_penalty=0.5 | 0.5412 | 0.4919 | 0.6938 | 0.1532 | 0.02 | 0.0 |
| 8 | window=0.1667, extend_penalty=0.05, time_penalty=0.5 | 0.5411 | 0.3728 | 0.6930 | 0.1633 | 0.01 | 0.0 |
| 9 | window=0.25, extend_penalty=0.5, time_penalty=0.5 | 0.5406 | 0.4919 | 0.6938 | 0.1570 | 0.02 | 0.0 |
| 10 | window=0.25, extend_penalty=0.05, time_penalty=0.5, dynamics=0.3 | 0.5385 | 0.4353 | 0.6943 | 0.1513 | 0.03 | 0.0 |
| 11 | time_penalty=0.5 | 0.5380 | 0.4104 | 0.6966 | 0.1382 | 0.04 | 0.0 |
| 12 | window=0.1667, time_penalty=0.5 | 0.5371 | 0.4315 | 0.6929 | 0.1602 | 0.05 | 0.0 |
| 13 | window=0.1667, extend_penalty=0.05, time_penalty=0.5, dynamics=0.3 | 0.5366 | 0.3729 | 0.6914 | 0.1586 | 0.01 | 0.0 |
| 14 | window=0.1667, extend_penalty=0.05, time_penalty=2 | 0.5361 | 0.3363 | 0.6925 | 0.1535 | 0.00 | 0.0 |
| 15 | window=0.1667, time_penalty=2 | 0.5352 | 0.3363 | 0.6939 | 0.1528 | 0.00 | 0.0 |
| 16 | window=0.25, time_penalty=0.5 | 0.5350 | 0.4504 | 0.6898 | 0.1456 | 0.02 | 0.0 |
| 17 | window=0.25, time_penalty=0.5, dynamics=0.3 | 0.5345 | 0.4507 | 0.6905 | 0.1471 | 0.03 | 0.0 |
| 18 | time_penalty=0.5, dynamics=0.3 | 0.5341 | 0.4507 | 0.6895 | 0.1405 | 0.03 | 0.0 |
| 19 | window=0.25, extend_penalty=0.5, time_penalty=2 | 0.5340 | 0.4542 | 0.6821 | 0.1502 | 0.00 | 0.0 |
| 20 | extend_penalty=0.5, time_penalty=2 | 0.5335 | 0.4542 | 0.6822 | 0.1490 | 0.00 | 0.0 |
| 236 | standard DTW | 0.4694 | 0.3971 | 0.6183 | 0.1467 | 0.32 | 0.1 |
| 2368 | euclidean | 0.3738 | 0.2324 | 0.5520 | 0.0413 | 0.08 | 0.0 |
| 2738 | pearson | 0.3730 | 0.2880 | 0.4584 | 0.1297 | 0.01 | 0.1 |

推荐：**window=0.25, time_penalty=2**  
`python3 dtw.py --input <centered.csv> --output outputs/<name> --window 0.25 --time-penalty 2 --scale-floor 0`

每个参数取某个值时，所有含这个值的组合的平均 ARI 和最好 ARI（看哪些值普遍有帮助）：

| 参数 | 值 | 组合数 | 平均 ARI | 最好 ARI |
| --- | ---: | ---: | ---: | ---: |
| window | 0 | 684 | 0.3975 | 0.5453 |
| window | 0.0834 | 684 | 0.3966 | 0.5129 |
| window | 0.1667 | 684 | 0.3975 | 0.5411 |
| window | 0.25 | 684 | 0.3976 | 0.5476 |
| step_penalty | 0 | 684 | 0.4157 | 0.5476 |
| step_penalty | 0.1 | 684 | 0.4072 | 0.5288 |
| step_penalty | 0.5 | 684 | 0.3903 | 0.4540 |
| step_penalty | 2 | 684 | 0.3760 | 0.3810 |
| open_penalty | 0 | 912 | 0.4153 | 0.5476 |
| open_penalty | 0.2 | 912 | 0.4004 | 0.5298 |
| open_penalty | 2 | 912 | 0.3762 | 0.3813 |
| extend_penalty | 0 | 912 | 0.3975 | 0.5476 |
| extend_penalty | 0.05 | 912 | 0.3975 | 0.5447 |
| extend_penalty | 0.5 | 912 | 0.3969 | 0.5412 |
| time_penalty | 0 | 144 | 0.4245 | 0.5298 |
| time_penalty | 0.5 | 864 | 0.4101 | 0.5447 |
| time_penalty | 2 | 864 | 0.3962 | 0.5476 |
| time_penalty | 10 | 864 | 0.3810 | 0.4747 |
| power | 1 | 1296 | 0.3830 | 0.5028 |
| power | 2 | 1296 | 0.4085 | 0.5476 |
| dynamics | 0 | 864 | 0.4030 | 0.5476 |
| dynamics | 0.3 | 864 | 0.3954 | 0.5422 |
| dynamics | 1 | 864 | 0.3889 | 0.5330 |
