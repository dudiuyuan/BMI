# DTW 网格搜索结果

数据：`../../wpj/predata/数据准备/outputs/simulated/baseline`；scale_floor=0.1；聚类：徐秋凡的 k-medoids，20 个种子取平均。

网格：`{"window": [0, 0.0834, 0.1667, 0.25], "step_penalty": [0, 0.1, 0.5, 2], "open_penalty": [0, 0.2, 2], "extend_penalty": [0, 0.05, 0.5], "time_penalty": [0, 0.5, 2, 10], "power": [1, 2], "dynamics": [0, 0.3, 1]}`

排序：ARI 平均 → ARI 最差种子 → NMI → 真标签轮廓系数 → 开启的参数越少越好 → 越快越好。前几名之间差距很小时，可能只是聚类的随机波动。

### 标签 true_cluster（k=6），共 2738 组

| 排名 | 配置 | ARI 平均 | ARI 最差 | NMI | 真标签轮廓 | 早晚峰同簇率 | 秒 |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | window=0.25, step_penalty=0.1, open_penalty=2, extend_penalty=0.5, time_penalty=10 | 0.9244 | 0.6941 | 0.9416 | 0.5824 | 0.00 | 0.0 |
| 2 | step_penalty=0.1, open_penalty=2, extend_penalty=0.5, time_penalty=10 | 0.9244 | 0.6941 | 0.9416 | 0.5820 | 0.00 | 0.0 |
| 3 | extend_penalty=0.5, time_penalty=10 | 0.9217 | 0.7077 | 0.9492 | 0.5756 | 0.00 | 0.0 |
| 4 | extend_penalty=0.5, time_penalty=2, dynamics=0.3 | 0.9208 | 0.7074 | 0.9502 | 0.5660 | 0.00 | 0.0 |
| 5 | step_penalty=0.1, extend_penalty=0.5, time_penalty=2, dynamics=0.3 | 0.9208 | 0.7088 | 0.9490 | 0.5653 | 0.00 | 0.0 |
| 6 | window=0.25, step_penalty=0.1, extend_penalty=0.5, time_penalty=0.5, dynamics=1 | 0.9203 | 0.7088 | 0.9474 | 0.5565 | 0.00 | 0.0 |
| 7 | window=0.25, extend_penalty=0.5, time_penalty=2, dynamics=0.3 | 0.9202 | 0.7074 | 0.9475 | 0.5710 | 0.00 | 0.0 |
| 8 | window=0.25, step_penalty=0.1, extend_penalty=0.5, time_penalty=2, dynamics=0.3 | 0.9201 | 0.7088 | 0.9431 | 0.5694 | 0.00 | 0.0 |
| 9 | window=0.25, extend_penalty=0.5, time_penalty=10 | 0.9196 | 0.7077 | 0.9432 | 0.5776 | 0.00 | 0.0 |
| 10 | window=0.25, step_penalty=0.1, extend_penalty=0.5, time_penalty=10 | 0.9191 | 0.7098 | 0.9412 | 0.5785 | 0.00 | 0.0 |
| 11 | step_penalty=0.1, extend_penalty=0.5, time_penalty=10 | 0.9191 | 0.7098 | 0.9412 | 0.5771 | 0.00 | 0.0 |
| 12 | window=0.25, step_penalty=0.5, open_penalty=0.2, extend_penalty=0.5, time_penalty=10 | 0.9136 | 0.6964 | 0.9369 | 0.5879 | 0.00 | 0.0 |
| 13 | step_penalty=0.5, open_penalty=0.2, extend_penalty=0.5, time_penalty=10 | 0.9136 | 0.6964 | 0.9369 | 0.5877 | 0.00 | 0.0 |
| 14 | window=0.1667, step_penalty=0.1, open_penalty=2, time_penalty=2 | 0.9126 | 0.6929 | 0.9347 | 0.5803 | 0.00 | 0.0 |
| 15 | window=0.1667, step_penalty=0.1, open_penalty=2, time_penalty=0.5, dynamics=0.3 | 0.9126 | 0.6929 | 0.9347 | 0.5798 | 0.00 | 0.0 |
| 16 | window=0.1667, step_penalty=0.1, open_penalty=2, extend_penalty=0.05, time_penalty=0.5 | 0.9126 | 0.6929 | 0.9347 | 0.5794 | 0.00 | 0.0 |
| 17 | step_penalty=0.5, open_penalty=2, extend_penalty=0.05, time_penalty=0.5 | 0.9125 | 0.6929 | 0.9351 | 0.5468 | 0.00 | 0.0 |
| 18 | window=0.25, step_penalty=0.1, open_penalty=2, extend_penalty=0.5, time_penalty=2, power=1 | 0.9124 | 0.6941 | 0.9356 | 0.5786 | 0.00 | 0.0 |
| 19 | step_penalty=0.1, open_penalty=2, extend_penalty=0.5, time_penalty=2, power=1 | 0.9124 | 0.6941 | 0.9356 | 0.5767 | 0.00 | 0.0 |
| 20 | window=0.1667, open_penalty=2, extend_penalty=0.05, time_penalty=0.5, power=1 | 0.9124 | 0.6929 | 0.9347 | 0.5796 | 0.00 | 0.0 |
| 1936 | euclidean | 0.8525 | 0.6943 | 0.8848 | 0.5830 | 0.00 | 0.0 |
| 2734 | pearson | 0.7324 | 0.7053 | 0.7879 | 0.6556 | 0.00 | 0.1 |
| 2738 | standard DTW | 0.6354 | 0.5227 | 0.7559 | 0.3361 | 0.33 | 0.2 |

推荐：**window=0.25, step_penalty=0.1, open_penalty=2, extend_penalty=0.5, time_penalty=10**  
`python3 dtw.py --input <centered.csv> --output outputs/<name> --window 0.25 --step-penalty 0.1 --open-penalty 2 --extend-penalty 0.5 --time-penalty 10`

每个参数取某个值时，所有含这个值的组合的平均 ARI 和最好 ARI（看哪些值普遍有帮助）：

| 参数 | 值 | 组合数 | 平均 ARI | 最好 ARI |
| --- | ---: | ---: | ---: | ---: |
| window | 0 | 684 | 0.8674 | 0.9244 |
| window | 0.0834 | 684 | 0.8544 | 0.8812 |
| window | 0.1667 | 684 | 0.8718 | 0.9126 |
| window | 0.25 | 684 | 0.8694 | 0.9244 |
| step_penalty | 0 | 684 | 0.8600 | 0.9217 |
| step_penalty | 0.1 | 684 | 0.8621 | 0.9244 |
| step_penalty | 0.5 | 684 | 0.8680 | 0.9136 |
| step_penalty | 2 | 684 | 0.8728 | 0.9118 |
| open_penalty | 0 | 912 | 0.8614 | 0.9217 |
| open_penalty | 0.2 | 912 | 0.8615 | 0.9136 |
| open_penalty | 2 | 912 | 0.8743 | 0.9244 |
| extend_penalty | 0 | 912 | 0.8637 | 0.9126 |
| extend_penalty | 0.05 | 912 | 0.8644 | 0.9126 |
| extend_penalty | 0.5 | 912 | 0.8692 | 0.9244 |
| time_penalty | 0 | 144 | 0.8554 | 0.9111 |
| time_penalty | 0.5 | 864 | 0.8683 | 0.9203 |
| time_penalty | 2 | 864 | 0.8674 | 0.9208 |
| time_penalty | 10 | 864 | 0.8633 | 0.9244 |
| power | 1 | 1296 | 0.8649 | 0.9124 |
| power | 2 | 1296 | 0.8677 | 0.9244 |
| dynamics | 0 | 864 | 0.8688 | 0.9244 |
| dynamics | 0.3 | 864 | 0.8657 | 0.9208 |
| dynamics | 1 | 864 | 0.8645 | 0.9203 |

### 标签 shape_cluster（k=5），共 2738 组

| 排名 | 配置 | ARI 平均 | ARI 最差 | NMI | 真标签轮廓 | 早晚峰同簇率 | 秒 |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | window=0.25, open_penalty=0.2, extend_penalty=0.5, time_penalty=0.5 | 0.6360 | 0.6200 | 0.7900 | 0.2600 | 0.00 | 0.0 |
| 2 | window=0.25, extend_penalty=0.5, time_penalty=0.5, dynamics=1 | 0.6332 | 0.4396 | 0.7915 | 0.3019 | 0.00 | 0.0 |
| 3 | window=0.25, step_penalty=0.1, extend_penalty=0.5, time_penalty=0.5 | 0.6331 | 0.4218 | 0.7931 | 0.2596 | 0.00 | 0.0 |
| 4 | window=0.25, extend_penalty=0.5, time_penalty=2 | 0.6319 | 0.4477 | 0.7915 | 0.2697 | 0.00 | 0.0 |
| 5 | extend_penalty=0.05, time_penalty=0.5, dynamics=1 | 0.6271 | 0.4309 | 0.7826 | 0.2960 | 0.00 | 0.0 |
| 6 | step_penalty=0.1, open_penalty=0.2, time_penalty=0.5, dynamics=1 | 0.6264 | 0.6081 | 0.7744 | 0.3055 | 0.00 | 0.0 |
| 7 | open_penalty=0.2, extend_penalty=0.5 | 0.6263 | 0.6157 | 0.7774 | 0.3234 | 0.00 | 0.0 |
| 8 | window=0.25, open_penalty=0.2, extend_penalty=0.5 | 0.6252 | 0.4045 | 0.7811 | 0.2536 | 0.00 | 0.0 |
| 9 | window=0.25, step_penalty=0.1, extend_penalty=0.05, time_penalty=0.5, power=1 | 0.6251 | 0.4209 | 0.7802 | 0.2473 | 0.00 | 0.0 |
| 10 | window=0.25, open_penalty=0.2, extend_penalty=0.05, time_penalty=0.5, power=1 | 0.6240 | 0.4309 | 0.7813 | 0.2379 | 0.00 | 0.0 |
| 11 | window=0.25, step_penalty=0.1, extend_penalty=0.5 | 0.6240 | 0.4083 | 0.7827 | 0.2532 | 0.00 | 0.0 |
| 12 | step_penalty=0.5, extend_penalty=0.5 | 0.6232 | 0.6021 | 0.7686 | 0.3545 | 0.00 | 0.0 |
| 13 | window=0.25, step_penalty=0.1, extend_penalty=0.5, time_penalty=0.5, dynamics=0.3 | 0.6227 | 0.4218 | 0.7830 | 0.2754 | 0.00 | 0.0 |
| 14 | extend_penalty=0.5, time_penalty=0.5, dynamics=1 | 0.6226 | 0.4396 | 0.7834 | 0.3092 | 0.00 | 0.0 |
| 15 | open_penalty=0.2, extend_penalty=0.05, time_penalty=0.5, dynamics=1 | 0.6223 | 0.4222 | 0.7791 | 0.2979 | 0.00 | 0.0 |
| 16 | step_penalty=0.1, extend_penalty=0.05, time_penalty=0.5, dynamics=1 | 0.6218 | 0.4244 | 0.7732 | 0.3035 | 0.00 | 0.0 |
| 17 | window=0.25, open_penalty=0.2, extend_penalty=0.05, time_penalty=0.5, dynamics=1 | 0.6212 | 0.4309 | 0.7768 | 0.2755 | 0.00 | 0.0 |
| 18 | standard DTW | 0.6209 | 0.4247 | 0.7809 | 0.2898 | 0.65 | 0.2 |
| 19 | open_penalty=0.2, extend_penalty=0.05, time_penalty=2 | 0.6209 | 0.4229 | 0.7794 | 0.2923 | 0.00 | 0.0 |
| 20 | step_penalty=0.5, open_penalty=0.2, extend_penalty=0.05, time_penalty=0.5 | 0.6207 | 0.6015 | 0.7634 | 0.3698 | 0.00 | 0.0 |
| 1419 | euclidean | 0.5829 | 0.3693 | 0.7011 | 0.3905 | 0.00 | 0.0 |
| 2087 | pearson | 0.5825 | 0.4761 | 0.7025 | 0.3758 | 0.00 | 0.1 |

推荐：**window=0.25, open_penalty=0.2, extend_penalty=0.5, time_penalty=0.5**  
`python3 dtw.py --input <centered.csv> --output outputs/<name> --window 0.25 --open-penalty 0.2 --extend-penalty 0.5 --time-penalty 0.5`

每个参数取某个值时，所有含这个值的组合的平均 ARI 和最好 ARI（看哪些值普遍有帮助）：

| 参数 | 值 | 组合数 | 平均 ARI | 最好 ARI |
| --- | ---: | ---: | ---: | ---: |
| window | 0 | 684 | 0.5893 | 0.6271 |
| window | 0.0834 | 684 | 0.5781 | 0.5964 |
| window | 0.1667 | 684 | 0.5838 | 0.6156 |
| window | 0.25 | 684 | 0.5871 | 0.6360 |
| step_penalty | 0 | 684 | 0.5853 | 0.6360 |
| step_penalty | 0.1 | 684 | 0.5857 | 0.6331 |
| step_penalty | 0.5 | 684 | 0.5862 | 0.6232 |
| step_penalty | 2 | 684 | 0.5810 | 0.5940 |
| open_penalty | 0 | 912 | 0.5852 | 0.6332 |
| open_penalty | 0.2 | 912 | 0.5853 | 0.6360 |
| open_penalty | 2 | 912 | 0.5832 | 0.6140 |
| extend_penalty | 0 | 912 | 0.5844 | 0.6264 |
| extend_penalty | 0.05 | 912 | 0.5849 | 0.6271 |
| extend_penalty | 0.5 | 912 | 0.5844 | 0.6360 |
| time_penalty | 0 | 144 | 0.5886 | 0.6263 |
| time_penalty | 0.5 | 864 | 0.5861 | 0.6360 |
| time_penalty | 2 | 864 | 0.5840 | 0.6319 |
| time_penalty | 10 | 864 | 0.5830 | 0.6119 |
| power | 1 | 1296 | 0.5827 | 0.6251 |
| power | 2 | 1296 | 0.5860 | 0.6360 |
| dynamics | 0 | 864 | 0.5854 | 0.6360 |
| dynamics | 0.3 | 864 | 0.5848 | 0.6227 |
| dynamics | 1 | 864 | 0.5828 | 0.6332 |
