# 新数据验证（true_cluster）

新数据：`../simulation_new`；参数组：baseline 上 ARI ≥ 0.9 的 97 组 + 参照。

baseline 排名与新数据排名的相关（Spearman）：**-0.03**。接近 1 表示 baseline 上好的参数在新数据上也好；接近 0 表示 baseline 上的排名在新数据上不成立。

| 新数据排名 | 配置 | 开启参数 | 新数据 ARI | 新数据 ARI 最差 | 新数据 NMI | baseline ARI | 变化 | baseline 排名 |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | window=0.25, open_penalty=0.2, extend_penalty=0.5, time_penalty=0.5, dynamics=1 | 5 | 0.6628 | 0.6191 | 0.7192 | 0.9047 | -0.2419 | 88 |
| 2 | window=0.25, extend_penalty=0.5, time_penalty=0.5, dynamics=1 | 4 | 0.6609 | 0.6288 | 0.7177 | 0.9117 | -0.2509 | 28 |
| 3 | window=0.25, extend_penalty=0.5, time_penalty=2, dynamics=0.3 | 4 | 0.6607 | 0.6256 | 0.7160 | 0.9202 | -0.2596 | 7 |
| 4 | window=0.25, extend_penalty=0.05, time_penalty=2, dynamics=0.3 | 4 | 0.6596 | 0.6315 | 0.7210 | 0.9066 | -0.2470 | 87 |
| 5 | window=0.25, step_penalty=0.1, extend_penalty=0.5, time_penalty=0.5, dynamics=1 | 5 | 0.6594 | 0.6196 | 0.7163 | 0.9203 | -0.2609 | 6 |
| 6 | extend_penalty=0.5, time_penalty=2, dynamics=0.3 | 3 | 0.6575 | 0.6256 | 0.7142 | 0.9208 | -0.2633 | 4 |
| 7 | window=0.25, step_penalty=0.1, extend_penalty=0.5, time_penalty=2, dynamics=0.3 | 5 | 0.6558 | 0.6196 | 0.7140 | 0.9201 | -0.2643 | 8 |
| 8 | window=0.25, open_penalty=0.2, extend_penalty=0.5, time_penalty=2, dynamics=0.3 | 5 | 0.6557 | 0.6196 | 0.7135 | 0.9020 | -0.2464 | 90 |
| 9 | open_penalty=0.2, extend_penalty=0.5, time_penalty=2, dynamics=0.3 | 4 | 0.6545 | 0.6196 | 0.7119 | 0.9040 | -0.2495 | 89 |
| 10 | step_penalty=0.1, extend_penalty=0.5, time_penalty=2, dynamics=0.3 | 4 | 0.6539 | 0.6196 | 0.7117 | 0.9208 | -0.2669 | 5 |
| 11 | window=0.25, extend_penalty=0.5, time_penalty=10 | 3 | 0.6522 | 0.6203 | 0.7096 | 0.9196 | -0.2674 | 9 |
| 12 | extend_penalty=0.5, time_penalty=10 | 2 | 0.6511 | 0.6203 | 0.7085 | 0.9217 | -0.2706 | 3 |
| 13 | window=0.1667, step_penalty=0.5, extend_penalty=0.5, time_penalty=2 | 4 | 0.6475 | 0.6054 | 0.7097 | 0.9003 | -0.2528 | 97 |
| 14 | window=0.25, step_penalty=0.1, extend_penalty=0.5, time_penalty=10 | 4 | 0.6469 | 0.6208 | 0.7048 | 0.9191 | -0.2721 | 10 |
| 15 | step_penalty=0.1, extend_penalty=0.5, time_penalty=10 | 3 | 0.6469 | 0.6208 | 0.7048 | 0.9191 | -0.2721 | 11 |
| 16 | window=0.1667, step_penalty=0.5, open_penalty=0.2, extend_penalty=0.05, time_penalty=0.5, dynamics=0.3 | 6 | 0.6406 | 0.6017 | 0.7016 | 0.9009 | -0.2604 | 93 |
| 17 | window=0.1667, step_penalty=0.5, open_penalty=2, extend_penalty=0.5 | 4 | 0.6390 | 0.6165 | 0.7021 | 0.9107 | -0.2717 | 44 |
| 18 | window=0.1667, step_penalty=0.5, open_penalty=2, extend_penalty=0.5, time_penalty=0.5 | 5 | 0.6388 | 0.6165 | 0.7020 | 0.9107 | -0.2720 | 43 |
| 19 | window=0.25, step_penalty=0.5, open_penalty=2, extend_penalty=0.5, time_penalty=2 | 5 | 0.6375 | 0.6170 | 0.7035 | 0.9094 | -0.2719 | 78 |
| 20 | window=0.25, step_penalty=0.5, open_penalty=2, time_penalty=0.5, dynamics=1 | 5 | 0.6369 | 0.6116 | 0.6985 | 0.9090 | -0.2721 | 83 |
| 21 | step_penalty=0.5, open_penalty=2, extend_penalty=0.5, time_penalty=2 | 4 | 0.6368 | 0.6170 | 0.7030 | 0.9094 | -0.2726 | 79 |
| 22 | window=0.25, step_penalty=0.5, open_penalty=2, extend_penalty=0.5, time_penalty=0.5 | 5 | 0.6367 | 0.6165 | 0.7007 | 0.9098 | -0.2730 | 70 |
| 23 | window=0.25, step_penalty=0.5, open_penalty=2, extend_penalty=0.5 | 4 | 0.6366 | 0.6167 | 0.7003 | 0.9098 | -0.2732 | 71 |
| 24 | step_penalty=0.5, open_penalty=2, extend_penalty=0.05, time_penalty=0.5, dynamics=1 | 5 | 0.6363 | 0.6116 | 0.6988 | 0.9100 | -0.2737 | 65 |
| 25 | window=0.25, step_penalty=0.5, open_penalty=2, extend_penalty=0.5, time_penalty=0.5, dynamics=0.3 | 6 | 0.6362 | 0.6170 | 0.7003 | 0.9093 | -0.2731 | 80 |
| 26 | step_penalty=0.5, open_penalty=2, extend_penalty=0.5, time_penalty=0.5, dynamics=0.3 | 5 | 0.6351 | 0.6170 | 0.6995 | 0.9093 | -0.2742 | 81 |
| 27 | step_penalty=0.5, open_penalty=2, time_penalty=0.5, dynamics=1 | 4 | 0.6347 | 0.6116 | 0.6972 | 0.9089 | -0.2742 | 86 |
| 28 | window=0.25, step_penalty=0.5, open_penalty=2 | 3 | 0.6345 | 0.6002 | 0.6970 | 0.9092 | -0.2747 | 82 |
| 29 | window=0.25, step_penalty=0.5, open_penalty=2, extend_penalty=0.05, time_penalty=0.5, dynamics=1 | 6 | 0.6343 | 0.6116 | 0.6976 | 0.9096 | -0.2753 | 77 |
| 30 | window=0.1667, step_penalty=0.1, open_penalty=2, time_penalty=2 | 4 | 0.6329 | 0.5982 | 0.6931 | 0.9126 | -0.2797 | 14 |
| 31 | window=0.1667, step_penalty=0.5, open_penalty=2, extend_penalty=0.05, time_penalty=0.5, dynamics=1 | 6 | 0.6323 | 0.6116 | 0.6943 | 0.9107 | -0.2784 | 45 |
| 32 | window=0.25, step_penalty=0.5, open_penalty=2, time_penalty=0.5 | 4 | 0.6322 | 0.6000 | 0.6942 | 0.9099 | -0.2777 | 68 |
| 33 | window=0.1667, step_penalty=0.1, open_penalty=2, extend_penalty=0.05, time_penalty=0.5, dynamics=0.3 | 6 | 0.6322 | 0.5993 | 0.6937 | 0.9096 | -0.2775 | 72 |
| 34 | window=0.1667, step_penalty=0.5, open_penalty=2, extend_penalty=0.05, time_penalty=0.5, power=1 | 5 | 0.6315 | 0.6116 | 0.6924 | 0.9107 | -0.2792 | 46 |
| 35 | window=0.1667, open_penalty=2, extend_penalty=0.05, time_penalty=0.5, power=1 | 4 | 0.6313 | 0.5982 | 0.6934 | 0.9124 | -0.2811 | 20 |
| 36 | window=0.1667, step_penalty=0.1, open_penalty=2, time_penalty=0.5, dynamics=0.3 | 5 | 0.6307 | 0.5993 | 0.6923 | 0.9126 | -0.2820 | 15 |
| 37 | window=0.1667, step_penalty=0.5, open_penalty=2, time_penalty=0.5, dynamics=1 | 5 | 0.6305 | 0.6116 | 0.6918 | 0.9101 | -0.2796 | 59 |
| 38 | window=0.25, step_penalty=0.5, open_penalty=2, time_penalty=2, dynamics=0.3 | 5 | 0.6302 | 0.6121 | 0.6938 | 0.9096 | -0.2794 | 75 |
| 39 | step_penalty=0.5, open_penalty=2, time_penalty=2, dynamics=0.3 | 4 | 0.6302 | 0.6121 | 0.6938 | 0.9096 | -0.2794 | 76 |
| 40 | window=0.1667, step_penalty=0.1, open_penalty=2, extend_penalty=0.05, time_penalty=0.5 | 5 | 0.6302 | 0.5993 | 0.6930 | 0.9126 | -0.2824 | 16 |
| 41 | window=0.1667, step_penalty=0.5, open_penalty=2, extend_penalty=0.05, time_penalty=0.5 | 5 | 0.6302 | 0.5966 | 0.6926 | 0.9103 | -0.2801 | 56 |
| 42 | window=0.1667, step_penalty=0.1, open_penalty=2, time_penalty=0.5, power=1 | 4 | 0.6298 | 0.5993 | 0.6920 | 0.9122 | -0.2824 | 21 |
| 43 | window=0.1667, step_penalty=0.5, open_penalty=2, time_penalty=2 | 4 | 0.6297 | 0.6002 | 0.6911 | 0.9104 | -0.2807 | 52 |
| 44 | window=0.1667, step_penalty=0.5, open_penalty=2 | 3 | 0.6289 | 0.5988 | 0.6916 | 0.9111 | -0.2822 | 35 |
| 45 | window=0.1667, step_penalty=0.1, open_penalty=2, extend_penalty=0.05, time_penalty=2 | 5 | 0.6287 | 0.5982 | 0.6904 | 0.9103 | -0.2816 | 54 |
| 46 | window=0.1667, step_penalty=0.5, open_penalty=2, extend_penalty=0.05, time_penalty=0.5, dynamics=0.3 | 6 | 0.6283 | 0.6002 | 0.6893 | 0.9104 | -0.2821 | 50 |
| 47 | window=0.1667, step_penalty=0.5, open_penalty=2, time_penalty=0.5, dynamics=0.3 | 5 | 0.6281 | 0.5967 | 0.6895 | 0.9103 | -0.2821 | 55 |
| 48 | window=0.1667, step_penalty=0.5, open_penalty=2, time_penalty=0.5, power=1 | 4 | 0.6280 | 0.5962 | 0.6883 | 0.9101 | -0.2821 | 60 |
| 49 | window=0.1667, step_penalty=0.5, open_penalty=2, extend_penalty=0.05 | 4 | 0.6265 | 0.5966 | 0.6887 | 0.9111 | -0.2846 | 34 |
| 50 | window=0.1667, step_penalty=0.5, open_penalty=2, time_penalty=0.5 | 4 | 0.6265 | 0.5966 | 0.6886 | 0.9110 | -0.2845 | 40 |
| 51 | window=0.1667, step_penalty=0.5, open_penalty=2, extend_penalty=0.05, time_penalty=2 | 5 | 0.6264 | 0.5967 | 0.6879 | 0.9104 | -0.2840 | 51 |
| 52 | window=0.25, step_penalty=0.1, open_penalty=2, time_penalty=2, power=1 | 4 | 0.6257 | 0.5947 | 0.6874 | 0.9111 | -0.2853 | 36 |
| 53 | step_penalty=0.5, extend_penalty=0.05, time_penalty=2, power=1 | 3 | 0.6257 | 0.5984 | 0.6876 | 0.9009 | -0.2752 | 94 |
| 54 | step_penalty=0.5, open_penalty=0.2, extend_penalty=0.5, time_penalty=10 | 4 | 0.6256 | 0.5971 | 0.6864 | 0.9136 | -0.2880 | 13 |
| 55 | window=0.1667, step_penalty=2, extend_penalty=0.5, time_penalty=0.5, power=1 | 4 | 0.6249 | 0.6069 | 0.6874 | 0.9111 | -0.2862 | 37 |
| 56 | window=0.25, step_penalty=0.5, open_penalty=0.2, extend_penalty=0.05, time_penalty=0.5, power=1, dynamics=0.3 | 6 | 0.6245 | 0.5971 | 0.6861 | 0.9008 | -0.2763 | 95 |
| 57 | window=0.25, step_penalty=0.5, open_penalty=2, extend_penalty=0.05, time_penalty=2, dynamics=0.3 | 6 | 0.6245 | 0.6133 | 0.6877 | 0.9096 | -0.2851 | 73 |
| 58 | step_penalty=0.5, open_penalty=2, extend_penalty=0.05, time_penalty=2, dynamics=0.3 | 5 | 0.6245 | 0.6133 | 0.6877 | 0.9096 | -0.2851 | 74 |
| 59 | step_penalty=0.1, open_penalty=2, extend_penalty=0.5, time_penalty=10 | 4 | 0.6244 | 0.6100 | 0.6867 | 0.9244 | -0.3000 | 2 |
| 60 | window=0.1667, step_penalty=0.1, open_penalty=2, extend_penalty=0.5, time_penalty=10 | 5 | 0.6238 | 0.6100 | 0.6871 | 0.9099 | -0.2861 | 67 |
| 61 | window=0.25, step_penalty=2, extend_penalty=0.5, time_penalty=0.5, power=1 | 4 | 0.6237 | 0.6069 | 0.6851 | 0.9111 | -0.2874 | 38 |
| 62 | step_penalty=2, extend_penalty=0.5, time_penalty=0.5, power=1 | 3 | 0.6237 | 0.6069 | 0.6851 | 0.9111 | -0.2874 | 39 |
| 63 | window=0.25, step_penalty=0.5, open_penalty=0.2, extend_penalty=0.5, time_penalty=10 | 5 | 0.6235 | 0.5971 | 0.6850 | 0.9136 | -0.2901 | 12 |
| 64 | window=0.25, step_penalty=0.1, open_penalty=2, extend_penalty=0.5, time_penalty=10 | 5 | 0.6235 | 0.6100 | 0.6858 | 0.9244 | -0.3009 | 1 |
| 65 | step_penalty=2, extend_penalty=0.5, time_penalty=0.5 | 3 | 0.6234 | 0.6159 | 0.6867 | 0.9118 | -0.2884 | 27 |
| 66 | window=0.25, step_penalty=2, extend_penalty=0.5, time_penalty=2 | 4 | 0.6231 | 0.6159 | 0.6855 | 0.9118 | -0.2887 | 22 |
| 67 | window=0.25, step_penalty=2, extend_penalty=0.5, time_penalty=0.5, dynamics=0.3 | 5 | 0.6231 | 0.6159 | 0.6855 | 0.9118 | -0.2887 | 24 |
| 68 | step_penalty=2, extend_penalty=0.5, time_penalty=2 | 3 | 0.6231 | 0.6159 | 0.6855 | 0.9118 | -0.2887 | 23 |
| 69 | step_penalty=2, extend_penalty=0.5, time_penalty=0.5, dynamics=0.3 | 4 | 0.6231 | 0.6159 | 0.6855 | 0.9118 | -0.2887 | 25 |
| 70 | window=0.25, step_penalty=2, extend_penalty=0.5, time_penalty=0.5 | 4 | 0.6229 | 0.6159 | 0.6851 | 0.9118 | -0.2889 | 26 |
| 71 | window=0.25, step_penalty=0.5, open_penalty=0.2, time_penalty=0.5, power=1, dynamics=0.3 | 5 | 0.6218 | 0.5971 | 0.6853 | 0.9012 | -0.2794 | 92 |
| 72 | window=0.1667, step_penalty=0.1, open_penalty=2, time_penalty=2, power=1 | 4 | 0.6218 | 0.5947 | 0.6841 | 0.9103 | -0.2886 | 53 |
| 73 | step_penalty=0.5, open_penalty=2, extend_penalty=0.05, time_penalty=0.5 | 4 | 0.6210 | 0.4204 | 0.6882 | 0.9125 | -0.2914 | 17 |
| 74 | window=0.1667, step_penalty=0.1, open_penalty=2, extend_penalty=0.05, time_penalty=2, power=1 | 5 | 0.6199 | 0.5947 | 0.6820 | 0.9101 | -0.2902 | 63 |
| 75 | window=0.1667, open_penalty=2, extend_penalty=0.5, time_penalty=10 | 4 | 0.6195 | 0.5958 | 0.6837 | 0.9101 | -0.2906 | 64 |
| 76 | window=0.25, step_penalty=0.5, open_penalty=0.2, extend_penalty=0.5, time_penalty=0.5, power=1, dynamics=0.3 | 6 | 0.6193 | 0.5971 | 0.6816 | 0.9015 | -0.2822 | 91 |
| 77 | step_penalty=0.5, open_penalty=0.2, extend_penalty=0.5, time_penalty=0.5, power=1, dynamics=0.3 | 5 | 0.6193 | 0.5971 | 0.6816 | 0.9006 | -0.2813 | 96 |
| 78 | window=0.25, step_penalty=0.1, open_penalty=2, extend_penalty=0.05, time_penalty=2, power=1 | 5 | 0.6192 | 0.5947 | 0.6813 | 0.9112 | -0.2920 | 33 |
| 79 | step_penalty=0.5, open_penalty=2, extend_penalty=0.05, time_penalty=10 | 4 | 0.6183 | 0.6155 | 0.6822 | 0.9103 | -0.2920 | 58 |
| 80 | window=0.25, step_penalty=0.5, open_penalty=2, extend_penalty=0.05, time_penalty=10 | 5 | 0.6181 | 0.6155 | 0.6821 | 0.9103 | -0.2921 | 57 |
| 81 | window=0.25, step_penalty=0.5, open_penalty=2, time_penalty=10 | 4 | 0.6156 | 0.6049 | 0.6786 | 0.9101 | -0.2944 | 61 |
| 82 | step_penalty=0.5, open_penalty=2, time_penalty=10 | 3 | 0.6156 | 0.6049 | 0.6786 | 0.9101 | -0.2944 | 62 |
| 83 | step_penalty=2, open_penalty=0.2, extend_penalty=0.5, time_penalty=0.5, dynamics=0.3 | 5 | 0.6154 | 0.6009 | 0.6791 | 0.9117 | -0.2963 | 29 |
| 84 | window=0.25, step_penalty=2, open_penalty=0.2, extend_penalty=0.5, time_penalty=0.5, dynamics=0.3 | 6 | 0.6151 | 0.6009 | 0.6787 | 0.9117 | -0.2965 | 30 |
| 85 | window=0.1667, open_penalty=2, extend_penalty=0.5, time_penalty=0.5, power=1, dynamics=0.3 | 5 | 0.6140 | 0.5999 | 0.6757 | 0.9099 | -0.2960 | 66 |
| 86 | window=0.1667, open_penalty=2, extend_penalty=0.5, time_penalty=2, power=1 | 4 | 0.6139 | 0.6025 | 0.6763 | 0.9099 | -0.2959 | 69 |
| 87 | step_penalty=0.1, open_penalty=2, extend_penalty=0.5, time_penalty=2, power=1 | 4 | 0.6133 | 0.5955 | 0.6770 | 0.9124 | -0.2992 | 19 |
| 88 | step_penalty=0.1, open_penalty=2, extend_penalty=0.5, time_penalty=0.5, power=1, dynamics=0.3 | 5 | 0.6132 | 0.5921 | 0.6769 | 0.9090 | -0.2958 | 85 |
| 89 | window=0.25, step_penalty=0.1, open_penalty=2, extend_penalty=0.5, time_penalty=0.5, power=1, dynamics=0.3 | 6 | 0.6127 | 0.5921 | 0.6768 | 0.9090 | -0.2963 | 84 |
| 90 | window=0.25, step_penalty=0.1, open_penalty=2, extend_penalty=0.5, time_penalty=2, power=1 | 5 | 0.6126 | 0.5955 | 0.6762 | 0.9124 | -0.2998 | 18 |
| 91 | window=0.1667, step_penalty=2, extend_penalty=0.5, time_penalty=10 | 4 | 0.6087 | 0.5925 | 0.6704 | 0.9105 | -0.3018 | 47 |
| 92 | window=0.25, step_penalty=2, extend_penalty=0.5, time_penalty=10 | 4 | 0.6087 | 0.5925 | 0.6704 | 0.9105 | -0.3018 | 49 |
| 93 | step_penalty=2, extend_penalty=0.5, time_penalty=10 | 3 | 0.6087 | 0.5925 | 0.6704 | 0.9105 | -0.3018 | 48 |
| 94 | open_penalty=2, time_penalty=0.5, power=1, dynamics=1 | 3 | 0.5890 | 0.5574 | 0.6547 | 0.9108 | -0.3218 | 42 |
| 95 | window=0.25, open_penalty=2, time_penalty=0.5, power=1, dynamics=1 | 4 | 0.5868 | 0.5574 | 0.6526 | 0.9108 | -0.3240 | 41 |
| 96 | open_penalty=2, extend_penalty=0.05, time_penalty=0.5, power=1, dynamics=1 | 4 | 0.5863 | 0.5574 | 0.6520 | 0.9113 | -0.3250 | 32 |
| 97 | window=0.25, open_penalty=2, extend_penalty=0.05, time_penalty=0.5, power=1, dynamics=1 | 5 | 0.5841 | 0.5574 | 0.6500 | 0.9113 | -0.3272 | 31 |
| 98 | standard DTW | 0 | 0.5161 | 0.4165 | 0.6256 | 0.6354 | -0.1192 | 2738 |
| 99 | pearson | 0 | 0.4973 | 0.4649 | 0.5664 | 0.7324 | -0.2352 | 2734 |
| 100 | euclidean | 0 | 0.4881 | 0.4342 | 0.5619 | 0.8525 | -0.3645 | 1936 |

新数据上最好的 DTW 参数组：**window=0.25, open_penalty=0.2, extend_penalty=0.5, time_penalty=0.5, dynamics=1**  
`python3 dtw.py --input <centered.csv> --output outputs/<name> --window 0.25 --open-penalty 0.2 --extend-penalty 0.5 --time-penalty 0.5 --dynamics 1`
