# 自动测试

在 dtw 根目录手动执行 `python3 -B -m unittest discover -s tests -v`。

- test_dtw.py：核心递推、边界、穷举路径核对、带下限的幅度归一化、单对与矩阵的一致性、旧实验复现。
- test_diagnostics.py：平均连接、PAM k-medoids 的目标函数与交换检查、评价、标签关联、归一化后诊断的一致性。

测试中的临时文件在临时目录创建并自动清理。日常运行实验和 diagnose 不会自动执行测试。
- test_baselines.py：欧氏、皮尔森距离的手算值和边界，欧氏等于 DTW 对角路径，以及与 DTW 相同的输入、预处理、输出和诊断流程。
- test_xqf_adapter.py：接入徐秋凡 k-medoids 的格式转换与错误检查。
- test_pruning.py：窗口范围写法与旧实现逐项一致；格子跳过不改变距离和路线；LB_Keogh 不超过 DTW；提前终止只在超过上限时发生；整对跳过等于 min(精确距离, 阈值)；命令行记录、复用和诊断。
- test_tuning.py：不需要 Jupyter，直接执行 tune_dtw.ipynb 的代码格；指标手算值，以及在小数据上完整跑一遍并检查输出文件。
- test_modes.py：四种固定变体的参数；每种模式与在 dtw.py 中写出同样参数的结果逐格相同；模式固定的参数不能被覆盖。
