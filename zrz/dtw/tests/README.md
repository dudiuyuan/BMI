# 自动测试

在 dtw 根目录手动执行 `python3 -B -m unittest discover -s tests -v`。

- test_dtw.py：核心递推、边界、穷举路径核对、带下限的幅度归一化、单对与矩阵的一致性、旧实验复现。
- test_diagnostics.py：平均连接、PAM k-medoids 的目标函数与交换检查、评价、标签关联、归一化后诊断的一致性。

测试中的临时文件在临时目录创建并自动清理。日常运行实验和 diagnose 不会自动执行测试。
