# 集成交接记录

交付版本：Member E integration v1，2026-10-26。没有现存 Git 仓库提交号；以 `code_manifest.json` 的 SHA256 锁定交付源码。原目录已有 `Time-series clustering.md` 的分工要求，未提供 A–D 源码；本次没有覆盖任何成员代码。

| 成员 | 需合并/确认内容 | 当前状态 |
|---|---|---|
| A | 正式 simulation.py、生成规则说明；共同确认 preprocessing.py | 已提供可运行接口与补充模拟；成员版本待提供 |
| B | 正式 DTW DP、局部代价、window 及边界测试 | 已提供独立 DP 实现；成员版本待提供 |
| C | 正式聚类、初始化、并列距离约定 | 已提供交替 k-medoids；成员版本待提供 |
| D | 正式 benchmark、独立调参 seed、五因素重复实验及规模测试 | 已提供补充集成 benchmark；正式评测及成员版本待提供 |
| E | GSE3635、Spellman/phase 验证、main.py、运行文档 | 本目录实现和结果 |

同一距离比较应锁定相同输入、预处理、k、优化器、seed 与起点次数。请勿把参考标签或 true_cluster 放入 main.py 输入；验证必须在拟合之后调用。若更换 B/C 实现，需重新运行 analyze_real.py，而不能沿用本次 assignments。

建议组员合并时填写：姓名、源文件版本/提交、Python 版本、依赖、接口变更、测试命令和结果路径。当前没有这些信息，保留待填状态。最终课程整体验收仍须包含 A/D 的正式模拟设计与完整 benchmark；E 的生物学结果不替代它。
