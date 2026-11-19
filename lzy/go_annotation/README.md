# Spellman 800 基因的 UniProt GO 长表

本目录独立运行，不改变既有表达分析、聚类结果或用户提供的 GO.csv。仅依赖 Python 3.10+ 标准库，无需 pip 安装。

## 运行

在 PowerShell 中运行（可从任意工作目录执行）：

```powershell
python D:\BMI\project11\go_annotation\fetch_spellman_go.py
```

第一次运行需要访问 `https://rest.uniprot.org`。后续默认复用有 SHA256 校验的缓存。严格离线重建：

```powershell
python D:\BMI\project11\go_annotation\fetch_spellman_go.py --offline
```

获取未来版本的数据时，必须用新的输出目录，避免与旧缓存混合：

```powershell
python D:\BMI\project11\go_annotation\fetch_spellman_go.py --output D:\BMI\project11\go_annotation\results_new
```

默认输入为 `../data/references/spellman800.csv`（相对于脚本所在目录）。结果目录内复制了一份输入快照，因此原项目搬走后仍可这样复现：

```powershell
python fetch_spellman_go.py --input results/spellman800_input.csv --output results --offline
```

## 输出

- `results/GO.csv`：与用户示例完全相同的四列 `GO_id,entry_id,GO_meaning,GO_code`。一行一个唯一 UniProt accession 与 GO 术语组合。
- `results/GO_long.csv`：可追溯主表，一行一个唯一 Spellman 基因、UniProt accession 与 GO 术语组合。前四列同上，并添加 gene_id、Spellman 名称和阶段、GO 类别、证据、参考文献、匹配字段、UniProt 链接与数据库版本。
- `results/gene_uniprot_mapping.csv`：所有匹配关系，一行一个基因与条目组合。
- `results/gene_audit.csv`：严格保留全部 800 个基因，一行一个基因。区分 unmatched（未匹配）、matched_without_GO（匹配但无 GO）、matched_with_GO（匹配且有 GO）。记录一对多匹配。
- `results/genes_needing_review.csv`：未匹配、无 GO、多条目匹配或多个历史基因共享条目的基因；没有问题时只有表头。
- `results/raw/*.json`：UniProt API 原始响应字节；`*.meta.json` 为请求 URL、响应头、抓取时间及 SHA256。
- `results/manifest.json`：输入、脚本、输出 SHA256、检索式、全部请求、UniProt 发布版本及数量统计。

`GO_id` 是连续行号，与示例一致；真正的 GO 标识是 `GO_code`。不同输出的行号不应作为连接键；使用 gene_id、entry_id、GO_code。

## 透明匹配与注释规则

1. 输入必须恰好包含 800 个唯一基因 ID。不是只取 GSE3635 中清洗后匹配的 784 个。
2. 下载酿酒酵母 S288C（taxon 559292）参考蛋白组 UP000002311，分页跟随服务器的 Link，核对总条目数和 accession 唯一性。每页保存原始响应及元信息。
3. 使用 UniProt 的 geneName、synonyms、orderedLocusNames、orfNames 中与 Spellman ORF **完全相等**的名称匹配。禁止模糊匹配、大小写猜测、擅自历史别名替换或跨菌株替代。
4. 对参考蛋白组中未匹配的 ID，追加该菌株的 `gene_exact` 批量检索，再次在返回的基因字段上核对完全相等。所有匹配条目都保留，不按第一个结果静默选择，也不把未审阅条目伪称为 Swiss-Prot。
5. 从 `uniProtKBCrossReferences` 中读取 GO 条目。`GoTerm` 的 C/F/P 前缀分别对应 cellular_component（细胞组分）、molecular_function（分子功能）、biological_process（生物过程）；其余部分写入 GO_meaning。
6. 保留所有证据类型，包括计算推断的 IEA 和实验等证据，不把全部 GO 当作实验验证。`GO_evidence` 保留 UniProt 提供的“证据代码:来源”；`ECO_codes` 和 `references` 保留可用的证据与文献。空白表示该 API 响应未提供，不补造。
7. 同一基因、条目、GO ID 合并重复证据，证据排序后用分号连接。没有 GO 的基因不伪造空 GO 行，而在 800 行审计表中明确记录。
8. 只输出 UniProt 条目返回的 GO 注释，不扩展 GO 祖先，不进行 GO 富集检验，不把 Spellman phase 推导为 GO。
9. 所有响应必须来自同一 UniProt release；缓存逐文件校验。表格排序固定、CSV 为 UTF-8 BOM，适合 Excel 打开。原始 API 仍是英文术语，不自动翻译。

同一 UniProt 条目可能同时登记多个基因位点或历史 ORF 名称。长表保留这些精确名称匹配，在 `spellman_genes_for_entry` 列列出共享该条目的全部输入基因，审计表用 `shared_entry` 标记。GO 是蛋白条目层面的注释，不能把这种多对一关系当作多个基因各有独立实验验证。四列 GO.csv 按 entry_id + GO_code 去重，因此可能比带基因 ID 的 GO_long.csv 行数少。

## 本次结果与验证

实际 UniProt release 为 **2026_04**。800 个输入基因中，743 个有 GO、52 个匹配但无 GO、5 个未匹配。长表 6753 行；按条目与 GO 去重的四列表 6684 行。5 个共享条目涉及 11 个输入基因，已明确标记。未匹配 ID 为 YCL060C、YCL062W、YCLX09W、YCRX05W、YML035C-A，不猜测替代 ID。无 GO 不代表没有生物学功能。

运行验证脚本可逐行核对原始 UniProt 响应、原始 orf800 名单、全部 800 基因去向，以及离线重建前后所有 CSV 的字节级 SHA256 一致性：

```powershell
python D:\BMI\project11\go_annotation\verify_go.py
```

验证结果写入 `results/verification.json`。本节数量描述随附快照，更新数据库后以新的 manifest 为准。

在线重新查询时，数据库可能更新；要逐字节复现本次 CSV，使用本次输入、脚本和完整 raw 缓存进行 `--offline` 重建。这个结果是当前 UniProt 注释，不是 1998 年已有的 GO 注释快照。

## 参考来源

- Spellman PT, Sherlock G, Zhang MQ, et al. (1998). Comprehensive identification of cell cycle-regulated genes of the yeast Saccharomyces cerevisiae by microarray hybridization. Molecular Biology of the Cell 9(12):3273–3297. https://doi.org/10.1091/mbc.9.12.3273
- 基因名单实际来源：Bioconductor yeastCC 1.50.0 的 orf800.txt 与 CellCycle98.txt；项目原始参考文件在 `../data/references/`，包内 orf800.Rd 追溯到原研究 CellCycle98.xls。https://bioconductor.org/packages/3.22/data/experiment/html/yeastCC.html
- UniProt REST API: https://www.uniprot.org/help/api_queries
- S288C 参考蛋白组: https://www.uniprot.org/proteomes/UP000002311

实际数据库版本、请求地址、下载时间和校验和以本次 `manifest.json` 为准。
