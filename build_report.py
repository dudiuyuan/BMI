"""Generate the deliverable summary directly from reproducible output tables."""
import hashlib
import json
import platform
from pathlib import Path
import pandas as pd


def markdown_table(frame):
    cols = list(frame.columns)
    lines = ['| '+' | '.join(map(str, cols))+' |', '| '+' | '.join(['---']*len(cols))+' |']
    for row in frame.itertuples(index=False, name=None):
        cells = [f'{x:.4g}' if isinstance(x, float) else str(x) for x in row]
        lines.append('| '+' | '.join(cells)+' |')
    return '\n'.join(lines)


def report():
    root = Path(__file__).resolve().parent
    real = root/'results/biological_validation'
    meta = json.loads((root/'data/gse3635/metadata.json').read_text())
    params = json.loads((real/'selected_parameters.json').read_text())
    summary = json.loads((real/'dtw/validation_summary.json').read_text())
    enrichment = pd.read_csv(real/'dtw/cluster_enrichment.csv')
    phase = pd.read_csv(real/'dtw/phase_enrichment.csv')
    selection = pd.read_csv(real/'selection.csv')
    benchmark = pd.read_csv(root/'results/benchmarks/benchmark_runs.csv')
    q = meta['preprocessing']
    references = pd.read_csv(root/'data/references/spellman800.csv')
    samples = pd.read_csv(root/'data/gse3635/samples.csv')
    rawqc = pd.read_csv(root/'data/gse3635/raw_gpr_qc.csv')
    missing_qc = pd.read_csv(root/'data/gse3635/sample_qc.csv')
    # Recover matched and dropped reference genes by exact historical ID.
    prior = pd.read_csv(root/'data/gse3635/gene_log10ratio_with_missing.csv')
    reference_before = len(set(prior.gene_id) & set(references.gene_id))
    phase_sig = phase[phase.q_value_bh < .05][['cluster_id','phase','phase_hits','fraction_among_matched','q_value_bh']]
    comparisons = []
    for method in ['dtw','euclidean','pearson']:
        cfg = json.loads((real/method/'run_config.json').read_text())
        val = json.loads((real/method/'validation_summary.json').read_text())
        comparisons.append(dict(method=method, k=cfg['k'], window=cfg['window'],
                                enriched_clusters=','.join(map(str,val['significant_clusters'])),
                                distance_seconds=cfg['distance_seconds'], total_seconds=cfg['total_seconds']))
    # Targeted descriptive examples; published mechanism is not inferred from these labels.
    annotation = pd.read_csv(root/'data/gse3635/platform_annotation.csv').rename(columns={'ORF':'gene_id','Gene':'symbol'})
    example_symbols = ['MCM1','YOX1','YHP1','MCM2','MCM3','MCM4','MCM5','MCM6','MCM7','CDC46','CDC47','CDC54','CDC6','SWI4','CLN3']
    examples = pd.read_csv(real/'dtw/peak_times.csv').merge(annotation[['gene_id','symbol']].drop_duplicates('gene_id'), on='gene_id', suffixes=('_spellman',''))
    # GPL1914 leaves this gene symbol blank; SGD explicitly supplies its identity.
    examples.loc[examples.gene_id == 'YDR451C', 'symbol'] = 'YHP1 (SGD)'
    examples = examples[examples.symbol.isin(example_symbols) | (examples.gene_id == 'YDR451C')][['gene_id','symbol','cluster_id','phase','observed_peak_time']]
    examples.to_csv(real/'dtw/mechanism_examples.csv', index=False)
    base = benchmark[(benchmark.scaling == 'zscore') & (benchmark.level == 0) & (benchmark.challenge == 'shift')]
    worst = benchmark[(benchmark.scaling == 'zscore') & (benchmark.challenge == 'shift') & (benchmark.level == .2)]
    bench_compare = pd.concat([base.assign(condition='shift=0'), worst.assign(condition='shift=0.2')])
    bench_compare = bench_compare.groupby(['condition','method','window']).ari.agg(['mean','std']).reset_index()
    checks_file = root/'results/acceptance.json'
    checks = json.loads(checks_file.read_text()) if checks_file.exists() else {'status':'尚未生成验收记录'}
    compact_checks = {k: v for k, v in checks.items() if k != 'commands'}
    lines = f'''# Project 11 · Member E 工作总结

交付日期：2026-10-26。分析对象：GSE3635，全部清洗后保留基因。项目目录：`D:\\BMI\\project11`。

## 1. 交付结论

已完成官方 metadata 与样本时间核查、原始 GPR 留存、GPL1914 映射、统一预处理、三种距离的相同 k-medoids 聚类、medoid/heatmap、Spellman 800 基因及 phase 验证、通用 CLI 和隔离环境运行验证。

本次主分析保留 **{q['retained_genes']:,} 个基因 × 13 时间点**；无监督搜索选择 **k={params['k']}，DTW window={params['window']}**（最多匹配相差 {params['window']*10} 分钟的采样点），固定 seed=11。800 个参考基因中匹配 **{summary['matched']} 个**，输入背景比例 **{summary['background_fraction']:.2%}**。BH 校正后 q<0.05 的已知周期基因富集簇：**{summary['significant_clusters']}**。

这些结果属于**真实数据的生物学支持**，不构成真实数据聚类“准确率”。模拟 ARI/NMI 在独立结果目录报告，真实富集不替代模拟 benchmark。

## 2. 来源与数据含义

- [GEO GSE3635](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE3635)：α-factor 同步后每 10 分钟采样，共 0–120 分钟 13 点；GEO 明确说明无重复；芯片平台 [GPL1914](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL1914)。
- Cy5 为同步时间点样本，Cy3 为非同步培养物参考；这是双通道 DNA/cDNA microarray，不是 RNA-seq。
- [GSM81064 处理说明](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM81064) 为 Rosetta Resolver normalized log ratio。程序对全样本 VALUE 与 log10(Ratio) 数值核对，最大残差 `{meta['max_log10_ratio_residual']:.3g}`，支持底数为 10 的推断；不把 raw GPR 的 Log Ratio 列混作 GEO VALUE。
- GEO series matrix 与独立 SOFT sample table 的 VALUE 最大误差 `{meta['series_matrix_max_error']}`。主分析采用已提交的归一化 VALUE，没有声称重现专有 Rosetta Resolver 的原始归一化过程。
- 用户提供 `GSE3635_RAW.tar` 的原始下载日期未知，已如实记录；本次补充下载日期 2026-10-26，文件大小与 SHA256 见 `data/gse3635/source_manifest.json`。

原始 GPR 共 {int(rawqc.spots.sum()):,} 个 spot 记录。原始 spot Flags 与 GEO 已汇总的 Flagged? 处于不同处理层级，分别保存在 `raw_gpr_qc.csv` 和 `geo_quality_flags.csv`；本次以 GEO Flagged? 屏蔽已归一化矩阵，未重复套用未经验证的 GPR 阈值。

### 时间点和样本顺序

{markdown_table(samples[['sample_id','time_minutes','source_scan_date']])}

样本标题中的扫描日期在 0–30 分钟与 40–120 分钟之间分成两批，提示潜在技术批次与时间混杂。不能凭这些日期认定发生了批次效应，也不能在没有重复、且与时间混杂的条件下直接回归掉批次。

## 3. 清洗与 ID provenance

GEO 提交矩阵有 {meta['raw_probe_rows']} 行 probe/汇总 ID，经 GPL1914 的 ID→ORF 对应后排除 {meta['unmapped_or_control_rows']} 个控制或非有效 ORF 行，得到 {q['unique_genes']} 个唯一基因。本数据映射后重复基因数为 {meta['duplicate_genes_before_aggregation']}；通用规则仍规定多探针/重复 ID 逐时间点取中位数。平台中的历史 ORF、Gene、SGDID 原样保留，不猜测旧名与新名的对应。

共屏蔽 {meta['quality_flags_masked']} 个 GEO Flagged?=TRUE 值。表达表另有 {int(missing_qc.missing_values.sum())} 个原始缺失（与标记可能重叠，不能简单相加）。单基因缺失比例大于 20% 或有效点少于 2 时剔除，本次剔除 {q['dropped_genes']} 个，保留 {q['retained_genes']} 个；对保留基因的 {q['imputed_cells']} 个缺失按实际分钟线性插值，端点采用最近观测。所有去留保存在 `gene_filter_audit.csv`。

同时提供 `expression_log10ratio.csv` 和 `expression_zscore.csv`。主分析采用 gene-wise z-score，目的是比较轨迹形状；它消除基因均值和幅度，会放大小幅噪声，因此不能用这套结果声称保留绝对表达量或准确恢复幅度。`--scaling none` 可保留幅度对照。模拟与真实输入都调用同一 `preprocessing.py`。

## 4. 算法与参数选择

DTW 核心 DP 自实现：平方差局部代价，累计代价从 up/left/diagonal 取最小，终点累计后开方，不按路径长度归一化；Numba 只负责编译加速。窗口为索引带 `|i-j|<=w`，w=0 等同 Euclidean。Pearson 使用 1-r；常值与非常值约定距离 1，两个常值距离 0。k-medoids 为 5 个固定种子派生起点的交替分配/簇内更新，最多 100 次迭代；不是全局最优保证。

在全部 {q['retained_genes']} 个基因上分别计算 window∈{{0,1,2}}、k∈{{3,4,5,6,7,8}} 的 18 组候选；在固定 seed=11 抽取的 800 个基因上计算 silhouette，要求每簇至少为输入 1%（向上取整，且不少于 5），取分数最大者，并列优先较小窗口、较小 k。评分样本 ID 保存在 `selection_score_genes.csv`。**调参过程不读取 Spellman 名单、phase 或富集分数。**

{markdown_table(selection[['window','k','silhouette','minimum_cluster_size']])}

选中 sampled silhouette={params['silhouette']:.4f}，整体分离度偏弱。不同窗口改变距离几何，跨窗口 silhouette 比较只是描述性启发式，不是生物学最优证明；候选 k 范围有限，k={params['k']} 也不是宣称真实存在 {params['k']} 个细胞周期阶段。window={params['window']} 保留了有限时间弹性，不能据此消除同步误差和技术混杂。

三种距离在相同输入、标准化、k、seed 和优化器下运行：

{markdown_table(pd.DataFrame(comparisons))}

耗时来自本机实际运行，受首次编译、并行进程与系统负载影响，不宜直接当作严谨性能比较；正式效率实验应由 D 在控制环境中完成。每个方法目录都有全部基因 assignments、真实 medoid 表、完整 heatmap 和配置。

![DTW medoids](results/biological_validation/dtw/medoid_curves.png)

![DTW heatmap](results/biological_validation/dtw/heatmap.png)

每簇轨迹图至多展示 100 条均匀抽取的基因曲线，红线为实际 medoid；heatmap 展示全部基因。颜色范围按绝对值 99% 分位数裁剪仅用于显示，未更改聚类数值。

## 5. Spellman 参考集合和统计验证

参考来源为 [Bioconductor yeastCC](https://bioconductor.org/packages/3.22/data/experiment/html/yeastCC.html) 官方包 v1.50.0 中 `vignettes/orf800.txt` 与 `CellCycle98.txt`，原始文件已保留。`orf800.Rd` 明确追溯至 Spellman 的 CellCycle98.xls；该名单仅作聚类后参考，不作训练标签。

解析确认 800 个唯一 ORF。原始表 YEL077C 的 ORF 单元格为空，但 YPD 与 SGD 两列均为 YEL077C 且存在于 orf800，因此使用这一明确一致的 ID 回填；`spellman800.csv` 的 mapping_rule 保留了此特例。未对无法追溯的历史基因别名作猜测。

输入清洗前可匹配 {reference_before}/800，清洗后 {summary['matched']}/800，覆盖率 {summary['reference_coverage']:.2%}。未匹配名单在 `reference_not_in_input.csv`。背景严格使用**实际进入聚类的全部基因**（N={summary['input_genes']}），其中参考命中 K={summary['matched']}，不使用全酵母基因组或仅800基因作为周期基因富集背景。

对大小为 n、命中 a 的簇，2×2 表为 `[[a,n-a],[K-a,N-n-K+a]]`，单侧 Fisher exact 检验 enrichment，再在本次方法的所有 cluster 上作 BH；原始 p、q、fold enrichment、odds ratio 和四格计数都保存。

{markdown_table(enrichment[['cluster_id','cluster_size','known_count','known_fraction','fold_enrichment','p_value','q_value_bh']])}

q<0.05 表示相对输入背景的参考集合富集，不能换算成聚类准确率；未列入 Spellman 的基因也不等于已证实“非周期”。检验依赖基因集合层面的抽样假设，共表达基因非完全独立，故 p/q 作为探索性支持解释。三种方法分别校正，跨方法比较仍属探索性，未宣称全局跨方法 FDR。

## 6. 阶段验证和峰值时序

保留原始五类标签：G1=300、S=71、S/G2=121、G2/M=195、M/G1=113，**没有擅自改成四类 G1/S/G2/M**。phase 来自历史 CellCycle98 的 Peak 列，不把其数值 Phase 列或 phase order 当作 GSE3635 的分钟数。

阶段富集以输入中有 Spellman 阶段注释的 {summary['matched']} 个基因为背景；对全部 cluster×phase 检验作为一个检验族做 BH，避免把“周期基因富集”与“某阶段比例偏高”混淆。显著结果如下：

{markdown_table(phase_sig) if len(phase_sig) else '本次无 phase 项达到 q<0.05。'}

![Phase distribution](results/biological_validation/dtw/phase_distribution.png)

`peak_times.csv` 列出每个基因在本数据中的全程最大值时间（如有并列取最早），是描述性观察，可能受插补影响；不代表每个周期的真实峰，也不将 t=30 min 等机械指认为 S phase。DTW 聚类允许时间变形，簇内阶段混合是可能且应报告的现象。

## 7. 可选机制解释

[Pramila et al. (2002)](https://pmc.ncbi.nlm.nih.gov/articles/PMC187489/) 提供 Mcm1 与 Yox1/Yhp1 参与 ECB、M/G1 转录调控的实验背景。本分析仅展示相关基因的簇归属与观测峰时间，不能由共簇、相关性或峰值顺序推出调控因果；参考机制不等于本次聚类已重现实验证据。

{markdown_table(examples)}

平台中的某些历史符号与现代名称可能不同，表中主要使用下载平台的符号；YDR451C 的平台 Gene 栏为空，YHP1 名称单独由 [SGD 的 YDR451C 条目](https://wiki.yeastgenome.org/index.php/YDR451C) 确认并标注，未用于改变聚类或参考匹配。`mechanism_examples.csv` 保存完整可追溯例表。

## 8. 模拟集成验证（与真实验证分开）

补充模拟包含早峰、晚峰、振荡、瞬时上调、瞬时下调和近似平坦噪声六类；N=120，T=25；seeds=101,102,103。分别改变 noise、shift、amplitude、missing、stretch 五因素，各三级；其余默认参数固定。共 {len(benchmark)} 次方法运行，包括幅度条件下的未标准化对照，完整配置见 `results/benchmarks/experiment_config.json`。指标为 ARI、NMI 和耗时，报告均值与样本标准差。

模拟生成器的默认 noise=0.1；noise 因素实验显式覆盖为 0/0.15/0.4。峰函数为 `exp(-((u-c)/s)^2)`，早/晚 c=0.2/0.8、s=0.12；瞬时峰 c=0.5、s=0.08，上下调正负相反；振荡为 `sin(4πu)`；平坦模板为 0。位移为均匀随机平移并裁剪端点；幅度乘以 lognormal 系数；stretch 用 `u**exp(U(-stretch,stretch))` 的单调映射；缺失独立随机产生。随机生成真值另存，不传入聚类。

本模拟使用 **phase-sensitive** 真值；早晚峰是同形位移，DTW 可将它们视为相似，可能损害 ARI；这不自动表示 DP 实现出错。未假定 DTW 总优于基线。下表给出相位扰动的例子：

{markdown_table(bench_compare)}

![Supplementary benchmark](results/benchmarks/benchmark_ari.png)

本表是 E 为运行验收补充的有限模拟，未提供 A/D 所需的完整研究性 benchmark、独立调参流程、shape-based 对照和规模实验，必须与组员正式版本合并；真实富集结果不能取代它们。

## 9. 运行验收与复现

环境：Python {platform.python_version()}，Windows；使用新建 `.venv` 从 requirements.txt 安装依赖，精确版本锁在 `requirements-lock.txt`。测试覆盖 DTW 手算/边界、window0=Euclidean、距离对称、常值 Pearson、重复/不等间隔插值、聚类并列处理、Fisher 与超几何一致性和 BH 校正。

验收记录：

```json
{json.dumps(compact_checks, ensure_ascii=False, indent=2)}
```

命令、退出码与逐条日志见 `results/acceptance.json` 和 `results/acceptance_logs/`；完整命令和数据协议见 [README.md](README.md)。正式真实数据可直接执行：

```powershell
cd D:\\BMI\\project11
.\\.venv\\Scripts\\python.exe main.py --input data/gse3635/gene_log10ratio_with_missing.csv --method dtw --k {params['k']} --window {params['window']} --output results/reproduce_dtw --seed 11
.\\.venv\\Scripts\\python.exe biological_validation.py --assignments results/reproduce_dtw/assignments.csv --reference data/references/spellman800.csv --expression results/reproduce_dtw/processed_input.csv --output results/reproduce_dtw
```

## 10. 交付清单和仍需组员协作的事项

- [x] `data/gse3635/`：原始 tar/GPR、官方 normalized matrix/SOFT、platform、samples、ID mapping、QC、非标准化/标准化矩阵、来源 SHA256。
- [x] `biological_validation.py`：精确 ID 匹配、参考覆盖、Fisher/BH、phase summary、observed peak。
- [x] `results/biological_validation/{{dtw,euclidean,pearson}}/`：全体 assignments、medoids、轨迹图、heatmap、富集和配置。
- [x] `main.py`、requirements、README、运行检查和源码哈希。
- [x] 生物学支持与模拟准确率严格区分，来源及统计表可追溯。
- [ ] A–D 真实成员姓名、源码版本、提交号与正式模拟 benchmark 尚未提供，详见 [integration_handover.md](docs/integration_handover.md)。当前相关模块是本次独立集成实现，不能标成组员已提交成果。

本次没有生物学重复，z-score 会丢失幅度，缺失值插补和历史 ID 不匹配均有影响；簇分离度较低、扫描日期混杂、基因间依赖和小候选参数范围限制了解释。所有结果均应作为探索性课程分析。

## 11. 参考文献与用途

1. [GEO GSE3635](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE3635)、[GPL1914](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL1914)：本次实验与平台的原始 metadata。
2. [Spellman et al. (1998), DOI 10.1091/mbc.9.12.3273](https://pmc.ncbi.nlm.nih.gov/articles/PMC25624/)：800 周期基因参考的原始研究，不是全部簇 Ground Truth。
3. [Bioconductor yeastCC v1.50.0](https://bioconductor.org/packages/3.22/data/experiment/html/yeastCC.html)：参考名单和历史阶段表的实际下载来源；包文件 URL 及 SHA256 见 manifest。
4. [Pramila et al. (2002), DOI 10.1101/gad.1034302](https://pmc.ncbi.nlm.nih.gov/articles/PMC187489/)：可选 ECB/Mcm1/Yox1/Yhp1 机制背景，与 GSE3635 关联论文一致。
5. [Aach & Church (2001), DOI 10.1093/bioinformatics/17.6.495](https://pubmed.ncbi.nlm.nih.gov/11395426/)：时间序列对齐方法背景；本项目的 DP 与聚类代码自行编写，未将该论文当作现成完整聚类程序。

来源声明：保留原始文件的使用限制；CellCycle98 明示非商业用途，yeastCC 包为 Artistic-2.0。本交付用于课程研究，未额外授予上游数据再分发权利。
'''
    (root/'Member_E_summary.md').write_text(lines, encoding='utf-8')
    manifest = {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in sorted(root.rglob('*')) if p.is_file() and '.venv' not in str(p)
                and '__pycache__' not in str(p) and p.suffix in ['.py','.ps1']}
    (root/'code_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print('Wrote Member_E_summary.md and code_manifest.json')


if __name__ == '__main__':
    report()
