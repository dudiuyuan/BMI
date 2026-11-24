"""Frozen team expression clustering, followed by conditional Spellman GO enrichment.

Run: python -B run_analysis.py
All new files stay below this script's directory. No network access is used.
"""
import os
os.environ.setdefault('MPLCONFIGDIR', str(__import__('pathlib').Path(__file__).resolve().parent / '.mplconfig'))
import csv
import hashlib
import json
import platform
import shutil
import sys
import time
from collections import Counter, defaultdict
from dataclasses import asdict
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np
import scipy
from scipy.spatial.distance import pdist, squareform
from scipy.stats import hypergeom, fisher_exact
from numba import set_num_threads
from accelerate_team import compile_modules, load, array_distance

HERE = Path(__file__).resolve().parent
TEAM = HERE.parent.parent
METHODS = ['standard','time','time-dynamics','diagonal-time-dynamics','diagonal-time','euclidean','pearson']
K, SEED, N_INIT = 5, 42, 5
SOURCES = ['metrics.py','benchmark.py','summarize_results.py','zrz/dtw/dtw.py','zrz/dtw/euclidean.py',
           'zrz/dtw/pearson.py','xqf/xqf_clustering/k_medoids.py','xqf/xqf_clustering/cluster_analysis.py','xqf/clustering']


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def dump(p, obj):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')


def read(p):
    with p.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


def write(p, fields, rows):
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open('w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def protected_files():
    result = {}
    for directory, dirs, files in os.walk(TEAM):
        dirs[:] = [d for d in dirs if d not in ['.git','.venv','__pycache__'] and Path(directory) / d != HERE.parent]
        for name in files:
            p = Path(directory) / name
            result[str(p.relative_to(TEAM))] = digest(p)
    return result


def bh(p):
    p = np.asarray(p, dtype=float)
    order = np.argsort(p, kind='stable')
    q = np.minimum.accumulate((p[order] * len(p) / np.arange(1,len(p)+1))[::-1])[::-1]
    result = np.empty(len(p)); result[order] = np.minimum(q, 1.)
    return result


def validate_acceleration(original, function, cfg_type, raw_cluster, fast_cluster, X, times):
    rng = np.random.default_rng(20261124)
    examples = [(X[a],X[b],times) for a,b in rng.integers(0,len(X),size=(40,2))]
    for T in [2,5,13,25]:
        t = np.cumsum(rng.uniform(.1,3.,T))
        examples += [(rng.normal(size=T),rng.normal(size=T),t), (np.zeros(T),np.ones(T),t)]
    errors = []
    for method in METHODS[:5]:
        c = original.Config(**original.MODES[method][1]); cfg = cfg_type(*asdict(c).values())
        for x,y,t in examples:
            expected = original.dtw_distance(list(x),list(y),list(t),c)
            observed = function(np.asarray(x),np.asarray(y),np.asarray(t),cfg)
            errors.append(abs(expected-observed))
            assert np.isclose(expected,observed,rtol=1e-12,atol=1e-12), (method,expected,observed)
            span=(t[-1]-t[0]) or 1.
            normalized=np.array([(v-t[0])/span for v in t])
            vx=np.array(original._speed(list(x),list(normalized)))
            vy=np.array(original._speed(list(y),list(normalized)))
            array_result=array_distance(np.asarray(x),np.asarray(y),normalized,vx,vy,cfg)
            errors.append(abs(expected-array_result))
            assert np.isclose(expected,array_result,rtol=1e-12,atol=1e-12)
    cluster_checks = 0
    for n in [6,12,30]:
        D = squareform(pdist(rng.normal(size=(n,5))))
        for seed in [0,11,42]:
            assert raw_cluster.k_medoids(D.tolist(),3,seed=seed) == fast_cluster.k_medoids(D,3,seed=seed)
            cluster_checks += 1
    for D in [np.zeros((8,8)), np.ones((8,8))-np.eye(8)]:
        assert raw_cluster.k_medoids(D.tolist(),3) == fast_cluster.k_medoids(D,3)
        cluster_checks += 1
    return dict(dtw_pairs=len(errors),max_absolute_error=max(errors),cluster_exact_checks=cluster_checks,
                arithmetic='Original-source JIT and equivalent array recurrence both checked; fastmath disabled; Python random.Random and main k_medoids unchanged')


def clustering():
    files = list((TEAM/'wpj').glob('predata/*/outputs/real/GSE3635/run01/processed/centered.csv'))
    files = [p for p in files if '__MACOSX' not in str(p)]
    assert len(files) == 1
    source = files[0]
    original, function, matrix, cfg_type, raw_cluster, fast_cluster = compile_modules(TEAM)
    ids,times,values = original.read_matrix(source)
    X,times = np.asarray(values),np.asarray(times)
    assert len(ids) == 6214 and X.shape[1] == 13
    # Do not read GO or phase labels until every expression clustering is frozen.
    validation = validate_acceleration(original,function,cfg_type,raw_cluster,fast_cluster,X,times)
    dump(HERE/'acceleration_validation.json',validation)
    set_num_threads(min(8,os.cpu_count() or 1))
    sys.modules['dtw'] = original
    euclidean = load('go_original_euclidean',TEAM/'zrz/dtw/euclidean.py')
    pearson = load('go_original_pearson',TEAM/'zrz/dtw/pearson.py')
    original_hashes = {name:digest(TEAM/name) for name in SOURCES}
    for name in SOURCES:
        dest = HERE/'source_snapshot'/name; dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(TEAM/name,dest)
    parameters = {m:asdict(original.Config(**original.MODES[m][1])) if m in original.MODES else {} for m in METHODS}
    configuration = dict(input=str(source),input_sha256=digest(source),genes=len(ids),timepoints=times.tolist(),
        preprocessing='wpj centered.csv unchanged; no additional scaling; differs from old E zscore/QC pipeline',
        k=K,seed=SEED,n_init=N_INIT,max_iter=100,parameters=parameters,source_sha256=original_hashes,
        choice='K=5 fixed before GO, matches xqf reference setting; not inferred true phase count; no GO parameter selection',
        python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)
    config_path = HERE/'configuration.json'
    if config_path.exists():
        assert json.loads(config_path.read_text(encoding='utf-8')) == configuration, 'Settings changed; use a new directory'
    dump(config_path,configuration)
    write(HERE/'expression_input.csv',['gene_id']+[f't{int(t)}' for t in times],
          [dict(zip(['gene_id']+[f't{int(t)}' for t in times],[gene,*row])) for gene,row in zip(ids,X)])
    assignments, runs = {}, []
    for method in METHODS:
        folder = HERE/'clustering'/method; folder.mkdir(parents=True,exist_ok=True)
        predictions = folder/'assignments.csv'
        if predictions.exists() and (folder/'run.json').exists():
            run_info=json.loads((folder/'run.json').read_text())
            assert digest(predictions)==run_info['assignment_sha256'], 'Saved assignment integrity failure'
            assignments[method] = read(predictions)
            assert [r['gene_id'] for r in assignments[method]]==ids
            assert set(int(r['cluster_id']) for r in assignments[method])==set(range(K))
            runs.append(run_info)
            print(f'{method}: using frozen saved clustering',flush=True)
            continue
        started = time.perf_counter()
        distance_file = folder/'distance_matrix.npy'
        if distance_file.exists():
            D = np.load(distance_file)
            distance_seconds = None
        else:
            if method in original.MODES:
                cfg = cfg_type(*parameters[method].values())
                D = matrix(X,times,cfg)
            elif method == 'euclidean':
                D = squareform(pdist(X,'euclidean'))
            else:
                centered = X-X.mean(axis=1,keepdims=True)
                norm = np.linalg.norm(centered,axis=1)
                flat = np.sqrt((centered**2).mean(axis=1)) <= 1e-12*np.maximum(1.,np.abs(X).max(axis=1))
                normalized = centered/np.where(norm>0,norm,1)[:,None]
                D = np.clip(1-normalized@normalized.T,0.,2.)
                D[flat,:]=1.; D[:,flat]=1.; D[np.ix_(flat,flat)]=0.
                np.fill_diagonal(D,0.)
            distance_seconds = time.perf_counter()-started
            np.save(distance_file,D)
        assert np.isfinite(D).all() and np.all(D>=0) and np.array_equal(D,D.T) and np.all(np.diag(D)==0)
        # Baselines are numerically verified against supplied original scalar code.
        if method in ['euclidean','pearson']:
            f = euclidean.euclidean_distance if method=='euclidean' else pearson.pearson_distance
            for a,b in np.random.default_rng(42).integers(0,len(X),size=(100,2)):
                assert np.isclose(D[a,b],f(list(X[a]),list(X[b]),times),rtol=1e-12,atol=1e-12)
        print(f'{method}: distance matrix ready ({distance_seconds}s); clustering',flush=True)
        start = time.perf_counter()
        labels, medoids = fast_cluster.k_medoids(D,K,seed=SEED,n_init=N_INIT,max_iter=100)
        cluster_seconds = time.perf_counter()-start
        assert len(labels)==len(ids) and set(labels)==set(range(K)) and len(set(medoids))==K
        rows=[dict(gene_id=g,cluster_id=c,is_medoid=int(i in medoids)) for i,(g,c) in enumerate(zip(ids,labels))]
        write(predictions,['gene_id','cluster_id','is_medoid'],rows)
        write(folder/'medoids.csv',['cluster_id','gene_id']+[f't{int(t)}' for t in times],
            [dict(zip(['cluster_id','gene_id']+[f't{int(t)}' for t in times],[c,ids[i],*X[i]])) for c,i in enumerate(medoids)])
        run=dict(method=method,distance_seconds=distance_seconds,clustering_seconds=cluster_seconds,
                 objective=float(sum(D[i,medoids[c]] for i,c in enumerate(labels))),
                 cluster_sizes=dict(Counter(labels)),assignment_sha256=digest(predictions),
                 distance_sha256=digest(distance_file),medoid_gene_ids=[ids[i] for i in medoids])
        dump(folder/'run.json',run)
        assignments[method]=rows; runs.append(run)
        print(f'{method}: clusters {dict(Counter(labels))}',flush=True)
        del D
    dump(HERE/'clustering_runs.json',runs)
    return ids,assignments


def enrich(ids,assignments,exclude_shared=False):
    go_root=HERE.parent/'go_annotation/results'
    manifest=json.loads((go_root/'manifest.json').read_text(encoding='utf-8'))
    assert digest(go_root/'GO_long.csv')==manifest['output_sha256']['GO_long.csv']
    go=read(go_root/'GO_long.csv'); reference=read(go_root/'spellman800_input.csv')
    refs={r['gene_id'] for r in reference}
    terms=defaultdict(set); names={}; aspects={}
    shared={r['gene_id'] for r in read(go_root/'gene_audit.csv') if r.get('shared_entry')=='True'}
    if not shared:
        shared={g for group in manifest['stats'].get('shared_entries',{}).values() for g in group}
    for row in go:
        g,code=row['gene_id'],row['GO_code']
        if g in ids and (not exclude_shared or g not in shared):
            terms[code].add(g); names[code]=row['GO_meaning']; aspects[code]=row['GO_aspect']
    annotated=set().union(*terms.values())
    N=len(annotated)
    tested=sorted(t for t,genes in terms.items() if len(genes)>=3 and t not in ['GO:0008150','GO:0003674','GO:0005575'])
    rows,coverage=[],[]
    for method,pred in assignments.items():
        for cluster in range(K):
            full={r['gene_id'] for r in pred if int(r['cluster_id'])==cluster}
            group=full & annotated; n=len(group)
            coverage.append(dict(method=method,cluster_id=cluster,total_cluster_genes=len(full),
                spellman_genes=len(full&refs),GO_annotated_genes=n,background_N=N,annotation_fraction=n/len(full)))
            for code in tested:
                term_genes=terms[code]; hits=sorted(group&term_genes); a=len(hits); total=len(term_genes)
                b,c,d=n-a,total-a,N-n-total+a
                assert min(a,b,c,d)>=0
                rows.append(dict(method=method,cluster_id=cluster,GO_code=code,GO_meaning=names[code],
                    GO_aspect=aspects[code],a=a,b=b,c=c,d=d,cluster_n=n,background_N=N,background_K=total,
                    fold_enrichment=(a/n)/(total/N) if n else 0.,hit_genes=';'.join(hits)))
    p=hypergeom.sf(np.array([r['a'] for r in rows])-1,N,
                   np.array([r['background_K'] for r in rows]),np.array([r['cluster_n'] for r in rows]))
    assert np.isfinite(p).all()
    global_q=bh(p)
    for row,pv,q in zip(rows,p,global_q): row.update(p_value=float(pv),q_global=float(q))
    for method in METHODS:
        subset=[r for r in rows if r['method']==method]
        for row,q in zip(subset,bh([r['p_value'] for r in subset])): row['q_within_method']=float(q)
    for row in sorted(rows,key=lambda r:r['p_value'])[:20]+rows[::max(1,len(rows)//20)]:
        expected=fisher_exact([[row['a'],row['b']],[row['c'],row['d']]],alternative='greater').pvalue
        assert np.isclose(expected,row['p_value'],rtol=1e-10,atol=1e-300)
    assert np.allclose(bh([.01,.04,.03,1.]),[.04,.05333333333333334,.05333333333333334,1.])
    folder=HERE/('enrichment_no_shared' if exclude_shared else 'enrichment')
    rows.sort(key=lambda r:(r['method'],r['cluster_id'],r['q_global'],r['p_value'],r['GO_code']))
    fields=['method','cluster_id','GO_code','GO_meaning','GO_aspect','a','b','c','d','cluster_n','background_N',
            'background_K','fold_enrichment','p_value','q_global','q_within_method','hit_genes']
    write(folder/'all_tests.csv',fields,rows)
    significant=[r for r in rows if r['q_global']<.05 and r['fold_enrichment']>1]
    write(folder/'significant.csv',fields,significant)
    write(folder/'cluster_coverage.csv',list(coverage[0]),coverage)
    for method in METHODS:
        for cluster in range(K):
            write(folder/method/f'cluster_{cluster}.csv',fields,[r for r in rows if r['method']==method and r['cluster_id']==cluster])
    for r in reference:
        g=r['gene_id']
        r['in_expression']=g in ids; r['has_usable_GO']=g in annotated
        r['excluded_shared_entry']=exclude_shared and g in shared
    write(folder/'spellman_800_coverage.csv',list(reference[0]),reference)
    summary=dict(background_N=N,expression_genes=len(ids),spellman_in_expression=len(set(ids)&refs),
        total_terms=len(terms),tested_terms=len(tested),tests=len(rows),significant_tests=len(significant),
        release=manifest['uniprot_release'],minimum_background_term_genes=3,
        global_BH_family='all seven methods x five clusters x all tested GO terms (including zero hits)',
        universe='actual clustered genes with GO in the downloaded Spellman-only annotation table',
        exclude_shared=exclude_shared,excluded_shared_genes=sorted(shared&set(ids)) if exclude_shared else [],
        significant_by_method=dict(Counter(r['method'] for r in significant)))
    dump(folder/'summary.json',summary)
    dump(folder/'input_manifest.json',dict(GO_long_sha256=digest(go_root/'GO_long.csv'),
        GO_manifest_sha256=digest(go_root/'manifest.json'),GO_source=str(go_root),
        assignments={m:digest(HERE/'clustering'/m/'assignments.csv') for m in METHODS}))
    return rows,coverage,summary


def plots(assignments,rows,coverage):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for method in METHODS:
        medoids=read(HERE/'clustering'/method/'medoids.csv')
        columns=[c for c in medoids[0] if c.startswith('t')]
        fig,ax=plt.subplots(figsize=(8,4))
        for r in medoids: ax.plot([float(c[1:]) for c in columns],[float(r[c]) for c in columns],marker='o',label=f"Cluster {r['cluster_id']}: {r['gene_id']}")
        ax.set(xlabel='Time (min)',ylabel='Centered expression',title=f'{method}: medoid expression trajectories')
        ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(HERE/'clustering'/method/'medoid_profiles.png',dpi=150); plt.close(fig)
    top=[]
    for method in METHODS:
        for cluster in range(K):
            hits=[r for r in rows if r['method']==method and r['cluster_id']==cluster and r['q_global']<.05 and r['GO_aspect']=='biological_process']
            if hits: top.append(min(hits,key=lambda r:(r['q_global'],r['p_value'])))
    if top:
        fig,ax=plt.subplots(figsize=(12,max(5,len(top)*.32)))
        labels=[f"{r['method']} C{r['cluster_id']}: {r['GO_meaning']}" for r in top]
        ax.barh(np.arange(len(top)),-np.log10([max(r['q_global'],1e-300) for r in top]))
        ax.set_yticks(np.arange(len(top)),labels,fontsize=8); ax.invert_yaxis()
        ax.set(xlabel='-log10(global BH q)')
        fig.suptitle('Top biological process per expression cluster\nConditional enrichment in GO-annotated Spellman genes',fontsize=12)
        fig.tight_layout(rect=(0,0,1,.92));fig.savefig(HERE/'top_GO_BP.png',dpi=180);plt.close(fig)


def report(rows,coverage,summary,sensitivity):
    lines=['# GSE3635 表达轨迹聚类与 GO 富集分析','',
        '分析日期：2026-11-24。所有新增内容均位于本报告所在目录；其他组员文件仅只读。','',
        '## 1. 代码核查结论','',
        '**是按表达轨迹聚类。** zrz/dtw/dtw.py 读取 gene_id 和全部时间点表达值，逐基因对计算端点到端点的全局 DTW；局部代价为表达平方差，部分模式再加时间差、变化速率和拉伸惩罚。没有取峰值时间作为聚类特征，没有读取 GO 或 Spellman phase。DTW 本身只输出距离矩阵，不执行聚类。',
        'xqf/xqf_clustering/k_medoids.py 对距离矩阵执行多起点交替 k-medoids：分配到最近 medoid，再选簇内总距离最小的基因。输入若为 DTW 表达距离，则输出是轨迹相似性分组。cluster_analysis.py 的 compare_k 用 silhouette 比较 K，其他函数绘图。',
        '**注意文件差异：** xqf/clustering（无 .py 后缀）的 assign_clusters 中 labels[gene] = closest 缩进在外层循环外，可能遗留 -1 标签；本次使用 benchmark.py 指向的 xqf/xqf_clustering/k_medoids.py，该 .py 文件缩进正确，不使用有问题的无后缀副本。','',
        '## 2. 三个 Python 文件分别做什么','',
        '| 文件 | 实际工作 | 与 GO 的关系 |','|---|---|---|',
        '| metrics.py | 不依赖第三方库实现 Hubert–Arabie ARI 与算术平均熵 NMI，按列联计数比较已有真标签和预测标签；不是聚类算法。 | 不读取 GO。 |',
        '| benchmark.py | 冻结并校验 zrz 三种距离源码和 xqf k-medoids；读取 wpj 模拟 centered.csv；运行五种 DTW、Euclidean、Pearson；phase6 使用 true_cluster/K=6，shape5 使用 shape_cluster/K=5；缓存距离和预测；聚类后连接真值计算 ARI/NMI，保存预测、列联表、环境、来源校验和。 | 只处理 smoke/baseline/challenge/sweep 模拟组，不自动运行真实数据或 GO 富集。 |',
        '| summarize_results.py | 汇总 results.csv 中每条件五个独立数据 seed 3100–3104 的 ARI/NMI 和耗时均值、样本 SD；相对 standard 计算配对差值与胜/平/负；输出 280 个条件汇总、配对表、baseline 表、六组 PNG/SVG 图和结果解读。 | 不重新聚类，不处理 GO。SD 不是置信区间。 |','',
        '完整 sweep 为 100 份数据 × 7 方法 × 2 任务 = 1400 行；默认再加 baseline 的 14 行。benchmark 默认 --team-root 指向旧目录，小组在当前目录重跑时应显式指定 --team-root D:\\BMI\\project11\\lzy。该说明不代表本次运行了模拟 benchmark。',
        'metrics.py 另进行了 160 次与 scikit-learn 的 ARI/NMI 对照，最大误差约 4.65e-16，见 metrics_validation.json。',
        '距离模式曾按模拟真标签 ARI 调参，所以距离超参数来源是模拟监督选择；本次真实数据的分簇和 GO 分析没有读取真实 phase 标签来调参。','',
        '## 3. 本次真实数据分簇','',
        '未发现这套 zrz/xqf 代码生成的真实样本分簇文件，因此本次新运行；没有把旧 E 实现的 assignments 冒充组员代码结果。',
        '输入为 wpj/predata/数据准备/outputs/real/GSE3635/run01/processed/centered.csv，全部 6214 基因 × 13 时间点；直接使用已居中表达值，保留幅度差异，不再 z-score，不仅限于 800 基因。该预处理与旧 E 的 6136 基因/质量屏蔽/z-score 分析不同，结果不可直接混同。',
        '统一 K=5、seed=42、n_init=5、max_iter=100。K=5 是在 GO 分析前固定的比较设置（对应 xqf 绘图参考值），不是从本次数据估计的最优 K，更不是先验断言五个真实阶段。未用 GO 分数或 phase 选择模式、K 或种子。',
        '为支持 6214 基因的全矩阵，用 Numba 执行与原始三状态 DTW 相同的数组递推（以数组复用替代逐对反复创建列表，并预计算相同时间坐标及局部速率），逐对对照原始函数验证；同时编译 k-medoids 的两个计算辅助函数。主 k_medoids 和 Python random.Random 保持原逻辑，源文件没有修改。原函数 JIT 另以等价二分替代不支持的 bisect 作对照。Euclidean/Pearson 使用向量计算，并与原始函数核对。精度核查见 acceleration_validation.json；此加速耗时不能冒充原纯 Python benchmark 耗时。','',
        '| 方法 | C0 | C1 | C2 | C3 | C4 |','|---|---:|---:|---:|---:|---:|']
    for method in METHODS:
        sizes=[next(r['total_cluster_genes'] for r in coverage if r['method']==method and r['cluster_id']==c) for c in range(K)]
        lines.append('| '+method+' | '+' | '.join(map(str,sizes))+' |')
    lines += ['', '## 4. GO 数据范围与统计方法','',
        f"使用此前下载的 UniProt {summary['release']} GO_long.csv，原始注释对象仅为 Spellman 800 名单。6214 个聚类基因中有 {summary['spellman_in_expression']} 个属于该名单，其中 {summary['background_N']} 个有可用 GO。其他基因不视作 GO 阴性。",
        '**本次回答的是：在有 GO 注释、且进入表达聚类的 Spellman 基因中，每个簇相对该参考子集富集哪些功能。不能据此声称已对全部 6214 个基因完成全背景功能验证。** 如需全体聚类基因背景，必须先为全部输入基因按统一规则获取 GO，再另立分析。',
        '按 gene_id 与 GO_code 去重，一个基因的同一 GO 只计一次。保留 BP/MF/CC 及所有原始证据类型，包括 IEA；不将它们全称为实验证据。只使用 UniProt 返回的直接注释，不扩展祖先；根术语排除。基因是统计单位，不按长表行数计数。',
        f"背景 N={summary['background_N']}；预先仅检验背景命中至少 3 基因的 {summary['tested_terms']} 个术语。所有方法与簇使用同一术语集合，零命中也保留为检验；总计 {summary['tests']} 项。",
        '对每个簇与 GO，a=簇内有该 GO，b=簇内无该 GO，c=簇外有该 GO，d=簇外无该 GO，均限于上述可注释背景。单侧 Fisher exact（greater）检验富集，向量化计算用等价超几何右尾，并抽样对照 scipy fisher_exact。',
        '主结论使用跨 7 方法 × 5 簇 × 所有术语的全局 BH 校正 q_global<0.05，且 fold_enrichment>1。另提供各方法内 BH q_within_method，不能与主阈值混用。GO 术语存在层级与统计依赖，结果作探索性生物学支持，不等于聚类准确率；不因显著数更多就认定某方法最优。',
        'GO 官方建议自定义背景为实际可能进入所测集合的基因：https://geneontology.org/docs/go-enrichment-analysis/','',
        '## 5. 逐簇结果（主分析）','',
        '下表每簇报告最显著术语；“簇大小”为全表达簇，“GO基因”为实际用于富集的参考子集。完整 BP/MF/CC 结果、四格计数与命中基因见 enrichment/<method>/cluster_*.csv。不同方法相同编号不是相同生物学簇。','',
        '![各簇最显著生物过程](top_GO_BP.png)','',
        '| 方法 | 簇 | 簇大小 | GO基因 | 显著GO数 | 最显著 GO | 命中/GO基因 | Fold | 全局q |',
        '|---|---:|---:|---:|---:|---|---:|---:|---:|']
    for cov in coverage:
        hits=[r for r in rows if r['method']==cov['method'] and r['cluster_id']==cov['cluster_id'] and r['q_global']<.05 and r['fold_enrichment']>1]
        best=min(hits,key=lambda r:(r['q_global'],r['p_value'])) if hits else None
        cells=[cov['method'],str(cov['cluster_id']),str(cov['total_cluster_genes']),str(cov['GO_annotated_genes']),str(len(hits))]
        cells += ([best['GO_code']+' '+best['GO_meaning'],f"{best['a']}/{best['cluster_n']}",f"{best['fold_enrichment']:.2f}",f"{best['q_global']:.3g}"] if best else ['无全局显著富集','—','—','—'])
        lines.append('| '+' | '.join(cells)+' |')
    lines += ['', '### 主要生物过程（BP）','']
    for method in METHODS:
        lines.append(f'**{method}**')
        lines.append('')
        lines.append(f'![{method} 实际 medoid 表达轨迹](clustering/{method}/medoid_profiles.png)')
        lines.append('')
        for c in range(K):
            hits=[r for r in rows if r['method']==method and r['cluster_id']==c and r['GO_aspect']=='biological_process' and r['q_global']<.05]
            hits.sort(key=lambda r:(r['q_global'],r['p_value']))
            text='；'.join(f"{r['GO_meaning']} ({r['GO_code']}, {r['a']}/{r['cluster_n']}, fold={r['fold_enrichment']:.2f}, q={r['q_global']:.3g})" for r in hits[:3])
            lines.append(f'- C{c}：'+(text or '没有达到全局 q<0.05 的 BP 术语。'))
        lines.append('')
    lines += ['### 生物学解读','',
        '这些固定参数结果在部分表达簇中提供了 DNA 修复、端粒重组维护、染色质和细胞分裂相关功能支持；未达到阈值的簇不能解释为没有相关功能，部分簇可用 GO 基因数仅为个位数或十余个，统计功效有限。',
        'time-dynamics 与 diagonal-time-dynamics 的 C2 均出现 MCM complex 与 DNA 复制前复合物组装相关注释，说明这些轨迹分组中有已知复制相关基因集中。但 GO 注释支持的是已知功能的一致性，不能从共簇推断调控方向、因果关系或恢复了真实细胞周期阶段。',
        '主分析的部分染色质/核小体注释涉及共享蛋白条目的组蛋白基因；应结合下面的排除共享条目敏感性分析解释，不将重复证据误当作独立验证。具体显著数以保存的完整表为准，不用显著项数量给方法排名。','',
        '## 6. 共享 UniProt 条目的敏感性分析','',
        '先前 UniProt 审计显示 11 个历史基因共享 5 个蛋白条目，可能产生非独立注释。主分析保留原有精确基因映射，补充分析排除进入表达数据的全部共享基因；两者使用相同表达簇，不重新聚类。',
        f"排除名单：{', '.join(sensitivity['excluded_shared_genes'])}。补充背景 N={sensitivity['background_N']}；相同全局 BH 原则下 {sensitivity['significant_tests']} 个显著组合。补充分析重新构建可检验术语集，因此 q 值不宜逐数值直接比较。完整结果在 enrichment_no_shared/。",'',
        '| 方法 | 主分析显著组合 | 排除共享后的显著组合 |','|---|---:|---:|']
    for method in METHODS:
        lines.append(f"| {method} | {summary['significant_by_method'].get(method,0)} | {sensitivity['significant_by_method'].get(method,0)} |")
    lines += ['', '## 7. 复现与文件','',
        '依赖 Python、NumPy、SciPy、Numba、Matplotlib。执行 `python -B run_analysis.py`，已有真实分簇和距离缓存会复用，配置或源码变化会拒绝静默混用。全部结果、配置、源代码快照、原始输入快照都在此目录。',
        'configuration.json 记录完整参数和组员源码 SHA256；enrichment/input_manifest.json 记录 GO 和 assignments 校验和；verification.json 包含精度、Fisher 等价、覆盖与其他组员文件未修改的检查。',
        '主要文件：clustering/<method>/assignments.csv、medoids.csv、medoid_profiles.png、distance_matrix.npy；enrichment/all_tests.csv、significant.csv、cluster_coverage.csv、各簇 CSV；enrichment_no_shared/；source_snapshot/。',
        '本次未运行根目录的模拟 benchmark.py/summarize_results.py，未修改这些文件的默认路径。','',
        '## 8. 解释边界与参考','',
        '结果为固定 K 与已选距离模式的一次真实数据探索；缺少生物学重复，GO 非完整真值，功能相似也不证明调控因果。同一参考集合的基因预选择、共表达依赖、注释覆盖和共享蛋白条目都影响解释。输入仅居中，幅度仍影响 DTW/Euclidean，Pearson 对幅度不敏感。',
        '- Spellman PT et al. (1998), Molecular Biology of the Cell 9:3273–3297. https://doi.org/10.1091/mbc.9.12.3273',
        '- GSE3635 真实数据：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE3635',
        '- UniProt 注释及原始 API 快照：../go_annotation/results/manifest.json。此次没有新联网下载 GO。',
        '- GO 背景与富集解释：https://geneontology.org/docs/go-enrichment-analysis/','']
    (HERE/'分析报告.md').write_text('\n'.join(lines),encoding='utf-8')


def main():
    before_path=HERE/'protected_sources_before.json'
    before=json.loads(before_path.read_text(encoding='utf-8')) if before_path.exists() else protected_files()
    dump(before_path,before)
    ids,assignments=clustering()
    rows,coverage,summary=enrich(ids,assignments)
    _,_,sensitivity=enrich(ids,assignments,True)
    plots(assignments,rows,coverage)
    report(rows,coverage,summary,sensitivity)
    after=protected_files()
    assert before==after, 'Files outside the authorized destination changed; inspect snapshot'
    dump(HERE/'verification.json',dict(status='passed',protected_files=len(before),other_members_unchanged=True,
        acceleration=json.loads((HERE/'acceleration_validation.json').read_text()),
        fisher_vs_hypergeometric_checked=True,BH_known_example_checked=True,main=summary,sensitivity=sensitivity,
        dependencies=dict(numpy=np.__version__,scipy=scipy.__version__),
        new_files={str(p.relative_to(HERE)):digest(p) for p in HERE.rglob('*') if p.is_file() and p.name not in ['verification.json'] and '.mplconfig' not in str(p)}))
    print(json.dumps(summary,ensure_ascii=False,indent=2),flush=True)
    print('Report: '+str(HERE/'分析报告.md'),flush=True)


if __name__=='__main__':
    main()
