"""Reproducible functional/phase validation; no writes outside this directory."""
import os
from pathlib import Path
HERE=Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR',str(HERE/'.mplconfig'))
import csv,hashlib,json,re,sys
from collections import Counter,defaultdict
sys.dont_write_bytecode=True
import numpy as np
from scipy.stats import fisher_exact,hypergeom,pearsonr
from fast_cluster import pair_score,null_scores
from run_seeds import METHODS,SEEDS,OLD,load_csv

TEAM=HERE.parent.parent
CATEGORIES=['DNA replication and repair','Chromatin and silencing','Mitosis and segregation',
 'Cell cycle regulation','Cytoskeleton','Cell wall and cytokinesis','Budding and polarity',
 'Mating and pheromone','Meiosis and sporulation','Transport','Secretion and trafficking',
 'Transcription and RNA processing','Protein processing and turnover','Metabolism']
N_PERM=10000


def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda:f.read(4*1024*1024),b''):h.update(chunk)
    return h.hexdigest()


def dump(p,obj):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')


def write(p,rows,fields=None):
    p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields or list(rows[0]),extrasaction='ignore');w.writeheader();w.writerows(rows)


def protected():
    files={}
    for directory,dirs,names in os.walk(TEAM):
        dirs[:]=[d for d in dirs if d not in ['.git','.venv','__pycache__'] and Path(directory)/d!=HERE]
        for name in names:
            p=Path(directory)/name;files[str(p.relative_to(TEAM))]=sha(p)
    return files


def bh(p):
    p=np.asarray(p);order=np.argsort(p,kind='stable')
    v=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]
    q=np.empty(len(p));q[order]=np.minimum(v,1.);return q


def clean(value):
    # Retain meaningful qualifiers, remove only '(putative)', preserve raw text.
    value=re.sub(r'\(putative\)','',value.lower()).replace('metbolism','metabolism').replace('acod','acid')
    return re.sub(r'\s+',' ',value).strip()


def categories(term):
    rules=[r'dna (replication|repair)|recombination',r'chromatin|silencing|telomere length',
      r'mitosis|chromosome (segregation|condensation|transmission)|sister chromatid',
      r'cell cycle|checkpoint',r'cytoskeleton',r'cell wall|cytokinesis',
      r'bud |polarity|polarized growth|pseudohyphal',r'mating|pheromone|karyogamy|nuclear fusion|homothallic',
      r'meiosis|sporulation|spore formation',r'transport|permease|uptake|export|import|protein targeting',
      r'secretion|membrane trafficking|golgi|endocytosis|vacuole biogenesis',
      r'transcription|rna processing|rna splicing|rna 3.*processing',
      r'protein (glycosylation|degradation|processing|folding|synthesis)',
      r'metabolism|biosynthesis|glycolysis|gluconeogenesis|respiration|oxidative phosphorylation|atp synthesis|tca cycle|diauxic shift|assimilation|one-carbon|ethanol utilization|sucrose utilization|glucose repression']
    result=[category for category,pattern in zip(CATEGORIES,rules) if re.search(pattern,term)]
    # Cell-wall biosynthesis is not classified as generic metabolic biosynthesis.
    if 'Cell wall and cytokinesis' in result and 'Metabolism' in result:result.remove('Metabolism')
    return result


def annotations(ids):
    refroot=HERE.parent/'data/references'
    known=set((refroot/'orf800.txt').read_text().split())
    raw=list(csv.reader((refroot/'CellCycle98.txt').read_text(encoding='latin1').splitlines(),delimiter='\t'))
    ref=[];fine={};coarse={};phase={};mapping={}
    for r in raw:
        if len(r)<9:continue
        gene=r[6].strip().upper();rule='ORF'
        if not gene and r[7].strip()==r[8].strip() and r[7].strip() in known:
            gene=r[7].strip();rule='blank ORF restored from agreeing YPD/SGD names and orf800'
        if gene not in known:continue
        terms={clean(t) for t in r[0].split(';')}-{'' ,'unknown'}
        groups=set()
        for term in terms:
            mapped=categories(term);groups.update(mapped);mapping[term]=mapped
        fine[gene]=terms;coarse[gene]=groups;phase[gene]={r[3].strip()} if r[3].strip() else set()
        ref.append(dict(gene_id=gene,Process_raw=r[0],Function_raw=r[1],phase=r[3],
             terms=';'.join(sorted(terms)),categories=';'.join(sorted(groups)),
             putative='putative' in r[0].lower(),mapping_rule=rule,in_expression=gene in ids,
             coarse_unmapped=bool(terms) and not groups))
    assert len(ref)==len(known)==800 and {r['gene_id'] for r in ref}==known
    write(HERE/'input/spellman_annotations.csv',ref)
    maprows=[dict(term=t,categories=';'.join(mapping[t]),unmapped_other=not bool(mapping[t])) for t in sorted(mapping)]
    target=HERE/'input/process_mapping.csv'
    if target.exists():
        assert load_csv(target)==[{k:str(v) for k,v in r.items()} for r in maprows], 'Mapping changed; use a new output directory'
    write(target,maprows)
    gofile=HERE.parent/'go_annotation/results/GO_long.csv'
    gm=json.loads((gofile.parent/'manifest.json').read_text(encoding='utf-8'))
    assert sha(gofile)==gm['output_sha256']['GO_long.csv']
    go=defaultdict(set);bp=defaultdict(set);termnames={};aspects={}
    for r in load_csv(gofile):
        if r['GO_code'] in ['GO:0008150','GO:0003674','GO:0005575']:continue
        go[r['gene_id']].add(r['GO_code']);termnames[r['GO_code']]=r['GO_meaning'];aspects[r['GO_code']]=r['GO_aspect']
        if r['GO_aspect']=='biological_process':bp[r['gene_id']].add(r['GO_code'])
    bases={'Process_coarse':coarse,'Process_fine':fine,'Phase':phase,'GO_all':go}
    universes={name:sorted(g for g,tags in values.items() if tags and g in ids) for name,values in bases.items()}
    stats=dict(original_genes=800,process_known=sum(bool(v) for v in fine.values()),
        process_unknown=sum(r['Process_raw'].strip().lower()=='unknown' for r in ref),
        process_blank=sum(not r['Process_raw'].strip() for r in ref),
        fine_terms=len(mapping),coarse_categories=len(CATEGORIES),
        expression_genes=len(ids),spellman_in_expression=len(set(ids)&known),
        universes={k:len(v) for k,v in universes.items()},GO_release=gm['uniprot_release'],
        restored_ORFs=[r['gene_id'] for r in ref if r['mapping_rule']!='ORF'],
        unmapped_coarse_terms=[r['term'] for r in maprows if r['unmapped_other']],
        Process_source_sha256=sha(refroot/'CellCycle98.txt'),GO_long_sha256=sha(gofile),
        mapping_sha256=sha(target))
    dump(HERE/'annotation_summary.json',stats)
    return bases,universes,bp,termnames,aspects,ref,stats


def enrichment(ids,bases,universes,termnames,aspects):
    lookup={g:i for i,g in enumerate(ids)};tables={};allrows=[]
    rawfunctions={r['gene_id']:r['Function_raw'] for r in load_csv(HERE/'input/spellman_annotations.csv')}
    for basis,ann in bases.items():
        universe=set(universes[basis]);N=len(universe)
        terms=defaultdict(set)
        for g in universe:
            for t in ann[g]:terms[t].add(g)
        tested=sorted(t for t,genes in terms.items() if len(genes)>=5 and len(genes)<N)
        records=[]
        for method in METHODS:
            labels=np.load(HERE/'seed_runs'/method/'seed_42.npz')['labels']
            for cluster in range(5):
                group={g for g in universe if labels[lookup[g]]==cluster};n=len(group)
                for t in tested:
                    hits=sorted(group&terms[t]);a=len(hits);K=len(terms[t]);b=n-a;c=K-a;d=N-n-K+a
                    denominator=np.sqrt(float(n*(N-n)*K*(N-K)))
                    phi=(a*d-b*c)/denominator if denominator else None
                    two=fisher_exact([[a,b],[c,d]],alternative='two-sided').pvalue
                    greater=hypergeom.sf(a-1,N,K,n)
                    records.append(dict(basis=basis,method=method,cluster=cluster,term=t,
                        term_name=termnames.get(t,t),GO_aspect=aspects.get(t,''),N=N,K=K,n=n,a=a,b=b,c=c,d=d,
                        phi=phi,fold=(a/n)/(K/N) if n else None,p_phi=float(two),p_enrichment=float(greater),
                        hit_genes=';'.join(hits),Functions=';'.join(g+':'+rawfunctions.get(g,'') for g in hits)))
        for r,q1,q2 in zip(records,bh([r['p_phi'] for r in records]),bh([r['p_enrichment'] for r in records])):
            r.update(q_phi=float(q1),q_enrichment=float(q2),significant_enrichment=bool(q2<.05 and r['a']>=3 and r['fold'] is not None and r['fold']>1))
        for r in records[::max(1,len(records)//10)]:
            assert np.isclose(r['p_enrichment'],fisher_exact([[r['a'],r['b']],[r['c'],r['d']]],alternative='greater').pvalue)
            if r['phi'] is not None:
                x=np.array([int(g in terms[r['term']]) for g in sorted(universe)])
                lab=np.load(HERE/'seed_runs'/r['method']/'seed_42.npz')['labels']
                y=np.array([int(lab[lookup[g]]==r['cluster']) for g in sorted(universe)])
                assert np.isclose(r['phi'],pearsonr(x,y).statistic,atol=1e-12)
        write(HERE/'tables'/f'{basis}_all_tests.csv',records)
        write(HERE/'tables'/f'{basis}_significant.csv',[r for r in records if r['significant_enrichment']],list(records[0]))
        cluster_summary=[]
        for method in METHODS:
            for cluster in range(5):
                subset=[r for r in records if r['method']==method and r['cluster']==cluster]
                write(HERE/'tables/per_cluster'/basis/method/f'cluster_{cluster}.csv',subset)
                significant=sorted([r for r in subset if r['significant_enrichment']],key=lambda r:r['q_enrichment'])
                best=significant[0] if significant else None
                cluster_summary.append(dict(method=method,cluster=cluster,background_N=N,annotated_cluster_n=subset[0]['n'],
                    tested_labels=len(subset),significant_phi=sum(r['q_phi']<.05 for r in subset),
                    significant_enrichment=len(significant),top_term=best['term_name'] if best else '',
                    top_term_q=best['q_enrichment'] if best else '',top_hit_genes=best['hit_genes'] if best else ''))
        write(HERE/'tables'/f'{basis}_cluster_summary.csv',cluster_summary)
        tables[basis]=records;allrows+=records
    # Known-cycle-gene enrichment has a different, full-expression background.
    refs=set(bases['Phase']);universe=set(ids);N=len(ids);K=len(refs&universe);membership=[]
    for method in METHODS:
        labels=np.load(HERE/'seed_runs'/method/'seed_42.npz')['labels']
        for c in range(5):
            group={g for i,g in enumerate(ids) if labels[i]==c};a=len(group&refs);n=len(group)
            membership.append(dict(method=method,cluster=c,N=N,K=K,n=n,a=a,fold=(a/n)/(K/N),p_value=float(hypergeom.sf(a-1,N,K,n))))
    for r,q in zip(membership,bh([r['p_value'] for r in membership])):r['q_global']=float(q)
    write(HERE/'tables/Spellman_membership.csv',membership)
    return tables,membership


def edges_for(genes,annotations):
    tags=sorted(set().union(*(annotations[g] for g in genes)))
    index={t:i for i,t in enumerate(tags)}
    B=np.zeros((len(genes),len(tags)),dtype=np.int32)
    for i,g in enumerate(genes):
        for t in annotations[g]:B[i,index[t]]=1
    A=B@B.T>0
    edges=np.column_stack(np.where(np.triu(A,1))).astype(np.int64)
    return edges


def similarity(ids,bases,universes,bp):
    lookup={g:i for i,g in enumerate(ids)}
    scenarios={name:(universes[name],bases[name]) for name in ['Process_coarse','Process_fine','Phase']}
    scenarios['GO_BP']=(sorted(g for g,v in bp.items() if v and g in lookup),bp)
    common=sorted(set(universes['Process_coarse'])&{g for g,v in bp.items() if v}&set(ids))
    scenarios['common_Process_coarse']=(common,bases['Process_coarse'])
    scenarios['common_GO_BP']=(common,bp)
    results=[];seedrows=[];per_cluster=[]
    for scenario,(genes,ann) in scenarios.items():
        idx=np.array([lookup[g] for g in genes]);edges=edges_for(genes,ann)
        baseline=len(edges)/(len(genes)*(len(genes)-1)/2)
        for mi,method in enumerate(METHODS):
            for seed in SEEDS:
                labels=np.load(HERE/'seed_runs'/method/f'seed_{seed}.npz')['labels'][idx]
                score=float(pair_score(labels,edges))
                seedrows.append(dict(scenario=scenario,method=method,seed=seed,N=len(genes),S=score,
                    analytic_null_mean=baseline,ratio=score/baseline if baseline else None))
            labels=np.load(HERE/'seed_runs'/method/'seed_42.npz')['labels'][idx]
            actual=float(pair_score(labels,edges))
            cache=HERE/'permutations'/scenario/f'{method}.npy';cache.parent.mkdir(parents=True,exist_ok=True)
            if cache.exists():null=np.load(cache)
            else:
                null=null_scores(labels,edges,N_PERM,20261127);np.save(cache,null)
            assert len(null)==N_PERM
            p=(1+int(np.count_nonzero(null>=actual-1e-15)))/(N_PERM+1)
            results.append(dict(scenario=scenario,method=method,N=len(genes),S=actual,null_mean=float(null.mean()),
                analytic_null_mean=baseline,ratio=actual/baseline,permutation_p=p,
                permutation_z=(actual-null.mean())/null.std(ddof=1) if null.std()>0 else None,
                permutations=N_PERM,seed=42))
            for c in range(5):
                group=np.flatnonzero(labels==c)
                hits=sum(labels[i]==c and labels[j]==c for i,j in edges)
                per_cluster.append(dict(scenario=scenario,method=method,cluster=c,n=len(group),
                    shared_pairs=int(hits),all_pairs=len(group)*(len(group)-1)//2,
                    S=hits/(len(group)*(len(group)-1)/2) if len(group)>1 else None))
            print(f'Similarity {scenario}/{method}: ratio={actual/baseline:.3f}, p={p:.5f}',flush=True)
    # Adjust all scenario x method tests; 20 seeds are stability replicates, not 20 independent experiments.
    for r,q in zip(results,bh([r['permutation_p'] for r in results])):r['permutation_q']=float(q)
    summary=[]
    for scenario in scenarios:
        for method in METHODS:
            values=[r['ratio'] for r in seedrows if r['scenario']==scenario and r['method']==method]
            summary.append(dict(scenario=scenario,method=method,n_seeds=len(values),ratio_mean=float(np.mean(values)),
                ratio_sd=float(np.std(values,ddof=1)),ratio_min=min(values),ratio_max=max(values)))
    write(HERE/'tables/similarity_permutations.csv',results)
    write(HERE/'tables/similarity_seeds.csv',seedrows)
    write(HERE/'tables/similarity_seed_summary.csv',summary)
    write(HERE/'tables/similarity_clusters_descriptive.csv',per_cluster)
    dump(HERE/'input/common_background.json',dict(genes=common,N=len(common)))
    return results,summary


def figures(tables,similarity_summary):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from textwrap import fill
    folder=HERE/'figures';folder.mkdir(exist_ok=True)
    for basis in ['Process_coarse','Phase']:
        records=tables[basis];terms=sorted({r['term'] for r in records})
        fig,axes=plt.subplots(7,1,figsize=(max(11,len(terms)*.8),15),sharex=True,layout='constrained')
        for ax,method in zip(axes,METHODS):
            for r in records:
                if r['method']!=method or r['phi'] is None:continue
                x=terms.index(r['term']);y=r['cluster']
                sc=ax.scatter(x,y,s=abs(r['phi'])*1250,c=[r['phi']],cmap='RdBu_r',vmin=-1,vmax=1,edgecolors='#777777',linewidth=.35)
                if r['q_phi']<.05:ax.text(x,y,'*',ha='center',va='center',fontsize=10,color='black')
            ax.set(yticks=range(5),ylabel='Cluster',title=method,xlim=(-.6,len(terms)-.4),ylim=(4.6,-.6))
            ax.grid(alpha=.15);ax.set_axisbelow(True)
        axes[-1].set_xticks(range(len(terms)),[fill(t,22) for t in terms],rotation=45,ha='right',fontsize=8)
        fig.suptitle(f'{basis}: binary cluster membership vs annotation (phi)\nArea = |phi|; red positive / blue negative; * two-sided Fisher BH q < 0.05',fontsize=13)
        fig.colorbar(sc,ax=axes,label='Phi coefficient',shrink=.55,pad=.02)
        for v in [.1,.3,.5]:axes[0].scatter([],[],s=v*1250,c='none',edgecolors='gray',label=f'|phi|={v}')
        axes[0].legend(loc='upper left',bbox_to_anchor=(1.08,1),fontsize=8)
        fig.savefig(folder/f'{basis}_phi_balloon.png',dpi=180);fig.savefig(folder/f'{basis}_phi_balloon.svg');plt.close(fig)
        for method in METHODS:
            fig,ax=plt.subplots(figsize=(13,5.4))
            for r in records:
                if r['method']!=method or r['phi'] is None:continue
                x=terms.index(r['term']);y=r['cluster']
                sc=ax.scatter(x,y,s=abs(r['phi'])*1450,c=[r['phi']],cmap='RdBu_r',vmin=-1,vmax=1,edgecolors='gray',linewidth=.4)
                if r['q_phi']<.05:ax.text(x,y,'*',ha='center',va='center',fontsize=12)
            ax.set(yticks=range(5),ylabel='Expression cluster',xlim=(-1.0,len(terms)-.3),ylim=(4.6,-.6))
            ax.set_xticks(range(len(terms)),[fill(t,17) for t in terms],rotation=40,ha='right',fontsize=8)
            ax.grid(alpha=.15);ax.set_axisbelow(True)
            fig.suptitle(f'{method}: {basis} association\nArea = |phi|; red positive, blue negative; * two-sided BH q<0.05',fontsize=12)
            fig.subplots_adjust(left=.08,right=.89,bottom=.43,top=.83)
            fig.colorbar(sc,ax=ax,label='Phi coefficient',fraction=.035,pad=.025)
            fig.savefig(folder/f'{basis}_phi_{method}.png',dpi=180);fig.savefig(folder/f'{basis}_phi_{method}.svg');plt.close(fig)
    # All five clusters remain visible, including those with no significant GO.
    go=tables['GO_all']
    for method in METHODS:
        selected=[]
        for c in range(5):
            hit=sorted([r for r in go if r['method']==method and r['cluster']==c and r['significant_enrichment'] and r['GO_aspect']=='biological_process'],key=lambda r:(r['q_enrichment'],r['p_enrichment']))
            for r in hit[:3]:
                if r['term'] not in selected:selected.append(r['term'])
        fig,ax=plt.subplots(figsize=(11,max(3.2,len(selected)*.5)),layout='constrained')
        if not selected:
            ax.text(.5,.5,'No significant biological-process terms',ha='center',transform=ax.transAxes)
        for r in go:
            if r['method']!=method or r['term'] not in selected:continue
            sc=ax.scatter(r['cluster'],selected.index(r['term']),s=r['a']*20,
                c=[-np.log10(max(r['q_enrichment'],1e-300))],cmap='viridis',vmin=0,vmax=max(1,max(-np.log10(max(x['q_enrichment'],1e-300)) for x in go)),
                edgecolor='black' if r['significant_enrichment'] else '#bbbbbb',linewidth=.8)
        ax.set_xticks(range(5),[f'C{c}' for c in range(5)])
        ax.set(xlim=(-.6,4.6),xlabel='Expression cluster',title=f'{method}: GO biological-process enrichment\nArea = gene hits; black outline = BH q<0.05 and hits>=3')
        names={r['term']:r['term_name'] for r in go}
        ax.set_yticks(range(len(selected)),[fill(names[t]+' '+t,65) for t in selected],fontsize=8);ax.invert_yaxis()
        ax.set_ylim(len(selected)-.65,-.65)
        if selected:fig.colorbar(sc,ax=ax,label='-log10(enrichment BH q)')
        for n in [3,10,20]:ax.scatter([],[],s=n*20,c='none',edgecolor='gray',label=f'{n} hits')
        ax.legend(fontsize=8,loc='upper left',bbox_to_anchor=(1.03,1))
        fig.savefig(folder/f'GO_BP_{method}.png',dpi=180);fig.savefig(folder/f'GO_BP_{method}.svg');plt.close(fig)
    scenarios=['Process_coarse','GO_BP','common_Process_coarse','common_GO_BP']
    fig,axes=plt.subplots(2,2,figsize=(14,9),layout='constrained')
    for ax,scenario in zip(axes.flat,scenarios):
        data=[next(r for r in similarity_summary if r['scenario']==scenario and r['method']==m) for m in METHODS]
        ax.bar(range(7),[r['ratio_mean'] for r in data],yerr=[r['ratio_sd'] for r in data],capsize=3,color='#4477aa')
        ax.axhline(1,color='black',linestyle='--',lw=1)
        ax.set_xticks(range(7),METHODS,rotation=45,ha='right',fontsize=8)
        ax.set(title=scenario,ylabel='S / expected random S')
    fig.suptitle('Within-cluster shared annotation: 20 initialization seeds\nError bars = sample SD across algorithm seeds, not biological uncertainty')
    fig.savefig(folder/'seed_stability.png',dpi=180);fig.savefig(folder/'seed_stability.svg');plt.close(fig)


def report(stats,tables,membership,perm,seedsummary):
    lines=['# Biological validation：表达簇与功能标签','',
     '日期：2026-11-27。所有新增过程、结果和图位于本目录。用户确认的气球图定义为：表达簇成员指示变量与功能标签指示变量的 φ 相关系数。另提供 GO 富集气泡图。','',
     '## 1. 数据、映射与研究范围','',
     '复用 2026-11-24 的 GSE3635 全部 6214 基因 × 13 时间点 centered.csv、七种距离矩阵、K=5、n_init=5、max_iter=100；不重新选择距离参数或 K。主展示用原 seed=42；稳定性增加 seed=42–61 共 20 个种子，每种方法用相同种子集合。',
     '聚类只使用表达距离；Process、Function、phase、GO 均在聚类后连接。Numba 数组版保持 xqf 的初始 random.Random、求和顺序、并列处理及交替更新逻辑，并对原始源码小样本及七种全 6214 基因 seed42 结果作完全一致核对。','',
     f"原始 Spellman 表解析得到 800 个唯一 ID；Process 已知 {stats['process_known']}、unknown {stats['process_unknown']}、空白 {stats['process_blank']}。空 ORF 的 YEL077C 依据一致的 YPD/SGD 名称及 orf800 恢复，没有直接删成799。",'',
     'Process 保留原文和 Function，仅移除明确的 (putative) 并另保留 putative 标记；修正 metbolism/acod；按分号拆分，多标签保留。Function 只作为命中基因的解释文字。',
     '自定义十四大类采用指引草案及明确关键词，术语到类别映射保存在 input/process_mapping.csv。其他、未分类不当作共同功能，也不用于检验；映射不依据簇的显著性调整。我们此前已看过表达簇/GO结果，因此本分析诚实标为探索性，不宣称独立盲法验证或事前注册。',
     '“Transcription and RNA processing”包含原文 RNA 加工；分泌类包含 Golgi、内吞及液泡生物发生；允许多类，例如 DNA 损伤 checkpoint 同时可属于 DNA 修复和周期调控。这些是本项目归类，不是标准 GO 标签。','',
     '| 注释体系 | 实际可分析背景基因数 |','|---|---:|']
    for basis,n in stats['universes'].items():lines.append(f'| {basis} | {n} |')
    lines+=['',f"其中 GO 为此前下载的 UniProt {stats['GO_release']} 直接注释，仅覆盖 Spellman 名单；背景不包含未获得注释的其余表达基因。GO 不扩展祖先，不把 IEA 等计算证据称为实验证据。Process_coarse 排除没有任何可分类标签的记录；Process_fine 只要求原 Process 已知。",'',
     '## 2. 统计定义与校正','',
     '每个簇×标签形成 [[a,b],[c,d]]，全部计数在对应可注释背景中：a=簇内有标签，b=簇内无标签，c=簇外有标签，d=簇外无标签。φ=(ad−bc)/sqrt((a+b)(c+d)(a+c)(b+d))，等于两个0/1变量的 Pearson 相关；正值为共同出现，负值为偏少。边际常数导致未定义时留空，不伪造0。',
     '气球图面积=|φ|，颜色=φ，*表示双侧 Fisher 检验经 BH 校正 q_phi<0.05。它不是表达曲线 Pearson 相关图。相关方向是描述性关联，不解释为调控因果。',
     '富集使用单侧 Fisher greater（等价超几何右尾），fold=(a/n)/(K/N)。先按背景基因频数≥5且非全背景筛选术语，所有零命中组合都参与BH；显著要求 q_enrichment<0.05、a≥3、fold>1。',
     '每个注释体系分别将七方法×五簇×全部合格标签作为一个 BH 检验族；q_phi 与 q_enrichment 分开校正。不同体系为不同问题，不把其显著数量直接做优劣排名。',
     '整体 S 为同簇基因对中至少共享一个标签的比例；它按基因对加权，因此大簇权重更高。背景内打乱簇标签保持簇大小不变，10000 次，p=(1+#(S_perm>=S))/(10001)。置换随机种子固定20261127。所有42个“场景×方法”置换p一起BH；它检验标签与当前簇的关联，不是实验生物学重复。',
     '随机 S 的精确期望为背景内所有基因对共享标签的比例，与当前簇大小无关（只要有同簇对）。同时保存实际置换均值，并报告 S/该精确期望。20种子均值±样本SD只衡量算法初始化稳定性，不把20次当20份独立实验，不平均p值。',
     '各簇的S仅作描述，不额外宣称逐簇置换显著。标签丰富度和类别粒度会影响S，因此另外在同一基因交集上比较 Process_coarse 和 GO_BP；即便如此，也不能把S数值当作统一的生物学准确率。','',
     '## 3. 主气球图','',
     '![Process φ关联气球图](figures/Process_coarse_phi_balloon.png)','',
     '![Phase φ关联气球图](figures/Phase_phi_balloon.png)','',
     '## 4. 富集结果','',
     '| 注释体系 | 检验数 | 显著簇×标签数 |','|---|---:|---:|']
    for b,rows in tables.items():lines.append(f"| {b} | {len(rows)} | {sum(r['significant_enrichment'] for r in rows)} |")
    lines+=['','下表列出每种方法最显著的自定义 Process 大类、以及最显著 GO 术语。所有逐簇表、命中基因和 Function 在 tables/ 中；未显著并不代表没有功能。','',
     '| 方法 | 最显著 Process 大类（簇、命中、fold、q） | 最显著 GO（簇、命中、fold、q） |','|---|---|---|']
    for method in METHODS:
        cells=[]
        for b in ['Process_coarse','GO_all']:
            rr=sorted([r for r in tables[b] if r['method']==method and r['significant_enrichment']],key=lambda r:r['q_enrichment'])
            if rr:
                r=rr[0];cells.append(f"C{r['cluster']} {r['term_name']}; {r['a']}/{r['n']}; fold={r['fold']:.2f}; q={r['q_enrichment']:.3g}")
            else:cells.append('无达到阈值的结果')
        lines.append('| '+method+' | '+' | '.join(cells)+' |')
    for method in METHODS:lines+=['',f'### {method} 的 GO 生物过程','',f'![GO {method}](figures/GO_BP_{method}.png)']
    lines+=['','## 5. 整体功能一致性与稳定性','','![20种子稳定性](figures/seed_stability.png)','',
     '| 场景 | 方法 | 背景N | seed42 S比值 | 置换q | 20种子比值均值±SD |','|---|---|---:|---:|---:|---:|']
    for r in perm:
        s=next(v for v in seedsummary if v['scenario']==r['scenario'] and v['method']==r['method'])
        lines.append(f"| {r['scenario']} | {r['method']} | {r['N']} | {r['ratio']:.3f} | {r['permutation_q']:.4g} | {s['ratio_mean']:.3f} ± {s['ratio_sd']:.3f} |")
    lines+=['','## 6. GO 是否更好','',
     'GO 注释覆盖更广、ID与术语标准化、可保留证据来源，适合作为主要功能解释；Process 大类更简洁，适合作为直观可视化和补充核查。但不能因为GO检验更多或显著项更多就认定更好。',
     '公平比较见 common_Process_coarse 与 common_GO_BP，它们具有完全相同的基因背景；仍需注意自定义大类与 GO 多标签的粒度不同。应该综合覆盖率、标签可靠性、对随机初始化的稳定性及命中基因解释判断，不用任何一项替代模拟ARI/NMI。',
     '本次相同背景交集为362个基因。以预先保留的seed42标准DTW为例，Process大类 S比值约1.028、置换BH q约0.170，未达阈值；GO_BP S比值约1.255、q约0.000162。此例说明GO能提供这套粗分类未显示的关联支持，不证明GO是真值，也不代表标准DTW在所有种子下都无功能一致性，20种子的分布应一起看。',
     '基于740个可注释基因的覆盖、可追溯标准术语和同背景对照，建议GO作为主要功能解释，Process大类气球图作为简洁补充；phase及周期名单验证单独保留。不同注释体系给方法的相对排序不完全相同，不能选择最有利的标签体系宣布某方法最优。',
     '前次2026-11-24 GO结果采用最低背景频数3，本次按指引统一采用5，因此GO显著数量可能变化，不能据差异断言结果冲突。','',
     '## 7. 对老师 biological validation 要求的补充','',
     'tables/Spellman_membership.csv 对全部6214输入基因、其中793个Spellman命中做周期基因富集，跨35项BH。Phase气球图和Phase富集使用有阶段注释的匹配基因背景。它们分别回答“是不是周期基因”和“哪一阶段偏多”，功能/GO回答“做什么”。三部分共同提供支持，模拟benchmark仍须独立报告。','',
     '## 8. 复现','',
     '在当前目录执行：','',
     '```powershell','python -B run_seeds.py','python -B analyze.py','python -B verify.py','```','',
     '依赖 NumPy、SciPy、Numba、Matplotlib。所有源输入、距离矩阵和原组员代码只读；源文件与输出校验和、配置、检验核对、文件保护记录见 configuration.json、verification.json。20种子结果在 seed_runs/；10000次置换原始数组在 permutations/。',
     'Process、GO和phase属于外部支持标签；历史注释不完整、基因间相关、GO层级依赖、固定K、无生物学重复等限制仍然存在。未知不代表阴性，关联不等于因果，初始化SD不代表生物学不确定性。','',
     '参考：Spellman et al. (1998), https://doi.org/10.1091/mbc.9.12.3273；GO 背景选择 https://geneontology.org/docs/go-enrichment-analysis/；置换p值加一原则 https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.permutation_test.html','']
    (HERE/'分析过程与结果.md').write_text('\n'.join(lines),encoding='utf-8')


def main():
    before_path=HERE/'protected_files_before.json'
    before=json.loads(before_path.read_text(encoding='utf-8')) if before_path.exists() else protected()
    dump(before_path,before)
    seedcheck=json.loads((HERE/'seed_validation.json').read_text())
    ids=[r['gene_id'] for r in load_csv(OLD/'expression_input.csv')]
    bases,universes,bp,termnames,aspects,ref,stats=annotations(ids)
    config=dict(date='2026-11-27',methods=METHODS,seeds=SEEDS,k=5,n_init=5,min_background=5,min_hits=3,
         n_permutations=N_PERM,permutation_seed=20261127,mapping_sha256=stats['mapping_sha256'],
         old_configuration_sha256=sha(OLD/'configuration.json'),
         guidance='User supplied 1.md; exploratory mapping frozen before this run, not claimed blinded',
         phi='signed Pearson correlation of two binary membership variables; area abs(phi), color phi',
         correction='within each annotation basis, all methods x clusters x tested terms; two-sided and greater separate',
         similarity_correction='BH over all six scenarios x seven methods',
         versions={name:__import__(name).__version__ for name in ['numpy','scipy','numba','matplotlib']})
    old=HERE/'configuration.json'
    if old.exists():assert json.loads(old.read_text(encoding='utf-8'))==config
    dump(old,config)
    tables,membership=enrichment(ids,bases,universes,termnames,aspects)
    perm,seedsummary=similarity(ids,bases,universes,bp)
    figures(tables,seedsummary)
    report(stats,tables,membership,perm,seedsummary)
    after=protected();assert before==after,'A file outside this new directory changed'
    dump(HERE/'verification.json',dict(status='passed',external_files_unchanged=True,protected_file_count=len(before),
         seed_validation=seedcheck,phi_matches_binary_pearson=True,Fisher_matches_hypergeometric=True,
         GO_counts={b:dict(tests=len(r),significant=sum(v['significant_enrichment'] for v in r)) for b,r in tables.items()},
         files={str(p.relative_to(HERE)):sha(p) for p in HERE.rglob('*') if p.is_file() and p.name!='verification.json' and '.mplconfig' not in str(p)}))
    print(json.dumps(stats,ensure_ascii=False,indent=2),flush=True)
    print('Report complete',flush=True)


if __name__=='__main__':main()
