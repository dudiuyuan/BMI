#!/usr/bin/env python3
"""Visualize one existing experiment. Never launch new parameter experiments."""
import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from dataclasses import asdict
from dtw import Config,align
from dtw import read_matrix,prepare_expression
from support.diagnostic_tools import average_linkage,k_medoids,evaluate
from support.xqf_adapter import xqf_k_medoids,DEFAULT_SEED

ROOT=Path(__file__).resolve().parent


def read_labels(path,ids,column):
    with path.open(encoding='utf-8-sig',newline='') as f:
        reader=csv.DictReader(f)
        if not {'gene_id',column} <= set(reader.fieldnames or []):
            raise ValueError(f'Truth table needs gene_id and {column}')
        mapping={}
        for row in reader:
            if not row['gene_id'] or not row[column] or row['gene_id'] in mapping:
                raise ValueError('Empty or duplicate ID/type in truth table')
            mapping[row['gene_id']]=row[column]
    if set(ids)-mapping.keys():raise ValueError('Some experiment genes have no truth label')
    return [mapping[g] for g in ids]


def load_experiment(folder):
    meta=json.loads((folder/'run.json').read_text(encoding='utf-8'))
    if meta.get('mode')!='matrix':raise ValueError('Choose an experiment containing distance_matrix.csv (matrix mode)')
    source=Path(meta['input'])
    if not source.is_absolute():source=(folder/source).resolve()
    if hashlib.sha256(source.read_bytes()).hexdigest()!=meta['input_sha256']:
        raise ValueError('Input file changed since this experiment; refusing to mix results')
    all_ids,t,all_x=read_matrix(source)
    with (folder/'distance_matrix.csv').open(encoding='utf-8-sig',newline='') as f:
        rows=list(csv.reader(f))
    if not rows or rows[0][0]!='gene_id':raise ValueError('Invalid distance matrix header')
    ids=rows[0][1:]
    if len(ids)<2 or len(set(ids))!=len(ids):raise ValueError('Need at least two unique gene IDs')
    if len(ids)>500:raise ValueError('Visual diagnosis supports up to 500 experiment genes')
    if len(rows)-1!=len(ids) or [r[0] for r in rows[1:]]!=ids:
        raise ValueError('Matrix row/column gene order does not match')
    d=[[float(v) for v in row[1:]] for row in rows[1:]]
    all_x=prepare_expression(all_x,meta.get('preprocessing'))
    lookup=dict(zip(all_ids,all_x))
    if set(ids)-lookup.keys():raise ValueError('Matrix contains IDs absent from original input')
    return meta,source,ids,t,[lookup[g] for g in ids],d


DISTANCE_NAMES={'dtw':'DTW','euclidean':'欧氏距离（Euclidean）','pearson':'皮尔森相关距离（1 − r）'}
CLUSTERING_NAMES={'kmedoids':'k-medoids（PAM）','average':'平均连接层次聚类','xqf':'k-medoids（徐秋凡版，随机多起点）'}


def build_report(folder,k,truth=None,label_column='true_cluster',clustering=None,seed=None):
    folder=folder.resolve()
    meta,source,ids,t,x,d=load_experiment(folder)
    # Cluster the saved distances before attaching truth; no distance matrix is recalculated.
    clustering=clustering or meta.get('clustering','average')
    if clustering not in CLUSTERING_NAMES:
        raise ValueError('Unknown clustering algorithm')
    seed=None if clustering!='xqf' else (meta.get('clustering_seed',DEFAULT_SEED) if seed is None else seed)
    if clustering=='xqf':assignment,groups=xqf_k_medoids(d,k,seed=seed)
    else:assignment,groups=(k_medoids if clustering=='kmedoids' else average_linkage)(d,k)
    medoids=[min(g,key=lambda i:(sum(d[i][j] for j in g),i)) for g in groups]
    truth=truth or (source.parent.parent/'ground_truth.csv' if source.parent.name=='processed' else source.parent/'ground_truth.csv')
    labels=read_labels(truth,ids,label_column);n=len(ids)
    # Experiments without a distance_method field were all written by the DTW runner.
    distance_method=meta.get('distance_method','dtw')
    if distance_method not in DISTANCE_NAMES:raise ValueError('Unknown distance method in run.json')
    cfg=Config(**meta['config']) if distance_method=='dtw' else None
    r={'name':folder.name,'config':asdict(cfg) if cfg else {'distance_method':distance_method},'matrix':d,'assignment':assignment,'clusters':groups,
       'medoids':medoids,'metrics':evaluate(labels,assignment,d),'composition':[dict(Counter(labels[i] for i in g)) for g in groups],
       'nearest':[],'paths':{},'config_file':str(folder/'run.json')}
    # Each stored DTW path is recomputed and kept only if its cost reproduces the saved
    # distance, so unchanged results still show paths after unrelated code edits.
    # Euclidean and Pearson compare equal time points directly and have no path.
    warning=''
    pruning=meta.get('pruning') or {}
    tau=pruning.get('threshold') if pruning.get('pairs') else None
    for i in range(n):
        ranked=sorted((j for j in range(n) if j!=i),key=lambda j:(d[i][j],j))
        same=next((j for j in ranked if labels[j]==labels[i]),None)
        other=next((j for j in ranked if labels[j]!=labels[i]),None)
        r['nearest'].append({'top':ranked[:5],'same':same,'other':other})
        targets=set(ranked[:5]+[medoids[assignment[i]]]+[j for j in (same,other) if j is not None])
        if cfg is None:continue
        for j in targets:
            a,b=sorted((i,j));key=f'{a}:{b}'
            if key in r['paths']:continue
            result=align(x[a],x[b],cfg,times_x=t,times_y=t)
            # With --prune pairs, distances above the threshold were saved as the threshold.
            truncated=tau is not None and abs(d[a][b]-tau)<=1e-12 and result.distance>=tau-1e-10
            if abs(result.distance-d[a][b])>1e-10 and not truncated:
                warning='当前代码的部分对齐费用与原实验不一致；这些基因对未展示路径，原距离矩阵保持不变。'
                continue
            r['paths'][key]={'path':result.path,'diagnostics':result.diagnostics,'cost':result.cost}
    r['mixed_nearest']=[i for i in range(n) if labels[i]!=labels[r['nearest'][i]['top'][0]]]
    early=[i for i,l in enumerate(labels) if l=='early_peak'];late=[i for i,l in enumerate(labels) if l=='late_peak']
    if early and late:r['metrics']['early_late_cocluster_rate']=sum(assignment[i]==assignment[j] for i in early for j in late)/(len(early)*len(late))
    preprocessing=meta.get('preprocessing',{'method':'none','version':1})
    expression_description=(f"逐基因减均值，再除以 max(自身标准差, {preprocessing['scale_floor']:g})"
                            if preprocessing['method']=='center_scale' else '输入表中的表达值，未经额外缩放')
    payload={'ids':ids,'times':t,'values':x,'labels':labels,'types':sorted(set(labels)),'k':k,'runs':[r],
             'input':str(source),'experiment':str(folder),'truth':str(truth.resolve()),'label_column':label_column,
             'original_n':n,'distance_method':distance_method,'distance_name':DISTANCE_NAMES[distance_method],
             'input_sha256':meta['input_sha256'],'truth_sha256':hashlib.sha256(truth.read_bytes()).hexdigest(),
             'clustering':clustering,'clustering_seed':seed,'pruning':pruning,
             'clustering_name':CLUSTERING_NAMES[clustering],
             'preprocessing':preprocessing,'expression_description':expression_description,
             'warning':warning}
    raw=json.dumps(payload,ensure_ascii=False,allow_nan=False,separators=(',',':'))
    template=(ROOT/'support/report_template.html').read_text(encoding='utf-8')
    # The sole output is the report, next to the selected experiment. No extra data or configs.
    report=folder/'diagnose.html'
    report.write_text(template.replace('__REPORT_DATA__',raw.replace('<','\\u003c')),encoding='utf-8')
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('experiment',type=Path,help='Existing result folder containing run.json and distance_matrix.csv')
    p.add_argument('--k',type=int,required=True,help='Number of clusters for this diagnosis')
    p.add_argument('--truth',type=Path)
    p.add_argument('--label-column',default='true_cluster')
    p.add_argument('--clustering',choices=['kmedoids','average','xqf'],
                   help='Override downstream clustering; default follows run.json (legacy runs: average). '
                        'xqf = Xu Qiufan k-medoids from xqf/xqf_clustering')
    p.add_argument('--seed',type=int,help='Random seed for --clustering xqf (default: run.json or 42)')
    args=p.parse_args()
    try:print(build_report(args.experiment,args.k,args.truth,args.label_column,args.clustering,args.seed))
    except (ValueError,OSError,KeyError,TypeError,OverflowError,IndexError) as exc:p.error(str(exc))


if __name__=='__main__':main()
