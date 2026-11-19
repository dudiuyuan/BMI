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
from run_dtw import read_matrix,prepare_expression
from support.diagnostic_tools import average_linkage,k_medoids,evaluate

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


def build_report(folder,k,truth=None,label_column='true_cluster',clustering=None):
    folder=folder.resolve()
    meta,source,ids,t,x,d=load_experiment(folder)
    # Cluster the saved distances before attaching truth; no distance matrix is recalculated.
    clustering=clustering or meta.get('clustering','average')
    if clustering not in ('average','kmedoids'):
        raise ValueError('Unknown clustering algorithm')
    assignment,groups=(k_medoids if clustering=='kmedoids' else average_linkage)(d,k)
    medoids=[min(g,key=lambda i:(sum(d[i][j] for j in g),i)) for g in groups]
    truth=truth or (source.parent.parent/'ground_truth.csv' if source.parent.name=='processed' else source.parent/'ground_truth.csv')
    labels=read_labels(truth,ids,label_column);cfg=Config(**meta['config']);n=len(ids)
    r={'name':folder.name,'config':asdict(cfg),'matrix':d,'assignment':assignment,'clusters':groups,
       'medoids':medoids,'metrics':evaluate(labels,assignment,d),'composition':[dict(Counter(labels[i] for i in g)) for g in groups],
       'nearest':[],'paths':{},'config_file':str(folder/'run.json')}
    kernel_hash=hashlib.sha256((ROOT/'dtw.py').read_bytes()).hexdigest()
    paths_allowed=meta.get('kernel_sha256',kernel_hash)==kernel_hash
    warning='' if paths_allowed else '核心代码已改变：仅展示原实验距离和分组，未重算匹配路径。'
    for i in range(n):
        ranked=sorted((j for j in range(n) if j!=i),key=lambda j:(d[i][j],j))
        same=next((j for j in ranked if labels[j]==labels[i]),None)
        other=next((j for j in ranked if labels[j]!=labels[i]),None)
        r['nearest'].append({'top':ranked[:5],'same':same,'other':other})
        targets=set(ranked[:5]+[medoids[assignment[i]]]+[j for j in (same,other) if j is not None])
        if not paths_allowed:continue
        for j in targets:
            a,b=sorted((i,j));key=f'{a}:{b}'
            if key in r['paths']:continue
            result=align(x[a],x[b],cfg,times_x=t,times_y=t)
            if abs(result.distance-d[a][b])>1e-10:
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
             'original_n':n,'kernel':str(ROOT/'dtw.py'),'kernel_sha256':kernel_hash,
             'input_sha256':meta['input_sha256'],'truth_sha256':hashlib.sha256(truth.read_bytes()).hexdigest(),
             'clustering':clustering,
             'clustering_name':'k-medoids（PAM）' if clustering=='kmedoids' else '平均连接层次聚类',
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
    p.add_argument('--clustering',choices=['kmedoids','average'],
                   help='Override downstream clustering; default follows run.json (legacy runs: average)')
    args=p.parse_args()
    try:print(build_report(args.experiment,args.k,args.truth,args.label_column,args.clustering))
    except (ValueError,OSError,KeyError,TypeError,OverflowError,IndexError) as exc:p.error(str(exc))


if __name__=='__main__':main()
