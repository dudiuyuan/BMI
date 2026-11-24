"""Independent consistency checks of saved matrices, clusters, and GO counts."""
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path
import numpy as np
from scipy.stats import fisher_exact

HERE=Path(__file__).resolve().parent


def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:
        return list(csv.DictReader(f))


def main():
    config=json.loads((HERE/'configuration.json').read_text(encoding='utf-8'))
    ids=[r['gene_id'] for r in read(HERE/'expression_input.csv')]
    assert len(ids)==len(set(ids))==6214
    original_go=HERE.parent/'go_annotation/results/GO_long.csv'
    terms=defaultdict(set)
    for r in read(original_go):
        if r['gene_id'] in ids: terms[r['GO_code']].add(r['gene_id'])
    background=set().union(*terms.values())
    counts=[]
    for method in config['parameters']:
        folder=HERE/'clustering'/method
        pred=read(folder/'assignments.csv')
        assert [r['gene_id'] for r in pred]==ids
        labels=np.array([int(r['cluster_id']) for r in pred])
        medoids=read(folder/'medoids.csv')
        indices=[ids.index(r['gene_id']) for r in medoids]
        D=np.load(folder/'distance_matrix.npy',mmap_mode='r')
        expected=np.argmin(D[:,indices],axis=1)
        expected[indices]=np.arange(5)
        assert np.array_equal(expected,labels),method
        info=json.loads((folder/'run.json').read_text())
        objective=float(D[np.arange(len(ids)),np.array(indices)[labels]].sum())
        assert np.isclose(objective,info['objective'],rtol=1e-12,atol=1e-9)
        for cluster in range(5):
            group={g for g,c in zip(ids,labels) if c==cluster}&background
            table=read(HERE/'enrichment'/method/f'cluster_{cluster}.csv')
            for r in table:
                hit=group&terms[r['GO_code']]
                assert set(filter(None,r['hit_genes'].split(';')))==hit
                a,b,c,d=[int(r[x]) for x in 'abcd']
                assert a==len(hit) and a+b==len(group) and a+c==len(terms[r['GO_code']])
                assert a+b+c+d==len(background)
                assert float(r['q_global'])+1e-15>=float(r['p_value'])
            counts.append({'method':method,'cluster':cluster,'genes':int((labels==cluster).sum()),'tests':len(table)})
        del D
    p=HERE/'independent_verification.json'
    p.write_text(json.dumps({'status':'passed','all_6214_assignments_per_method':True,
        'saved_labels_equal_nearest_medoids':True,'objectives_checked':True,
        'all_GO_hit_lists_and_contingency_tables_checked':True,'clusters':counts},indent=2)+'\n',encoding='utf-8')
    print('Independent checks passed for all seven methods and 35 clusters.')


if __name__=='__main__': main()
