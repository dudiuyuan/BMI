"""Independent verification of enrichment, signed phi, FDR, and pair scores."""
import itertools,json,csv,sys
from pathlib import Path
sys.dont_write_bytecode=True
import numpy as np
from scipy.stats import false_discovery_control,fisher_exact
from fast_cluster import pair_score,null_scores

HERE=Path(__file__).resolve().parent


def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))


def main():
    ids=[r['gene_id'] for r in read(HERE.parent/'go_enrichment_20261124/expression_input.csv')]
    idx={g:i for i,g in enumerate(ids)}
    annotation=read(HERE/'input/spellman_annotations.csv')
    assert len(annotation)==800 and len({r['gene_id'] for r in annotation})==800
    assert next(r for r in annotation if r['gene_id']=='YEL077C')['mapping_rule']!='ORF'
    assert sum(bool(r['Process_raw']) and r['Process_raw'].strip().lower()!='unknown' for r in annotation)==398
    split=lambda text:set(filter(None,text.split(';')))
    bases={'Process_coarse':{r['gene_id']:split(r['categories']) for r in annotation},
           'Process_fine':{r['gene_id']:split(r['terms']) for r in annotation},
           'Phase':{r['gene_id']:{r['phase']} for r in annotation}}
    go={};bp={}
    for r in read(HERE.parent/'go_annotation/results/GO_long.csv'):
        if r['GO_code'] in ['GO:0008150','GO:0003674','GO:0005575']:continue
        go.setdefault(r['gene_id'],set()).add(r['GO_code'])
        if r['GO_aspect']=='biological_process':bp.setdefault(r['gene_id'],set()).add(r['GO_code'])
    bases['GO_all']=go
    checked=0
    for basis,ann in bases.items():
        bg={g for g,tags in ann.items() if tags and g in idx}
        rows=read(HERE/'tables'/f'{basis}_all_tests.csv')
        for name,pname in [('q_phi','p_phi'),('q_enrichment','p_enrichment')]:
            expected=false_discovery_control(np.array([float(r[pname]) for r in rows]),method='bh')
            assert np.allclose(expected,[float(r[name]) for r in rows],rtol=1e-12,atol=1e-15)
        labels={m:np.load(HERE/'seed_runs'/m/'seed_42.npz')['labels'] for m in {r['method'] for r in rows}}
        for r in rows:
            group={g for g in bg if labels[r['method']][idx[g]]==int(r['cluster'])}
            term={g for g in bg if r['term'] in ann[g]}
            hits=group&term
            assert split(r['hit_genes'])==hits
            a,b,c,d=[int(r[x]) for x in 'abcd']
            assert a==len(hits) and a+b==len(group) and a+c==len(term) and a+b+c+d==len(bg)
            if r['phi']:
                denom=np.sqrt(float(len(group)*(len(bg)-len(group))*len(term)*(len(bg)-len(term))))
                assert np.isclose(float(r['phi']),(a*d-b*c)/denom)
            checked+=1
        for r in rows[::max(1,len(rows)//15)]:
            table=[[int(r['a']),int(r['b'])],[int(r['c']),int(r['d'])]]
            assert np.isclose(fisher_exact(table,alternative='two-sided').pvalue,float(r['p_phi']))
            assert np.isclose(fisher_exact(table,alternative='greater').pvalue,float(r['p_enrichment']))
    # Exhaustive toy assignment check: random cluster labels retain group sizes.
    labels=np.array([0,0,1,1,2,2]);edges=np.array([[0,1],[0,2],[2,3]],dtype=np.int64)
    values=[pair_score(np.array(p),edges) for p in set(itertools.permutations(labels))]
    assert np.isclose(np.mean(values),len(edges)/15)
    assert np.array_equal(null_scores(labels,edges,100,27),null_scores(labels,edges,100,27))
    perm=read(HERE/'tables/similarity_permutations.csv')
    expected=false_discovery_control(np.array([float(r['permutation_p']) for r in perm]))
    assert np.allclose(expected,[float(r['permutation_q']) for r in perm])
    common=set(json.loads((HERE/'input/common_background.json').read_text())['genes'])
    for r in perm:
        name=r['scenario']
        ann=bases.get(name,bp if name in ['GO_BP','common_GO_BP'] else bases['Process_coarse'])
        genes=sorted(common if name.startswith('common_') else {g for g,t in ann.items() if t and g in idx})
        labels=np.load(HERE/'seed_runs'/r['method']/'seed_42.npz')['labels']
        num=den=0
        for i,g in enumerate(genes):
            for h in genes[i+1:]:
                if labels[idx[g]]==labels[idx[h]]:
                    den+=1;num+=bool(ann[g]&ann[h])
        assert np.isclose(num/den,float(r['S']))
        null=np.load(HERE/'permutations'/name/(r['method']+'.npy'))
        assert len(null)==10000
        p=(1+sum(null>=float(r['S'])-1e-15))/10001
        assert np.isclose(p,float(r['permutation_p']))
    summary=dict(status='passed',all_contingency_and_hit_rows_checked=checked,
        BH_matches_scipy=True,exact_toy_permutation_expectation=True,
        all_42_S_values_recomputed_by_independent_pair_loops=True,permutations_reproducible=True)
    (HERE/'independent_verification.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
