"""Re-use frozen full expression distances and vary k-medoids initialization."""
import csv,hashlib,importlib.util,json,sys,time
from pathlib import Path
sys.dont_write_bytecode=True
import numpy as np
from fast_cluster import fit

HERE=Path(__file__).resolve().parent
OLD=HERE.parent/'go_enrichment_20261124'
METHODS=['standard','time','time-dynamics','diagonal-time-dynamics','diagonal-time','euclidean','pearson']
SEEDS=list(range(42,62))


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def load_csv(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))


def run():
    folder=HERE/'seed_runs';folder.mkdir(parents=True,exist_ok=True)
    spec=importlib.util.spec_from_file_location('validation_original',OLD/'source_snapshot/xqf/xqf_clustering/k_medoids.py')
    original=importlib.util.module_from_spec(spec);spec.loader.exec_module(original)
    rng=np.random.default_rng(27)
    for n in [6,15,30]:
        X=rng.normal(size=(n,5));D=np.sqrt(((X[:,None]-X[None,:])**2).sum(axis=2))
        for seed in [0,42,61]:
            y,m,_=fit(D,3,seed)
            yr,mr=original.k_medoids(D.tolist(),3,seed)
            assert y.tolist()==yr and m.tolist()==mr
    for D in [np.zeros((8,8)),np.ones((8,8))-np.eye(8)]:
        y,m,_=fit(D,3);yr,mr=original.k_medoids(D.tolist(),3)
        assert y.tolist()==yr and m.tolist()==mr
    timings=[];checks=[]
    for method in METHODS:
        source=OLD/'clustering'/method
        matrix_file=source/'distance_matrix.npy'
        meta=json.loads((source/'run.json').read_text(encoding='utf-8'))
        assert sha(matrix_file)==meta['distance_sha256']
        assert sha(source/'assignments.csv')==meta['assignment_sha256']
        D=np.load(matrix_file)
        expected=load_csv(source/'assignments.csv')
        ids=[r['gene_id'] for r in expected]
        for seed in SEEDS:
            output=folder/method/f'seed_{seed}.npz';output.parent.mkdir(parents=True,exist_ok=True)
            if output.exists():
                saved=np.load(output);labels=saved['labels'];medoids=saved['medoids'];objective=float(saved['objective']);seconds=float(saved['seconds'])
            else:
                start=time.perf_counter();labels,medoids,objective=fit(D,5,seed);seconds=time.perf_counter()-start
                np.savez_compressed(output,labels=labels,medoids=medoids,objective=objective,seconds=seconds)
            assert len(labels)==6214 and set(labels)==set(range(5))
            if seed==42:
                assert labels.tolist()==[int(r['cluster_id']) for r in expected]
                assert [ids[i] for i in medoids]==meta['medoid_gene_ids']
                checks.append(method)
            timings.append(dict(method=method,seed=seed,seconds=seconds,objective=objective))
            print(f'{method}: seed {seed} complete ({seconds:.2f}s)',flush=True)
        del D
    (HERE/'seed_validation.json').write_text(json.dumps(dict(status='passed',small_exact_checks=11,
        full_6214_seed42_exact_matches=checks,seeds=SEEDS,n_init=5,k=5,max_iter=100,
        rule='Python random.Random sampled initial medoids; identical ordered sums/ties; no fastmath',runs=timings),indent=2),encoding='utf-8')


if __name__=='__main__':run()
