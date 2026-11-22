import csv
import tempfile
import unittest
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from support.diagnostic_tools import average_linkage,k_medoids,evaluate
from diagnose import read_labels

from diagnose import load_experiment,build_report
from dtw import preprocessing_config,prepare_expression
from dtw import Config,align
from unittest.mock import patch
import hashlib
import json


class DiagnosticTests(unittest.TestCase):
    def test_pam_matches_exhaustive_optimum_on_hand_checked_points(self):
        from itertools import combinations
        points=[0,1,2,8,9,14]
        matrix=[[abs(a-b) for b in points] for a in points]
        assignment,groups=k_medoids(matrix,2)
        self.assertEqual(groups,[[0,1,2],[3,4,5]])
        actual=sum(min(sum(matrix[i][m] for i in g) for m in g) for g in groups)
        optimal=min(sum(min(matrix[i][m] for m in medoids) for i in range(6))
                    for medoids in combinations(range(6),2))
        self.assertEqual(actual,optimal)

    def test_pam_ties_invalid_inputs_and_one_swap_local_optimum(self):
        import random
        rng=random.Random(19)
        for n,k in [(8,3),(9,2),(7,1)]:
            matrix=[[0.0]*n for _ in range(n)]
            for i in range(n):
                for j in range(i):matrix[i][j]=matrix[j][i]=rng.random()
            assignment,groups=k_medoids(matrix,k)
            self.assertEqual((assignment,groups),k_medoids(matrix,k))
            medoids=[min(g,key=lambda m:sum(matrix[i][m] for i in g)) for g in groups]
            objective=lambda centers:sum(min(matrix[i][m] for m in centers) for i in range(n))
            cost=objective(medoids)
            for pos in range(k):
                for candidate in set(range(n))-set(medoids):
                    swapped=medoids[:];swapped[pos]=candidate
                    self.assertGreaterEqual(objective(swapped)+1e-10,cost)
        for k in (1,3,5):
            assignment,groups=k_medoids([[0]*5 for _ in range(5)],k)
            self.assertEqual(len(groups),k)
            self.assertEqual(sorted(i for g in groups for i in g),list(range(5)))
        for matrix,k in [([[0,1],[2,0]],1),([[0]],2),([[1]],1)]:
            with self.assertRaises(ValueError):k_medoids(matrix,k)

    def test_report_replays_normalization_and_clustering(self):
        import re
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=root/'input.csv'
            source.write_text('gene_id,t0,t1,t2\na,1,2,3\nb,12,14,16\nc,0,0,0\n')
            spec=preprocessing_config(.1)
            x=prepare_expression([[1,2,3],[12,14,16],[0,0,0]],spec)
            with (root/'distance_matrix.csv').open('w',newline='') as f:
                writer=csv.writer(f);writer.writerow(['gene_id','a','b','c'])
                for i,gene in enumerate(['a','b','c']):
                    writer.writerow([gene]+[align(x[i],y).distance for y in x])
            (root/'run.json').write_text(json.dumps(dict(mode='matrix',input=str(source),
                input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),config={},
                preprocessing=spec,clustering='kmedoids')))
            truth=root/'ground_truth.csv'
            truth.write_text('gene_id,true_cluster\na,shape\nb,shape\nc,flat\n')
            report=build_report(root,2)
            payload=json.loads(re.search(r'<script id="report-data" type="application/json">(.*?)</script>',
                report.read_text(),re.S).group(1))
            self.assertEqual(payload['values'],x)
            self.assertEqual(payload['clustering'],'kmedoids')
            self.assertEqual(payload['runs'][0]['metrics']['ari'],1)
            self.assertEqual(payload['warning'],'')
            self.assertGreater(len(payload['runs'][0]['paths']),0)
            # An outdated code fingerprint no longer hides paths: each path is
            # kept because its recomputed cost matches the saved distance.
            meta=json.loads((root/'run.json').read_text());meta['kernel_sha256']='0'*64
            (root/'run.json').write_text(json.dumps(meta))
            payload=json.loads(re.search(r'<script id="report-data" type="application/json">(.*?)</script>',
                build_report(root,2).read_text(),re.S).group(1))
            self.assertEqual(payload['warning'],'')
            self.assertGreater(len(payload['runs'][0]['paths']),0)
            self.assertEqual(load_experiment(root)[4],x)

    def test_saved_experiment_is_read_without_recomputing(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp)
            source=folder/'centered.csv'
            source.write_text('gene_id,t0,t1\na,0,1\nb,1,0\n')
            (folder/'run.json').write_text(json.dumps({'mode':'matrix','input':str(source),
                'input_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'config':{}}))
            (folder/'distance_matrix.csv').write_text('gene_id,a,b\na,0,1.41421356237\nb,1.41421356237,0\n')
            with patch('diagnose.align',side_effect=AssertionError('must not recalculate matrix')):
                meta,path,ids,times,values,matrix=load_experiment(folder)
            self.assertEqual(ids,['a','b'])
            self.assertEqual(matrix[0][1],1.41421356237)
            source.write_text('changed')
            with self.assertRaises(ValueError):load_experiment(folder)

    def test_hand_checked_average_linkage(self):
        # First merge 0/1 at 1; then 2/3 at 2; cross distances are all larger.
        d=[[0,1,9,8],[1,0,7,6],[9,7,0,2],[8,6,2,0]]
        assignment,groups=average_linkage(d,2)
        self.assertEqual(groups,[[0,1],[2,3]])
        self.assertEqual(assignment,[0,0,1,1])
        metrics=evaluate(['a','a','b','b'],assignment,d)
        self.assertEqual(metrics['ari'],1)
        self.assertEqual(metrics['purity'],1)
        self.assertEqual(metrics['nearest_same_type'],1)

    def test_ari_is_permutation_invariant_and_can_be_negative(self):
        d=[[0,1,2,3],[1,0,3,2],[2,3,0,1],[3,2,1,0]]
        labels=['a','a','b','b']
        self.assertEqual(evaluate(labels,[5,5,2,2],d)['ari'],1)
        self.assertAlmostEqual(evaluate(labels,[0,1,0,1],d)['ari'],-.5)

    def test_truth_join_uses_id_and_kept_subset(self):
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp)/'truth.csv'
            p.write_text('gene_id,true_cluster\nc,b\na,a\nb,a\n')
            self.assertEqual(read_labels(p,['b','c'],'true_cluster'),['a','b'])
            with self.assertRaises(ValueError):read_labels(p,['missing'],'true_cluster')

    def test_reject_bad_matrix(self):
        with self.assertRaises(ValueError):average_linkage([[0,1],[2,0]],1)
        with self.assertRaises(ValueError):average_linkage([[0]],2)

    def test_linkage_matches_direct_group_averages(self):
        import random
        rng=random.Random(17)
        d=[[0.0]*8 for _ in range(8)]
        for i in range(8):
            for j in range(i+1,8):d[i][j]=d[j][i]=rng.random()
        groups=[[i] for i in range(8)]
        while len(groups)>3:
            pairs=[(sum(d[i][j] for i in a for j in b)/(len(a)*len(b)),ia,ib)
                   for ia,a in enumerate(groups) for ib,b in enumerate(groups) if ia<ib]
            _,a,b=min(pairs)
            merged=groups[a]+groups[b]
            groups=[g for i,g in enumerate(groups) if i not in (a,b)]+[merged]
        _,actual=average_linkage(d,3)
        self.assertEqual({frozenset(g) for g in actual},{frozenset(g) for g in groups})


if __name__=='__main__':unittest.main()
