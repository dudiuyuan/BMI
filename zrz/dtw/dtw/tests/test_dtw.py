import sys
import unittest
import json
import subprocess
import tempfile
from pathlib import Path
from itertools import product
from statistics import pstdev
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from dtw import Config, align
from run_dtw import read_matrix,preprocessing_config,prepare_expression


def paths(m,n,i=0,j=0):
    if (i,j)==(m-1,n-1):
        yield [(i,j)]
    for di,dj in [(1,1),(0,1),(1,0)]:
        if i+di<m and j+dj<n:
            for suffix in paths(m,n,i+di,j+dj):
                yield [(i,j)]+suffix


class Tests(unittest.TestCase):
    def test_against_exhaustive_paths(self):
        # Independent enumeration verifies affine state history and unequal lengths.
        x,y=[0,2,1],[0,1]
        for a,b,lam in product([0,.2],[0,.05],[0,2]):
            cfg=Config(open_penalty=a,extend_penalty=b,time_penalty=lam)
            costs=[]
            for p in paths(len(x),len(y)):
                c=sum((x[i]-y[j])**2+lam*abs(i/2-j)**2 for i,j in p)
                last=None
                for (i,j),(k,l) in zip(p,p[1:]):
                    move=(k-i,l-j)
                    if move!=(1,1):
                        c+=b if move==last else a
                    last=move
                costs.append(c)
            r=align(x,y,cfg)
            self.assertAlmostEqual(r.cost,min(costs))
            self.assertAlmostEqual(align(x,y,cfg,return_path=False).cost,r.cost)

    def test_symmetry_and_time_units(self):
        x,y=[0,2,0],[1,0,1]
        cfg=Config(time_penalty=2,dynamics=.3,open_penalty=.4,extend_penalty=.1)
        r=align(x,y,cfg,times_x=[0,10,30],times_y=[0,10,30])
        self.assertAlmostEqual(r.cost,align(y,x,cfg,times_x=[0,1,3],times_y=[0,1,3]).cost)

    def test_tiny_positive_window_is_euclidean(self):
        x,y=[1,4,2],[3,0,2]
        r=align(x,y,Config(window=1e-12))
        self.assertEqual(r.cost,20)
        self.assertEqual(r.path,[(0,0),(1,1),(2,2)])

    def test_linear_is_affine_equal_rates(self):
        x,y=[0,0,3,0],[0,3,0,0]
        self.assertEqual(align(x,y,Config(step_penalty=.3)).cost,
                         align(x,y,Config(open_penalty=.3,extend_penalty=.3)).cost)

    def test_invalid_inputs(self):
        for x in [[],[float('nan')],[float('inf')]]:
            with self.assertRaises(ValueError): align(x,[0])
        with self.assertRaises(ValueError): Config(step_penalty=-1)
        with self.assertRaises(ValueError): align([0,1],[0,1],times_x=[0,0],times_y=[0,1])
        with self.assertRaises(ValueError): align([0,1],[0,1,2],Config(window=1e-12))

    def test_constant_and_singletons(self):
        self.assertEqual(align([0]*4,[0]*4).distance,0)
        self.assertEqual(align([2],[5]).distance,3)

    def test_zero_switches_disable_and_allow_unused_dynamics(self):
        x,y=[0,0,2,0],[0,2,0,0]
        self.assertEqual(align(x,y,Config(window=0)).cost,align(x,y,Config(window=None)).cost)
        self.assertEqual(align(x,y,Config(dynamics=5,time_penalty=0)).cost,align(x,y).cost)
        self.assertEqual(Config(open_penalty=None,extend_penalty=None).open_penalty,0)
        Config(open_penalty=1)  # Independently configurable; extension can be free.

    def test_combined_transition_penalties_against_exhaustive_paths(self):
        x,y=[0,2,1],[0,1]
        cfg=Config(step_penalty=.1,open_penalty=.2,extend_penalty=.05,time_penalty=1)
        costs=[]
        for p in paths(len(x),len(y)):
            c=sum((x[i]-y[j])**2+abs(i/2-j)**2 for i,j in p)
            last=None
            for (i,j),(k,l) in zip(p,p[1:]):
                move=(k-i,l-j)
                if move!=(1,1): c+=.1+(.05 if move==last else .2)
                last=move
            costs.append(c)
        self.assertAlmostEqual(align(x,y,cfg).cost,min(costs))

    def test_cli_direct_parameters_normalization_and_legacy_replay(self):
        script=Path(__file__).resolve().parents[1]/'run_dtw.py'
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            source=root/'input.csv'
            source.write_text('gene_id,t0,t1\na,10,10\nb,13,13\n')
            def run(name,*options):
                output=root/name
                result=subprocess.run([sys.executable,'-B',str(script),'--input',str(source),
                    '--output',str(output),'--pair','a','b',*options],capture_output=True,text=True)
                self.assertEqual(result.returncode,0,result.stderr)
                return json.loads((output/'run.json').read_text())
            default=run('default')
            self.assertEqual(default['config']['window'],0)
            for key in ('step_penalty','time_penalty','dynamics','open_penalty','extend_penalty'):
                self.assertEqual(default['config'][key],0)
            self.assertNotIn('method',default)
            self.assertEqual(default['cost'],0)  # Both constant curves center to zero.
            self.assertEqual(default['preprocessing'],preprocessing_config(.1))
            self.assertEqual(default['clustering'],'kmedoids')
            unscaled=run('unscaled','--scale-floor','0')
            self.assertEqual(unscaled['cost'],18)
            combined=run('combined','--step-penalty','0.1','--open-penalty','0.2',
                '--extend-penalty','0.05','--time-penalty','1','--dynamics','0.3',
                '--window','0.5','--power','3')
            self.assertEqual(combined['config'],dict(window=.5,step_penalty=.1,
                open_penalty=.2,extend_penalty=.05,time_penalty=1,dynamics=.3,power=3))
            legacy=root/'legacy.json'
            legacy.write_text(json.dumps({'method':'affine','config':dict(window=None,
                open_penalty=.2,extend_penalty=.05)}))
            loaded=run('loaded','--config-file',str(legacy))
            self.assertEqual(loaded['config']['open_penalty'],.2)
            self.assertEqual(loaded['config']['window'],0)
            self.assertEqual(loaded['cost'],18)
            self.assertEqual(loaded['preprocessing']['method'],'none')
            self.assertEqual(loaded['clustering'],'average')
            replay=run('replayed','--config-file',str(root/'default/run.json'),'--scale-floor','0')
            self.assertEqual(replay['cost'],default['cost'])
            self.assertEqual(replay['preprocessing'],default['preprocessing'])

    def test_guarded_scaling_removes_amplitude_without_inflating_small_noise(self):
        spec=preprocessing_config(.1)
        a,b,constant,noise=prepare_expression([[1,2,3],[12,14,16],[7,7,7],[-.01,0,.01]],spec)
        for x,y in zip(a,b):self.assertAlmostEqual(x,y)
        self.assertAlmostEqual(pstdev(a),1)
        self.assertEqual(constant,[0,0,0])
        self.assertLess(pstdev(noise),.1)
        for x,y in zip(noise,[-.1,0,.1]):self.assertAlmostEqual(x,y)
        raw=[[1,2,3],[12,14,16]]
        self.assertEqual(prepare_expression(raw,preprocessing_config(0)),raw)
        self.assertEqual(prepare_expression([raw[1]],spec)[0],b)
        for floor in (-1,float('inf'),float('nan')):
            with self.assertRaises(ValueError):preprocessing_config(floor)
        with self.assertRaises(ValueError):prepare_expression(raw,{'method':'bad'})
        with self.assertRaises(ValueError):prepare_expression(raw,{'method':'none','version':99})

    def test_pair_and_matrix_use_identical_normalized_values(self):
        script=Path(__file__).resolve().parents[1]/'run_dtw.py'
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=root/'input.csv'
            source.write_text('gene_id,t0,t1,t2\na,1,2,4\nb,2,8,3\nc,0,0,0\n')
            for name,options in [('matrix',[]),('pair',['--pair','a','b'])]:
                result=subprocess.run([sys.executable,'-B',str(script),'--input',str(source),
                    '--output',str(root/name),'--open-penalty','.2','--extend-penalty','.05',*options],
                    capture_output=True,text=True)
                self.assertEqual(result.returncode,0,result.stderr)
            with (root/'matrix/distance_matrix.csv').open() as f:
                import csv
                rows=list(csv.reader(f))
            pair=json.loads((root/'pair/run.json').read_text())
            self.assertEqual(float(rows[1][2]),pair['distance'])


if __name__=='__main__': unittest.main()
