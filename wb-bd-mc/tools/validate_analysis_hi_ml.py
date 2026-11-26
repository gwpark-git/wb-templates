#!/usr/bin/env python3
"""Executed observer/HI -> grouped Dataset -> saved/reloaded ML smoke.

This is interface verification on short nonequilibrium paths, not adoption of
a scientific training target or evidence of learned transport/rheology.
"""
from datetime import datetime, timezone
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time

CODE=Path(__file__).resolve().parents[1]
# A new labelled engineering ensemble, prospectively fixed; never retry/reseed.
SEEDS=tuple(range(810729,810741))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--package',type=Path,required=True)
    p.add_argument('--run-dir',type=Path,default=CODE/'runs/repair-analysis-hi-20261111/ml-integration')
    p.add_argument('--budget',type=Path,default=CODE/'DATA_BUDGET.json')
    args=p.parse_args()
    package=args.package.resolve(strict=True)
    # Use the validator's staged sitecustomize in this and spawned processes.
    sys.path.insert(0,str(package.parent))
    import importlib.util
    spec=importlib.util.spec_from_file_location('_integration_site',package.parent/'sitecustomize.py')
    finder=importlib.util.module_from_spec(spec);spec.loader.exec_module(finder)
    import numpy as np
    import bdmc
    from bdmc import _analysis,_hydrodynamics,_reference
    from bdmc.analysis import pair_correlation,directional_pair_correlation,coordination_number
    from bdmc.reference import RepulsiveSpec,RepulsiveSimulation
    from bdmc.hydrodynamics import HydrodynamicSpec,HydrodynamicSimulation
    from bdmc.ml import Dataset,RegressionTask,train,load_artifact
    from bdmc.io import save_dataset,load_dataset
    from bdmc.storage import DataBudget
    if Path(bdmc.__file__).resolve().parent!=package: raise RuntimeError('stale package')
    binary={name:hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest()
            for name,module in (('_analysis',_analysis),('_hydrodynamics',_hydrodynamics),('_reference',_reference))}
    budget=DataBudget(args.budget)
    with budget.writer(args.run_dir,8_000_000):
        args.run_dir.mkdir(parents=True,exist_ok=True)
        run=Path(tempfile.mkdtemp(prefix='run-',dir=args.run_dir))
        report=dict(schema='bdmc.analysis-hi-ml-smoke.v1',created_at=datetime.now(timezone.utc).isoformat(),
                    seeds=list(SEEDS),native_sha256=binary,status='running',models={},
                    scope='Engineering interface smoke, finite paths; no scientific target adoption, equilibrated dataset or learned transport claim.')
        (run/'seed_manifest.json').write_text(json.dumps({'schema':'bdmc.engineering-seeds.v1','seeds':SEEDS,
            'policy':'New labelled ensemble; source groups are seed lineages shared across conditions; no retries.'},indent=2)+'\n')
        start=time.perf_counter()
        try:
            for model in ('periodic_structure','free_space_hi'):
                X=[];y=[];groups=[];conditions=[];simulation_seconds=0.;analysis_seconds=0.
                for condition in (.75,1.,1.25):
                    for seed in SEEDS:
                        now=time.perf_counter()
                        if model=='periodic_structure':
                            q=np.random.default_rng(seed).uniform(0,6,(32,3))
                            sim=RepulsiveSimulation(RepulsiveSpec(n_particles=32,box_length=6,
                                epsilon=condition,sigma=1.25,dt=1e-4,diffusion=.1,seed=seed),q)
                            batch=sim.sample(30,5,True)
                        else:
                            q=np.array([[0.,0.,0.],[2.,0.,0.],[0.,2.,0.],[0.,0.,2.]])
                            forces=np.zeros((4,3));forces[:,0]=.1
                            sim=HydrodynamicSimulation(HydrodynamicSpec(n_particles=4,radius=.2,
                                viscosity=condition/(6*np.pi*.2),dt=1e-4,thermal_energy=.01,
                                external_forces=forces,seed=seed),q)
                            batch=sim.sample(30,5,True)
                        simulation_seconds+=time.perf_counter()-now
                        now=time.perf_counter()
                        if model=='periodic_structure':
                            rdf=pair_correlation(batch['positions'],[6.]*3,np.linspace(0,1.25,11))
                            direct=coordination_number(rdf);integrated=coordination_number(rdf,source='g')
                            np.testing.assert_allclose(direct,integrated,rtol=1e-14,atol=1e-14)
                            sectors=directional_pair_correlation(batch['positions'],[6.]*3,
                                np.linspace(0,1.25,11),[-1.,0.,1.],[0.,np.pi,2*np.pi])
                            np.testing.assert_array_equal(sectors['counts'].sum(axis=(1,2)),rdf['counts'])
                            target=float(direct[-1])
                        else:
                            # No periodic RDF is applied to this free-space ensemble.
                            target=float(np.mean(batch['positions'][-1,:,0]-q[:,0]))
                        analysis_seconds+=time.perf_counter()-now
                        X.append([condition]);y.append([target]);groups.append(str(seed));conditions.append(str(condition))
                dataset=Dataset(X,y,groups,feature_names=('epsilon' if model=='periodic_structure' else 'viscosity_scale',),
                    target_names=('coordination_r_le_1p25' if model=='periodic_structure' else 'finite_path_mean_dx',),
                    group_keys={'condition':conditions},task=RegressionTask(name=model,target_kind='realization'),
                    metadata={'source':'native short paths','native_sha256':binary,'seeds':list(SEEDS),
                              'model':model,'clock':{'dt':1e-4,'steps':30,'stride':5},
                              'limits':'Not equilibrated; grouped by seed lineage; no claim of condition extrapolation.'})
                now=time.perf_counter();save_dataset(run/model/'dataset',dataset,budget=budget)
                loaded=load_dataset(run/model/'dataset')
                if dataset.fingerprint!=loaded.fingerprint: raise RuntimeError('dataset reload mismatch')
                dataset_io=time.perf_counter()-now
                now=time.perf_counter();artifact=train(loaded,method='ridge',seed=150)
                training=time.perf_counter()-now
                now=time.perf_counter();artifact.save(run/model/'artifact')
                reloaded=load_artifact(run/model/'artifact')
                predicted=artifact.predict(loaded.X)
                np.testing.assert_allclose(reloaded.predict(loaded.X),predicted,rtol=0,atol=0)
                artifact_io=time.perf_counter()-now
                np.savez_compressed(run/model/'inference-input.npz',X=loaded.X,prediction=predicted)
                child=('import pathlib,numpy as np,bdmc; from bdmc.ml import load_artifact; '
                       f'\nif pathlib.Path(bdmc.__file__).resolve().parent!=pathlib.Path({str(package)!r}): '
                       'raise RuntimeError("stale inference package")\n'
                       f'a=load_artifact({str(run/model/"artifact")!r}); '
                       f'd=np.load({str(run/model/"inference-input.npz")!r}); '
                       'np.testing.assert_allclose(a.predict(d["X"]),d["prediction"],rtol=0,atol=0); '
                       'print("PASS separate-process saved inference")')
                now=time.perf_counter()
                process=subprocess.run([sys.executable,'-c',child],capture_output=True,text=True,timeout=30,
                    env={**os.environ,'PYTHONPATH':str(package.parent),'PYTHONDONTWRITEBYTECODE':'1',
                         'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'})
                (run/model/'inference.log').write_text(process.stdout+process.stderr)
                if process.returncode: raise RuntimeError('separate inference failed')
                report['models'][model]=dict(rows=len(X),independent_seed_groups=len(set(groups)),
                    fingerprint=dataset.fingerprint,metrics=artifact.metrics,
                    timings=dict(simulation_seconds=simulation_seconds,analysis_seconds=analysis_seconds,
                                 dataset_io_seconds=dataset_io,training_seconds=training,artifact_io_seconds=artifact_io,
                                 separate_inference_seconds=time.perf_counter()-now),
                    exact_reload=True,separate_process_inference=True)
            report['status']='passed'
        except Exception as exc:
            report.update(status='failed',error=f'{type(exc).__name__}: {exc}');raise
        finally:
            report.update(elapsed_seconds=time.perf_counter()-start,storage_after=budget.check(),
                          completed_at=datetime.now(timezone.utc).isoformat())
            (run/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
            print(run/'report.json',report['status'],flush=True)


if __name__=='__main__': main()
