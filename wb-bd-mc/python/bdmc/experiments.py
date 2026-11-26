"""Small, recorded calibration experiments; these do not select the research task."""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np

from .ml import Dataset, RegressionTask, compare, grouped_split, load_artifact


def dumbbell_second_moment(*, time_steps: int, dt: float, dimension: int,
                           mobility: float, diffusivity: float, spring_alpha: float,
                           initial_q_squared: float = 1.) -> float:
    """Exact finite-Euler second moment before minimum-image branch crossings.

    q is partner minus particle. Independent bead noise has variance 2Ddt;
    relative noise has variance 4Ddt per coordinate for either supported policy.
    """
    if (isinstance(time_steps, (bool,np.bool_)) or not isinstance(time_steps, (int, np.integer))
            or not 0 <= time_steps <= 2**64-1):
        raise ValueError('time_steps must be a nonnegative uint64 integer')
    if not np.isfinite([dt,mobility,diffusivity,spring_alpha,initial_q_squared]).all() or initial_q_squared < 0:
        raise ValueError('finite parameters and nonnegative initial second moment required')
    if dt <= 0 or dimension not in (1,2,3) or mobility < 0 or diffusivity < 0 or spring_alpha <= 0:
        raise ValueError('invalid finite-step moment parameters')
    kappa = float(dimension)*float(spring_alpha)*float(spring_alpha)
    drift = 2*float(mobility)*float(dt)*kappa
    noise_variance = 4*dimension*(float(diffusivity)*float(dt))
    if not all(math.isfinite(x) for x in (kappa,drift,noise_variance)):
        raise ValueError('derived moment scales exceed finite floating-point range')
    if drift >= 2:
        raise ValueError('unstable Euler pair relaxation')
    a2 = (1-drift)**2
    # Stable O(1) geometric sum, including the force-free limit and long horizons.
    if a2 == 1:
        geometric = float(time_steps)
    elif a2 == 0:
        geometric = float(time_steps > 0)
    else:
        geometric = -math.expm1(int(time_steps)*math.log(a2))/(1-a2)
    result = float(a2**int(time_steps)*initial_q_squared + noise_variance*geometric)
    if not math.isfinite(result):
        raise ValueError('second moment exceeds finite floating-point range')
    return result


def calibration_dataset(model: str, *, particles: int = 512, steps: int = 100,
                        stride: int = 20, seeds_per_condition: int = 3,
                        conditions: int = 9) -> tuple[Dataset, np.ndarray, dict]:
    """Generate actual modern-core runs and analytic moment references.

    Whole conditions are held together later. Repeated times from one run never
    cross partitions. Independent runs use distinct seeds, which are metadata.
    """
    from . import Simulation, SimulationSpec, _core
    if model not in ('brownian','dumbbell'):
        raise ValueError('calibration supports brownian and dumbbell only')
    if particles < 2 or particles % 2 or steps <= 0 or stride <= 0 or seeds_per_condition < 1 or conditions < 3:
        raise ValueError('positive even particle count, steps, stride and >=3 conditions required')
    X, y, reference, groups, condition_keys, run_records = [], [], [], [], [], []
    for c, parameter in enumerate(np.linspace(.3,1.8,conditions)):
        for s in range(seeds_per_condition):
            seed = 70000+c*100+s
            D = float(parameter) if model=='brownian' else 1.
            alpha = float(parameter) if model=='dumbbell' else 1.
            spec = SimulationSpec(model=model,n_particles=particles,dimension=3,
                                  box_lengths=[10000.]*3,dt=.004,mobility=1.,
                                  diffusivity=D,spring_alpha=alpha,seed=seed,
                                  noise_policy='legacy_uniform')
            initial = np.zeros((particles,3))
            if model=='dumbbell':
                initial[particles//2:,0] = 1.
            sim = Simulation(spec,initial_positions=initial)
            batch = sim.sample(steps,stride=stride,include_initial=True)
            run_id = f'{model}-condition{c}-seed{seed}'
            run_records.append({'run_id':run_id,'spec':spec.to_dict(),'metadata':dict(sim.metadata)})
            for frame,(t,n) in enumerate(zip(batch['times'],batch['steps'])):
                if n==0:
                    continue
                positions = batch['unwrapped_positions'][frame]
                if model=='brownian':
                    observed = np.square(positions-initial).sum(axis=1).mean()
                    exact = 2*3*D*float(t)
                    # The Dt feature is declared in advance from Brownian theory.
                    features = [D,float(t),D*float(t)]
                else:
                    q = positions[particles//2:]-positions[:particles//2]
                    observed = np.square(q).sum(axis=1).mean()
                    exact = dumbbell_second_moment(time_steps=int(n),dt=.004,dimension=3,
                                                   mobility=1.,diffusivity=D,spring_alpha=alpha)
                    features = [alpha,float(t)]
                X.append(features); y.append(float(observed)); reference.append(exact)
                groups.append(run_id); condition_keys.append(f'{model}-condition{c}')
    metadata = {'schema':'bdmc.calibration.v1','model':model,
                'purpose':'engineering/learning-method calibration; research target remains undecided',
                'label':'one finite-ensemble realization of a moment, not a noiseless expectation',
                'units':'normalised length and time; D and mobility explicit',
                'noise_policy':'legacy_uniform','boundary':'periodic box10000; q theory assumes same image branch',
                'feature_engineering':'Brownian includes D*time; dumbbell uses alpha and time',
                'core_binary_sha256':hashlib.sha256(Path(_core.__file__).read_bytes()).hexdigest(),
                'runs':run_records}
    names = ('diffusivity','time','diffusivity_times_time') if model=='brownian' else ('spring_alpha','time')
    dataset = Dataset(np.asarray(X),np.asarray(y),np.asarray(groups),feature_names=names,
                      target_names=('msd' if model=='brownian' else 'connector_second_moment',),
                      group_keys={'condition':np.asarray(condition_keys)},metadata=metadata,
                      task=RegressionTask(name=f'{model}_calibration',target_kind='realization',target_units=('length_squared',)))
    return dataset,np.asarray(reference)[:,None],{'n_runs':len(run_records),'n_rows':len(X)}


def train_calibration(output: str | Path, *, model: str, include_torch: bool = True,
                      **generation_options) -> dict:
    """Generate -> condition split -> fit fixed baselines -> save/reload/infer."""
    output = Path(output)
    output.mkdir(parents=True,exist_ok=False)
    start = time.perf_counter()
    dataset,theory,counts = calibration_dataset(model,**generation_options)
    from .io import save_dataset
    save_dataset(output/'training_data',dataset)
    generation_seconds = time.perf_counter()-start
    split = grouped_split(dataset,seed=19,group_by=('groups','condition'))
    methods = {'ridge':{'alpha':1e-6},'gaussian_process':{'noise_variance':.01}}
    if include_torch:
        methods['mlp'] = {'hidden_layers':[16,16],'epochs':200,'learning_rate':.01}
    # Keep data and metadata inspectable independently of the trained artifacts.
    np.savez_compressed(output/'dataset.npz',X=dataset.X,y=dataset.y,groups=dataset.groups,
                        condition=dataset.group_keys['condition'],analytic_reference=theory,
                        train=split.train,validation=split.validation,test=split.test)
    (output/'dataset.json').write_text(json.dumps({'feature_names':dataset.feature_names,
        'target_names':dataset.target_names,'metadata':dataset.metadata,'fingerprint':dataset.fingerprint},indent=2)+'\n')
    start = time.perf_counter()
    artifacts = compare(dataset,methods,split=split,seed=23)
    fit_seconds = time.perf_counter()-start
    records = {}
    for method,artifact in artifacts.items():
        artifact.save(output/method)
        reloaded = load_artifact(output/method)
        prediction = artifact.predict(dataset.X[split.test])
        np.testing.assert_allclose(reloaded.predict(dataset.X[split.test]),prediction,rtol=1e-12,atol=1e-12)
        np.save(output/f'{method}_test_predictions.npy',prediction,allow_pickle=False)
        records[method] = {'metrics':artifact.metrics,'save_load_prediction_match':True}
    truth = dataset.y[split.test]
    mean_baseline = np.full_like(truth,dataset.y[split.train].mean())
    report = {'created_at':datetime.now(timezone.utc).isoformat(),'model':model,**counts,
              'dataset_fingerprint':dataset.fingerprint,'split':'whole condition and source run',
              'split_rows':{name:len(getattr(split,name)) for name in ('train','validation','test')},
              'generation_seconds':generation_seconds,'fit_seconds_all_methods':fit_seconds,
              'analytic_test_rmse':float(np.sqrt(np.mean((truth-theory[split.test])**2))),
              'train_mean_test_rmse':float(np.sqrt(np.mean((truth-mean_baseline)**2))),
              'methods':records,'limitations':['small calibration only; fixed configurations, no hyperparameter search',
              'no research target selected; no scientific claim of ML outperforming known theory',
              'analytic references are finite-step moment laws, not whole-distribution or legacy trajectory equivalence',
              'model uncertainty is not calibrated physical sampling error; no GPU/MPI backend']}
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


def train_network_demo(output: str | Path, *, executable: str | Path, template: str | Path,
                       include_torch: bool = True) -> dict:
    """Compatibility-path demo, not validation/migration of the MC kinetics.

    The ten uncorrected source-audit candidates remain visible in run metadata.
    Only one topology thread, default friction and no shear are accepted.
    """
    from .legacy import LegacyRunner
    output = Path(output)
    output.mkdir(parents=True,exist_ok=False)
    runner = LegacyRunner(executable,template)
    X,y,groups,condition_keys,records = [],[],[],[],[]
    start = time.perf_counter()
    for c,coefficient in enumerate(np.linspace(10,40,9)):
        for s in range(3):
            seed=81000+c*100+s
            result=runner.run(output/'runs',{'Np':64,'Nt':101,'dt/tauR':.001,'Rt':.001,
                        'box_dimension':6,'LBx':6,'LBy':6,'LBz':6,'N_THREADS_BD':1,
                        'N_skip_ener':10,'N_skip_file':10,'record_RDIST':False,
                        'repulsion_coefficient':float(coefficient),
                        'basic_random_seed':seed,'basic_random_seed_SS':seed+37})
            run_id=result.directory.name
            # Verify reported bridge count against independent labelled chains.
            chain=np.loadtxt(result.directory/'case.chain',ndmin=2).astype(int)
            n_chains=64*int(result.parameters['N_chains_per_particle'])
            counted=(chain[:,:n_chains]!=chain[:,n_chains:]).sum(axis=1)
            np.testing.assert_array_equal(counted,result.bridge_counts)
            for frame in range(1,len(result.times)):
                X.append([float(coefficient),float(result.times[frame])])
                y.append(float(result.bridge_counts[frame])/n_chains)
                groups.append(run_id);condition_keys.append(f'coefficient{c}')
            records.append({'run_id':run_id,'record':str(result.directory/'run.json'),
                            'executable_sha256':result.metadata['executable_sha256']})
    metadata={'schema':'bdmc.network-demo.v1','purpose':'engineering interoperability demo only',
              'backend':'legacy_subprocess; original topology-MC schedule preserved',
              'label':'recorded bridge fraction, independently checked against chain endpoints',
              'scientific_status':'legacy source-audit issues remain; not a physically validated training dataset',
              'historical_regression':'deferred pending SSD data','runs':records}
    dataset=Dataset(np.asarray(X),np.asarray(y),np.asarray(groups),
            feature_names=('repulsion_coefficient','legacy_tauR_time'),target_names=('bridge_fraction',),
            group_keys={'condition':np.asarray(condition_keys)},metadata=metadata,
            task=RegressionTask(name='network_adapter_demo',target_kind='realization',target_units=('dimensionless',)))
    split=grouped_split(dataset,seed=19,group_by=('groups','condition'))
    from .io import save_dataset
    save_dataset(output/'training_data',dataset)
    np.savez_compressed(output/'dataset.npz',X=dataset.X,y=dataset.y,groups=dataset.groups,
                        condition=dataset.group_keys['condition'],train=split.train,
                        validation=split.validation,test=split.test)
    (output/'dataset.json').write_text(json.dumps({'metadata':metadata,'fingerprint':dataset.fingerprint,
        'feature_names':dataset.feature_names,'target_names':dataset.target_names},indent=2)+'\n')
    methods={'ridge':{'alpha':1e-6},'gaussian_process':{'noise_variance':.05}}
    if include_torch:
        methods['mlp']={'hidden_layers':[16,16],'epochs':200}
    artifacts=compare(dataset,methods,split=split,seed=23)
    summary={}
    for method,artifact in artifacts.items():
        artifact.save(output/method)
        original=artifact.predict(dataset.X[split.test])
        restored=load_artifact(output/method).predict(dataset.X[split.test])
        np.testing.assert_allclose(original,restored,rtol=1e-12,atol=1e-12)
        summary[method]={'metrics':artifact.metrics,'save_load_prediction_match':True}
    truth=dataset.y[split.test]
    report={'model':'transient_network_legacy_demo','n_runs':len(records),'n_rows':len(y),
            'elapsed_seconds':time.perf_counter()-start,'dataset_fingerprint':dataset.fingerprint,
            'methods':summary,'train_mean_test_rmse':float(np.sqrt(np.mean((truth-dataset.y[split.train].mean())**2))),
            'scope':metadata['scientific_status'],'split':'whole condition and source run',
            'split_rows':{name:len(getattr(split,name)) for name in ('train','validation','test')}}
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    return report
