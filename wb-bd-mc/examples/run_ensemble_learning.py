"""Bounded local-parallel generation through grouped learning and saved inference."""
import argparse
import json
from pathlib import Path


def main():
    from bdmc import SimulationSpec
    from bdmc.ensemble import EnsembleJob,execute_ensemble
    from bdmc.ensemble_data import ensemble_dataset
    from bdmc.io import save_dataset
    from bdmc.ml import Dataset,RegressionTask,grouped_split,train,load_artifact
    from bdmc.storage import DataBudget
    import numpy as np

    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--data-budget',type=Path,required=True)
    args=parser.parse_args()
    budget=DataBudget(args.data_budget)
    budget.require_output(args.output,planned_bytes=2_000_000)
    args.output.mkdir(parents=True,exist_ok=False)
    jobs=[]
    for c,D in enumerate(np.linspace(.4,1.6,6)):
        for s in range(2):
            run_id=f'condition{c}-seed{s}'
            jobs.append(EnsembleJob(run_id=run_id,
                spec=SimulationSpec(n_particles=256,dimension=3,diffusivity=float(D),seed=95000+c*10+s,dt=.004),
                steps=40,stride=10,observables=('msd',),condition_id=f'condition{c}',lineage_id=run_id))
    report=execute_ensemble(jobs,workers=2,timeout_seconds=60,
                            output_directory=args.output/'ensemble',budget=budget)
    original=ensemble_dataset(report,feature_names=('diffusivity','time'),target_names=('msd',),
        task=RegressionTask(name='parallel_brownian_calibration',target_kind='realization',target_units=('length_squared',)))
    # Declare the known Brownian Dt feature before model fitting/evaluation.
    data=Dataset(np.column_stack([original.X,original.X[:,0]*original.X[:,1]]),original.y,original.groups,
         feature_names=(*original.feature_names,'diffusivity_times_time'),target_names=original.target_names,
         group_keys=original.group_keys,task=original.task,
         metadata={**original.metadata,'declared_feature_engineering':'diffusivity*time'})
    save_dataset(args.output/'dataset',data,budget=budget)
    split=grouped_split(data,seed=91,group_by=('condition','continuation'))
    artifact=train(data,method='ridge',split=split,config={'alpha':1e-6},seed=91)
    with budget.writer(args.output/'ridge',planned_bytes=65_536):
        artifact.save(args.output/'ridge')
    expected=artifact.predict(data.X[split.test])
    restored=load_artifact(args.output/'ridge').predict(data.X[split.test])
    np.testing.assert_array_equal(expected,restored)
    theory=6*data.X[split.test,2:3]
    result={'purpose':'bounded execution/learning integration, research target remains open',
        'runs':len(jobs),'rows':len(data.X),'workers':2,'raw_trajectories_saved':False,
        'metrics':artifact.metrics,'analytic_test_rmse':float(np.sqrt(np.mean((data.y[split.test]-theory)**2))),
        'saved_inference_match':True,'grouping':split.group_by,
        'storage_after':budget.check(),'limits':'local CPU ensemble only; no single-trajectory MPI/GPU acceleration'}
    text=json.dumps(result,indent=2)+'\n'
    with budget.writer(args.output/'report.json',planned_bytes=len(text.encode())+4096):
        (args.output/'report.json').write_text(text)
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
