"""A common command line for typed simulation, learning and saved inference."""
import argparse
from contextlib import nullcontext
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import uuid

import numpy as np


def main(argv=None):
    parser=argparse.ArgumentParser(prog='bdmc')
    commands=parser.add_subparsers(dest='command',required=True)
    run=commands.add_parser('run',help='run the direct Brownian or fixed-dumbbell core')
    run.add_argument('--spec',type=Path,required=True,help='JSON SimulationSpec fields')
    run.add_argument('--initial-positions',type=Path)
    resume=commands.add_parser('resume',help='continue an exact default-core checkpoint')
    resume.add_argument('--checkpoint',type=Path,required=True)
    for command in (run,resume):
        command.add_argument('--steps',type=int,required=True)
        command.add_argument('--stride',type=int,default=1)
        command.add_argument('--output',type=Path,required=True)
        command.add_argument('--checkpoint-output',type=Path)
        command.add_argument('--data-budget',type=Path,help='fixed additional-data baseline; required for unattended work')
    train=commands.add_parser('train',help='train a registered regression method')
    train.add_argument('--dataset',type=Path,required=True)
    train.add_argument('--method',choices=['ridge','gaussian_process','mlp'],default='ridge')
    train.add_argument('--group-by',nargs='+',default=['groups'])
    train.add_argument('--seed',type=int,default=0)
    train.add_argument('--config',type=Path,help='JSON method configuration')
    train.add_argument('--output',type=Path,required=True)
    predict=commands.add_parser('predict',help='infer from an existing artifact')
    predict.add_argument('--artifact',type=Path,required=True)
    predict.add_argument('--features',type=Path,required=True,help='numeric NPY with columns in artifact order')
    predict.add_argument('--output',type=Path,required=True)
    args=parser.parse_args(argv)
    if args.command in ('run','resume'):
        from . import Simulation,SimulationSpec,_core
        if args.steps<0 or args.stride<=0:
            raise ValueError('steps must be nonnegative and stride positive')
        if args.output.exists():
            raise FileExistsError(args.output)
        if args.checkpoint_output is not None:
            cp=args.checkpoint_output.resolve();out=args.output.resolve()
            if cp==out or cp in out.parents:
                raise ValueError('checkpoint output must not equal or contain the trajectory directory')
            if cp.exists():
                raise FileExistsError(cp)
        if args.command=='run':
            spec=SimulationSpec(**json.loads(args.spec.read_text()))
            initial=None if args.initial_positions is None else np.load(args.initial_positions,allow_pickle=False)
            simulation=Simulation(spec,initial_positions=initial)
            lineage_id=str(uuid.uuid4())
        else:
            from .io import load_checkpoint
            simulation=load_checkpoint(args.checkpoint)
            spec=simulation.spec
            lineage_id=json.loads((args.checkpoint/'checkpoint.json').read_text())['lineage_id']
        record={'schema':'bdmc.simulation-run.v1','created_at':datetime.now(timezone.utc).isoformat(),
                'spec':spec.to_dict(),'metadata':dict(simulation.metadata),'steps':args.steps,'stride':args.stride,
                'core_sha256':hashlib.sha256(Path(_core.__file__).read_bytes()).hexdigest(),
                'trajectory_sha256':'0'*64,'run_id':str(uuid.uuid4()),'lineage_id':lineage_id,
                'recording':'initial and post-completed-step states; full final partial stride included'}
        if args.command=='resume':
            record['resumed_from_checkpoint']=str(args.checkpoint.resolve())
        budget=None
        if args.data_budget is not None:
            from .storage import DataBudget
            budget=DataBudget(args.data_budget)
            frames=1+(args.steps+args.stride-1)//args.stride
            raw=frames*(3*spec.n_particles*spec.dimension*8+24)+spec.n_particles*64
            planned=raw+raw//1000+131_072+len((json.dumps(record,indent=2,allow_nan=False)+'\n').encode())
            if args.checkpoint_output is not None:
                planned+=spec.n_particles*(2*spec.dimension*8+40)+131_072
                planned+=len(json.dumps({'lineage_id':lineage_id},indent=2).encode())
                budget.require_output(args.checkpoint_output)
            budget.require_output(args.output,planned_bytes=planned)
        batch=simulation.sample(args.steps,stride=args.stride,include_initial=True)
        with budget.writer(args.output,planned_bytes=planned) if budget is not None else nullcontext():
            args.output.mkdir(parents=True,exist_ok=False)
            arrays={key:value for key,value in batch.items() if isinstance(value,np.ndarray)}
            np.savez_compressed(args.output/'trajectory.npz',**arrays)
            record['trajectory_sha256']=hashlib.sha256((args.output/'trajectory.npz').read_bytes()).hexdigest()
            (args.output/'run.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
            if args.checkpoint_output is not None:
                from .io import save_checkpoint
                save_checkpoint(args.checkpoint_output,simulation,budget=budget,lineage_id=lineage_id)
    elif args.command=='train':
        from .io import load_dataset
        from .ml import grouped_split,train as fit
        dataset=load_dataset(args.dataset)
        split=grouped_split(dataset,seed=args.seed,group_by=tuple(args.group_by))
        config={} if args.config is None else json.loads(args.config.read_text())
        artifact=fit(dataset,method=args.method,split=split,config=config,seed=args.seed)
        artifact.save(args.output)
    else:
        from .ml import load_artifact
        if args.output.exists():
            raise FileExistsError(args.output)
        artifact=load_artifact(args.artifact)
        values=artifact.predict(np.load(args.features,allow_pickle=False))
        args.output.parent.mkdir(parents=True,exist_ok=True)
        # Passing a handle avoids NumPy silently adding another filename suffix.
        with args.output.open('xb') as stream:
            np.save(stream,values,allow_pickle=False)
    print(str(args.output.resolve()))
    return 0
