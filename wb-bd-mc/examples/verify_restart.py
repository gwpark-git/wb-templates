"""Record bounded same-backend checkpoint continuation checks."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from bdmc import Simulation,SimulationSpec,_core
from bdmc.io import save_checkpoint,load_checkpoint
from bdmc.storage import DataBudget


parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,required=True)
parser.add_argument('--data-budget',type=Path,required=True)
args=parser.parse_args()
budget=DataBudget(args.data_budget)
budget.require_output(args.output,planned_bytes=1_000_000)
args.output.mkdir(parents=True,exist_ok=False)
records=[]
for model in ('brownian','dumbbell'):
    for noise in ('legacy_uniform','gaussian'):
        n=3 if model=='brownian' else 4
        simulation=Simulation(SimulationSpec(model=model,n_particles=n,dimension=1,
                            box_lengths=[10.],dt=.01,seed=317,noise_policy=noise))
        simulation.advance(3)
        name=f'{model}-{noise}'
        save_checkpoint(args.output/name,simulation,budget=budget,lineage_id=name)
        resumed=load_checkpoint(args.output/name)
        reference=simulation.sample(13,stride=4)
        continued=resumed.sample(13,stride=4)
        digest=hashlib.sha256()
        for key in ('positions','unwrapped_positions','forces','steps','times','potential_energy'):
            np.testing.assert_array_equal(reference[key],continued[key])
            digest.update(reference[key].tobytes())
        save_checkpoint(args.output/f'{name}-final',resumed,budget=budget,lineage_id=name)
        records.append({'model':model,'noise':noise,'split_step':3,'additional_steps':13,
                        'array_equality':'exact','continued_batch_sha256':digest.hexdigest()})
report={'schema':'bdmc.restart-demo.v1','core_sha256':hashlib.sha256(Path(_core.__file__).read_bytes()).hexdigest(),
        'cases':records,'limits':'same backend/version; no legacy MC restart or cross-build guarantee',
        'storage_after':budget.check()}
with budget.writer(args.output/'report.json',planned_bytes=65_536):
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(records,indent=2))
