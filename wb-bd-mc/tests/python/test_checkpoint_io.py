import json
import subprocess
import sys

import numpy as np
import pytest

from bdmc import Simulation,SimulationSpec
from bdmc.io import load_checkpoint,read_checkpoint,save_checkpoint
from bdmc.storage import DataBudget,DataBudgetExceeded,create_data_baseline


@pytest.mark.parametrize('noise',['legacy_uniform','gaussian'])
@pytest.mark.parametrize('model',['brownian','dumbbell'])
def test_checkpoint_file_exact_continuation(tmp_path,noise,model):
    # Brownian uses odd draw counts to exercise cached Gaussian variates.
    n=3 if model=='brownian' else 4
    sim=Simulation(SimulationSpec(model=model,n_particles=n,dimension=1,
        box_lengths=[10.],noise_policy=noise,seed=97,dt=.01))
    sim.advance(3)
    checkpoint=tmp_path/'checkpoint'
    save_checkpoint(checkpoint,sim)
    restored=load_checkpoint(checkpoint)
    expected=sim.sample(13,stride=4)
    resumed=restored.sample(13,stride=4)
    for key in ('positions','unwrapped_positions','forces','steps','times','potential_energy'):
        np.testing.assert_array_equal(resumed[key],expected[key])
    if noise=='gaussian' and model=='brownian':
        assert read_checkpoint(checkpoint)['rng']['has_spare']
    with pytest.raises(FileExistsError):
        save_checkpoint(checkpoint,sim)


def test_checkpoint_numeric_corruption_and_version_rejected(tmp_path):
    sim=Simulation(SimulationSpec())
    save_checkpoint(tmp_path/'a',sim)
    arrays=tmp_path/'a/state.npz'
    arrays.write_bytes(arrays.read_bytes()+b'corrupt')
    with pytest.raises(ValueError,match='checksum'):
        load_checkpoint(tmp_path/'a')
    save_checkpoint(tmp_path/'b',sim)
    manifest=tmp_path/'b/checkpoint.json'
    data=json.loads(manifest.read_text());data['snapshot']['version']=999
    manifest.write_text(json.dumps(data))
    with pytest.raises(ValueError):
        load_checkpoint(tmp_path/'b')


def test_checkpoint_budget_checked_before_output_creation(tmp_path):
    root=tmp_path/'runs';root.mkdir()
    baseline=tmp_path/'budget.json';create_data_baseline([root],baseline)
    budget=DataBudget(baseline)
    sim=Simulation(SimulationSpec())
    with pytest.raises(ValueError,match='outside'):
        save_checkpoint(tmp_path/'outside',sim,budget=budget)
    assert not (tmp_path/'outside').exists()
    save_checkpoint(root/'good',sim,budget=budget)
    assert budget.check()['additional_bytes']>0


def test_cli_resume_preserves_trajectory_and_lineage(tmp_path):
    root=tmp_path/'runs';root.mkdir()
    baseline=tmp_path/'budget.json';create_data_baseline([root],baseline)
    specfile=tmp_path/'spec.json'
    spec=dict(n_particles=3,dimension=1,box_lengths=[10.],noise_policy='gaussian',seed=55,dt=.01)
    specfile.write_text(json.dumps(spec))
    def call(*args):
        subprocess.run([sys.executable,'-m','bdmc',*map(str,args)],check=True,capture_output=True,timeout=30)
    call('run','--spec',specfile,'--steps',3,'--output',root/'first',
         '--checkpoint-output',root/'checkpoint','--data-budget',baseline)
    call('resume','--checkpoint',root/'checkpoint','--steps',7,'--output',root/'resumed',
         '--checkpoint-output',root/'checkpoint2','--data-budget',baseline)
    expected=Simulation(SimulationSpec(**spec));expected.advance(3)
    reference=expected.sample(7,stride=1)
    with np.load(root/'resumed/trajectory.npz',allow_pickle=False) as arrays:
        np.testing.assert_array_equal(arrays['unwrapped_positions'],reference['unwrapped_positions'])
        np.testing.assert_array_equal(arrays['steps'],reference['steps'])
    first=json.loads((root/'first/run.json').read_text())
    second=json.loads((root/'resumed/run.json').read_text())
    checkpoint=json.loads((root/'checkpoint2/checkpoint.json').read_text())
    assert first['run_id']!=second['run_id']
    assert first['lineage_id']==second['lineage_id']==checkpoint['lineage_id']


def test_checkpoint_lineage_json_bytes_are_budgeted(tmp_path):
    root=tmp_path/'runs';root.mkdir()
    baseline=tmp_path/'budget.json';create_data_baseline([root],baseline)
    config=json.loads(baseline.read_text());config.update(soft_limit_bytes=250_000,additional_limit_bytes=300_000)
    baseline.write_text(json.dumps(config))
    with pytest.raises(DataBudgetExceeded):
        save_checkpoint(root/'large-lineage',Simulation(SimulationSpec()),
                        lineage_id='가'*50_000,budget=DataBudget(baseline))
    assert not (root/'large-lineage').exists()


def test_cli_resume_budgets_escaped_lineage_in_raw_run_record(tmp_path):
    source=tmp_path/'existing-checkpoint'
    save_checkpoint(source,Simulation(SimulationSpec(n_particles=3,dimension=1)),lineage_id='가'*50_000)
    root=tmp_path/'runs';root.mkdir()
    baseline=tmp_path/'budget.json';create_data_baseline([root],baseline)
    config=json.loads(baseline.read_text());config.update(soft_limit_bytes=200_000,additional_limit_bytes=250_000)
    baseline.write_text(json.dumps(config))
    process=subprocess.run([sys.executable,'-m','bdmc','resume','--checkpoint',str(source),
        '--steps','1','--output',str(root/'new'),'--data-budget',str(baseline)],
        capture_output=True,text=True,timeout=30)
    assert process.returncode!=0 and 'DataBudgetExceeded' in process.stderr
    assert not (root/'new').exists()
