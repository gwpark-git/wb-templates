import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from bdmc.io import load_dataset,save_dataset
from bdmc.ml import Dataset,RegressionTask,load_artifact
from bdmc.storage import DataBudget,DataBudgetExceeded,create_data_baseline


def dataset():
    x=np.arange(18,dtype=float)[:,None]
    return Dataset(x,2*x+1,np.repeat(np.arange(6),3),
        feature_names=('parameter',),target_names=('response',),
        group_keys={'condition':np.repeat(np.arange(6),3)},
        task=RegressionTask(name='linear_exchange_check',target_kind='ensemble_mean',target_units=('unit',)),
        metadata={'units':'normalised','fixture':'known linear response'})


def cli(*args):
    return subprocess.run([sys.executable,'-m','bdmc',*map(str,args)],check=True,
                          capture_output=True,text=True,timeout=30)


def test_dataset_roundtrip_full_identity(tmp_path):
    original=dataset()
    save_dataset(tmp_path/'dataset',original)
    restored=load_dataset(tmp_path/'dataset')
    assert restored.fingerprint==original.fingerprint
    assert restored.task==original.task
    assert restored.feature_names==original.feature_names
    np.testing.assert_array_equal(restored.group_keys['condition'],original.group_keys['condition'])
    with pytest.raises(FileExistsError):
        save_dataset(tmp_path/'dataset',original)


def test_dataset_corruption_and_semantic_edit_rejected(tmp_path):
    save_dataset(tmp_path/'data',dataset())
    manifest=tmp_path/'data/dataset.json'
    record=json.loads(manifest.read_text());record['task']['target_kind']='realization'
    manifest.write_text(json.dumps(record))
    with pytest.raises(ValueError,match='fingerprint'):
        load_dataset(tmp_path/'data')
    data=tmp_path/'data/data.npz'
    data.write_bytes(data.read_bytes()+b'tamper')
    with pytest.raises(ValueError,match='checksum'):
        load_dataset(tmp_path/'data')


def test_cli_real_core_run(tmp_path):
    spec=tmp_path/'spec.json'
    spec.write_text(json.dumps({'model':'brownian','n_particles':8,'dimension':3,'dt':.01,'seed':91}))
    cli('run','--spec',spec,'--steps',7,'--stride',3,'--output',tmp_path/'run')
    with np.load(tmp_path/'run/trajectory.npz',allow_pickle=False) as result:
        np.testing.assert_array_equal(result['steps'],[0,3,6,7])
        np.testing.assert_allclose(result['times'],[0,.03,.06,.07])
        assert result['positions'].shape==(4,8,3)
    record=json.loads((tmp_path/'run/run.json').read_text())
    assert record['spec']['seed']==91 and record['core_sha256']


def test_cli_training_and_prediction(tmp_path):
    original=dataset()
    save_dataset(tmp_path/'data',original)
    cli('train','--dataset',tmp_path/'data','--group-by','condition',
        '--method','ridge','--output',tmp_path/'artifact')
    features=tmp_path/'x.npy';np.save(features,[[4.],[8.]])
    output=tmp_path/'prediction.array'
    cli('predict','--artifact',tmp_path/'artifact','--features',features,'--output',output)
    np.testing.assert_allclose(np.load(output,allow_pickle=False),[[9.],[17.]],rtol=1e-5)
    artifact=load_artifact(tmp_path/'artifact')
    assert 'groups' in artifact.split.group_by
    assert artifact.provenance['dataset_fingerprint']==original.fingerprint


def test_complete_manifest_is_included_in_data_budget(tmp_path):
    root=tmp_path/'runs';root.mkdir()
    baseline=tmp_path/'budget.json';create_data_baseline([root],baseline)
    config=json.loads(baseline.read_text());config.update(soft_limit_bytes=100_000,additional_limit_bytes=150_000)
    baseline.write_text(json.dumps(config))
    data=Dataset([[0.],[1.],[2.]],[0.,1.,2.],['a','b','c'],feature_names=('x'*200_000,))
    with pytest.raises(DataBudgetExceeded):
        save_dataset(root/'oversized-metadata',data,budget=DataBudget(baseline))
    assert not (root/'oversized-metadata').exists()
