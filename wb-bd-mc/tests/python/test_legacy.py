import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pytest

from bdmc.legacy import LegacyRunner, read_input

CODE = Path(__file__).resolve().parents[2]
EXE = Path(os.environ.get('BDMC_LEGACY_EXECUTABLE',
                          str(CODE/'build/g++/stochastic_simulation')))


def test_input_duplicate_rejected(tmp_path):
    p = tmp_path/'input.inp'
    p.write_text('Np=4\nNp=6\n')
    with pytest.raises(ValueError,match='duplicate'):
        read_input(p)


@pytest.fixture
def runner():
    if not EXE.exists():
        pytest.skip('build legacy executable with make before integration tests')
    return LegacyRunner(EXE, CODE/'short_test.inp')


@pytest.mark.parametrize('override',[
    {'Method':'UNDEFINED'}, {'friction_junction':'LOOP_DISSOCIATION_TIME'},
    {'Integrator':'RK4'}, {'Step':'not_a_step'},
    {'N_THREADS_SS':2}, {'MECHANICAL_PERTURBATION':'SIMPLE_SHEAR'},
    {'CONTINUATION_TRAJ':True}, {'bogus_key':3}, {'Np':0},
    {'Method':'DUMBBELL','Np':3,'cell_list':False},
    {'filename_base':'elsewhere'}, {'Method':'BROWNIAN','dt/tauB':float('nan')},
])
def test_unsupported_configuration_rejected(runner,tmp_path,override):
    with pytest.raises(ValueError):
        runner.run(tmp_path,override)
    assert not list(tmp_path.iterdir())


def test_network_real_execution_and_provenance(runner,tmp_path):
    result = runner.run(tmp_path,{'Np':64,'Nt':101,'dt/tauR':.001,'Rt':.001,
                                'box_dimension':6,'LBx':6,'LBy':6,'LBz':6,
                                'N_THREADS_BD':1,'N_skip_ener':10,'N_skip_file':10,
                                'record_RDIST':False})
    assert result.positions.shape == (10,64,3)
    assert result.energy_raw.shape == (10,31)
    assert result.bridge_counts.max() > 0
    assert np.allclose(result.times,np.arange(10)*.01)
    assert np.isfinite(result.positions).all()
    record = json.loads((result.directory/'run.json').read_text())
    assert record['status'] == 'passed'
    assert record['executable_sha256'] == hashlib.sha256(EXE.read_bytes()).hexdigest()
    assert record['parameters']['Np'] == '64'
    chains = np.loadtxt(result.directory/'case.chain',ndmin=2).astype(int)
    nc = 64*10
    independent_counts = (chains[:,:nc] != chains[:,nc:]).sum(axis=1)
    np.testing.assert_array_equal(independent_counts,result.bridge_counts)


def test_failed_runtime_is_recorded(tmp_path):
    executable = tmp_path/'fail.sh'
    executable.write_text('#!/bin/sh\nexit 9\n')
    executable.chmod(0o755)
    runner = LegacyRunner(executable,CODE/'short_test.inp')
    with pytest.raises(RuntimeError):
        runner.run(tmp_path/'runs')
    records = list((tmp_path/'runs').glob('legacy-*/run.json'))
    assert len(records)==1
    assert json.loads(records[0].read_text())['returncode']==9


def test_explicit_template_friction_is_rejected_before_execution(tmp_path):
    template = tmp_path/'explicit.inp'
    template.write_text((CODE/'short_test.inp').read_text()+'\nfriction_junction=LOOP_DISSOCIATION_TIME\n')
    runner = LegacyRunner('/bin/false',template)
    with pytest.raises(ValueError,match='friction'):
        runner.run(tmp_path/'runs')
    assert not (tmp_path/'runs').exists()
