import numpy as np
import pytest

from bdmc.experiments import calibration_dataset, dumbbell_second_moment
from bdmc.ml import grouped_split


def test_relative_noise_and_simultaneous_pair_drift():
    # One coordinate, kappa=2^2=4, q0=1, a=1-2*mu*k*dt=.84.
    expected=.84**2+4*3*.02
    assert dumbbell_second_moment(time_steps=1,dt=.02,dimension=1,mobility=1,
                                 diffusivity=3,spring_alpha=2)==pytest.approx(expected)


def test_force_free_relative_diffusion():
    assert dumbbell_second_moment(time_steps=20,dt=.01,dimension=3,mobility=0,
        diffusivity=2,spring_alpha=1)==pytest.approx(1+4*2*3*.2)


@pytest.mark.parametrize('override',[{'dt':float('nan')},{'diffusivity':float('inf')},
                                    {'initial_q_squared':-1},{'time_steps':True}])
def test_invalid_analytic_parameters(override):
    values=dict(time_steps=2,dt=.01,dimension=3,mobility=1,diffusivity=1,spring_alpha=1)
    with pytest.raises(ValueError):
        dumbbell_second_moment(**{**values,**override})


@pytest.mark.parametrize('override',[{'spring_alpha':1e154,'mobility':0},
    {'diffusivity':1e308,'dt':1.,'mobility':0},{'time_steps':2**64}])
def test_derived_moment_overflow_rejected(override):
    values=dict(time_steps=1,dt=.01,dimension=3,mobility=1,diffusivity=1,spring_alpha=1)
    with pytest.raises(ValueError):
        dumbbell_second_moment(**{**values,**override})


def test_long_horizon_moment_has_constant_work():
    # The stable finite-step limit has an independent closed form.
    expected=1/(1-.03)  # d*D / (mu*kappa*(1-mu*kappa*dt)), kappa=3.
    value=dumbbell_second_moment(time_steps=10**15,dt=.01,dimension=3,
                                 mobility=1,diffusivity=1,spring_alpha=1)
    assert value==pytest.approx(expected)


def test_real_dataset_has_condition_and_run_provenance():
    dataset,theory,counts=calibration_dataset('brownian',particles=512,steps=24,stride=6,
                                           seeds_per_condition=1,conditions=6)
    assert counts=={'n_runs':6,'n_rows':24}
    assert dataset.metadata['core_binary_sha256']
    assert len(dataset.metadata['runs'])==6
    # Ensemble moments fluctuate; this checks normalization against known theory.
    assert np.max(np.abs(dataset.y/theory-1)) < .2
    split=grouped_split(dataset,seed=19,group_by=('condition',))
    labels=[set(dataset.group_keys['condition'][indices]) for indices in
            (split.train,split.validation,split.test)]
    assert not labels[0]&labels[1] and not labels[0]&labels[2] and not labels[1]&labels[2]
