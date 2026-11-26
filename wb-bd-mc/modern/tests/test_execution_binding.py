"""Explicit backend selection and the public model/metadata boundary."""
import numpy as np
import pytest
from bdmc.reference import (DumbbellSpec,DumbbellSimulation,RepulsiveSpec,
                            RepulsiveSimulation,execution_backends)


def test_capabilities_and_selection_never_silently_fall_back():
    caps={c['name']:c for c in execution_backends()}
    assert set(caps)=={'serial','openmp','cuda','metal','accelerate'}
    assert caps['serial']['available'] and caps['serial']['precision']=='float64'
    for name,c in caps.items():
        if not c['available']:
            with pytest.raises(ValueError,match='unavailable'):
                DumbbellSpec(execution_backend=name)
    with pytest.raises(ValueError):
        DumbbellSpec(execution_backend='auto')
    with pytest.raises(ValueError):
        DumbbellSpec(backend_threads=2)


@pytest.mark.parametrize('law',['hookean','fene'])
@pytest.mark.parametrize('threads',[1,2])
def test_openmp_same_noise_and_explicit_metadata(law,threads):
    if not next(c for c in execution_backends() if c['name']=='openmp')['available']:
        pytest.skip('OpenMP backend not compiled in this edition')
    kwargs=dict(n_dumbbells=17,spring_law=law,q0=4.,seed=931501,protocol='simple_shear',shear_rate=-.4)
    a=DumbbellSimulation(DumbbellSpec(**kwargs))
    b=DumbbellSimulation(DumbbellSpec(**kwargs,execution_backend='openmp',backend_threads=threads))
    aa=a.sample(23,stride=7);bb=b.sample(23,stride=7)
    np.testing.assert_array_equal(aa['steps'],bb['steps'])
    np.testing.assert_allclose(aa['connectors'],bb['connectors'],rtol=1e-12,atol=1e-13)
    assert b.spec.to_dict()['execution_backend']=='openmp'
    assert bb['metadata']['execution_backend']=='openmp'
    assert bb['metadata']['arithmetic_precision']=='float64'


def test_nonserial_repulsive_requires_explicit_all_pairs_choice():
    if not next(c for c in execution_backends() if c['name']=='openmp')['available']:
        pytest.skip('OpenMP backend not compiled in this edition')
    with pytest.raises(ValueError,match='brute_force'):
        RepulsiveSpec(execution_backend='openmp')
    kwargs=dict(n_particles=8,box_length=2.,seed=931503,epsilon=2.)
    a=RepulsiveSimulation(RepulsiveSpec(**kwargs))
    b=RepulsiveSimulation(RepulsiveSpec(**kwargs,execution_backend='openmp',force_search='brute_force',backend_threads=2))
    a.advance(20);b.advance(20)
    np.testing.assert_allclose(a.observe()['forces'],b.observe()['forces'],rtol=1e-11,atol=1e-12)
