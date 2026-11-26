from copy import deepcopy

import numpy as np
import pytest

from bdmc.ensemble_data import ensemble_dataset
from bdmc.ml import grouped_split


def report():
    results=[]
    for condition in range(6):
        for part in range(2):
            results.append({'run_id':f'run{condition}-{part}','condition_id':f'condition{condition}',
                'lineage_id':f'lineage{condition}','seed':condition,'spec':{'diffusivity':condition+1.},
                'status':'success','core_hash_verified':True,'core_binary_sha256':'fixture',
                'summary':{'steps':[0,2],'times':[0.,.2],
                           'values':{'msd':[0.,1.2*(condition+1.)]}},
                'observable_definitions':{'msd':'unwrapped finite-realization mean square displacement'}})
    return {'schema':'bdmc.ensemble.v1','metadata':{},'results':results}


def test_adapter_preserves_continuation_and_condition_groups():
    dataset=ensemble_dataset(report(),feature_names=('diffusivity','time'),target_names=('msd',))
    assert dataset.X.shape==(24,2)
    np.testing.assert_allclose(dataset.y[:2,0],[0.,1.2])
    split=grouped_split(dataset,group_by=('condition','continuation'),seed=91)
    for key in ('groups','condition','continuation'):
        sets=[set(dataset.labels(key)[idx]) for idx in (split.train,split.validation,split.test)]
        assert not sets[0]&sets[1] and not sets[0]&sets[2] and not sets[1]&sets[2]
    assert dataset.metadata['source_report_sha256']


def test_failure_exclusion_is_explicit_and_recorded():
    r=report();r['results'][0]['status']='timeout'
    with pytest.raises(ValueError,match='unsuccessful'):
        ensemble_dataset(r,feature_names=('diffusivity','time'),target_names=('msd',))
    data=ensemble_dataset(r,feature_names=('diffusivity','time'),target_names=('msd',),failure_policy='exclude')
    assert len(data.X)==22
    assert data.metadata['excluded_jobs']==[{'run_id':'run0-0','status':'timeout','error':None}]


@pytest.mark.parametrize('feature',['seed','missing','model'])
def test_seed_identity_and_nonnumeric_features_rejected(feature):
    with pytest.raises(ValueError):
        ensemble_dataset(report(),feature_names=(feature,),target_names=('msd',))


def test_unknown_vector_observable_and_unverified_core_rejected():
    r=report()
    with pytest.raises(ValueError,match='missing'):
        ensemble_dataset(r,feature_names=('time',),target_names=('unknown',))
    r['results'][0]['summary']['values']['msd']=[[1,2],[3,4]]
    with pytest.raises(ValueError,match='scalar'):
        ensemble_dataset(r,feature_names=('time',),target_names=('msd',))
    r=report();r['results'][0]['core_hash_verified']=False
    with pytest.raises(ValueError,match='verified'):
        ensemble_dataset(r,feature_names=('time',),target_names=('msd',))
