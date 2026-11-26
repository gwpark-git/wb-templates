"""Explicit scalar-observable adapter from recorded ensembles to ML datasets."""
from __future__ import annotations

import hashlib
import json

import numpy as np

from .ml import Dataset,RegressionTask


def ensemble_dataset(report, *, feature_names, target_names, task=None,
                     failure_policy='raise',metadata=None):
    """Keep run, condition and continuation identity across every sampled row.

    Numeric specification fields plus time/step may be features. Seeds remain
    provenance. Failed jobs require an explicit exclusion policy and stay listed
    in metadata; excluding them does not establish unbiased physical sampling.
    """
    if failure_policy not in ('raise','exclude'):
        raise ValueError('failure_policy must be raise or exclude')
    payload=report.to_dict() if hasattr(report,'to_dict') else report
    if payload.get('schema')!='bdmc.ensemble.v1':
        raise ValueError('unsupported ensemble schema')
    features=tuple(feature_names);targets=tuple(target_names)
    if not features or not targets or set(features)&{'seed','run_id','condition_id','lineage_id'}:
        raise ValueError('explicit physical feature/target names required; seeds/identities are provenance')
    results=payload['results']
    ids=[r['run_id'] for r in results]
    if len(ids)!=len(set(ids)):
        raise ValueError('duplicate ensemble run_id')
    rows,labels,groups,conditions,lineages,provenance,excluded=[],[],[],[],[],[],[]
    for result in sorted(results,key=lambda r:r['run_id']):
        # Failure policy belongs to dataset construction, not an implicit retry
        # in the simulator. Preserve exclusions so later analysis can assess the
        # selection; finite-realization targets are not exact ensemble means.
        if result['status']!='success':
            if failure_policy=='raise':
                raise ValueError(f"ensemble contains an unsuccessful job: {result['run_id']}")
            excluded.append({'run_id':result['run_id'],'status':result['status'],'error':result.get('error')})
            continue
        if not result['core_hash_verified']:
            raise ValueError('successful job has no verified core identity')
        for key in ('run_id','condition_id','lineage_id'):
            if not isinstance(result[key],str) or not result[key]:
                raise ValueError('ensemble identities must be nonempty strings')
        summary=result['summary']
        times=np.asarray(summary['times'],dtype=float)
        steps=np.asarray(summary['steps'])
        if (times.ndim!=1 or not len(times) or steps.shape!=times.shape or
            not np.isfinite(times).all() or np.any(np.diff(times)<=0) or
            not np.issubdtype(steps.dtype,np.integer) or np.any(steps<0) or np.any(steps[1:]<=steps[:-1])):
            raise ValueError('successful job needs nonempty increasing sampled times/steps')
        columns={}
        for target in targets:
            if target not in summary['values']:
                raise ValueError(f'missing requested observable {target}')
            column=np.asarray(summary['values'][target],dtype=float)
            if column.shape!=times.shape or not np.isfinite(column).all():
                raise ValueError('dataset adapter requires finite scalar observables')
            columns[target]=column
        for frame,(t,n) in enumerate(zip(times,steps)):
            row=[]
            for feature in features:
                value=t if feature=='time' else n if feature=='step' else result['spec'].get(feature)
                if isinstance(value,bool) or not isinstance(value,(int,float,np.number)) or not np.isfinite(value):
                    raise ValueError(f'feature {feature} is not a finite numeric spec field/time/step')
                row.append(float(value))
            rows.append(row);labels.append([columns[target][frame] for target in targets])
            groups.append(result['run_id']);conditions.append(result['condition_id']);lineages.append(result['lineage_id'])
        provenance.append({key:result[key] for key in ('run_id','condition_id','lineage_id','seed','spec','core_binary_sha256')})
    if not rows:
        raise ValueError('ensemble contains no usable successful rows')
    info={'schema':'bdmc.ensemble-dataset.v1','source_schema':payload['schema'],
          'source_report_sha256':hashlib.sha256(json.dumps(payload,sort_keys=True,allow_nan=False).encode()).hexdigest(),
          'runs':provenance,'excluded_jobs':excluded,'failure_policy':failure_policy,
          'observable_definitions':{name:next(r['observable_definitions'][name] for r in results if r['status']=='success')
                                    for name in targets},
          'sampling_note':'finite-realization observables; physical training validity is a separate question',
          'user_metadata':dict(metadata or {})}
    return Dataset(np.asarray(rows),np.asarray(labels),np.asarray(groups),
        feature_names=features,target_names=targets,
        group_keys={'condition':np.asarray(conditions),'continuation':np.asarray(lineages)},
        task=task or RegressionTask(name='ensemble_observable',target_kind='realization'),metadata=info)
