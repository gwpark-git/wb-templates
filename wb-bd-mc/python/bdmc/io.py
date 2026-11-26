"""Versioned numeric dataset exchange without executable serialization."""
from dataclasses import asdict
from contextlib import nullcontext
import hashlib
import json
from pathlib import Path
import uuid

import numpy as np


def _estimated_npz_bytes(arrays):
    raw=sum(value.nbytes for value in arrays.values())
    # Conservative deflate expansion, per-member NPY/ZIP headers and filesystem
    # block headroom. Array names are fixed short keys, not user feature names.
    return raw+raw//1000+4096*len(arrays)+65_536


def save_dataset(directory, dataset, *, budget=None):
    """Save a Dataset, including its scientific task and complete group identity."""
    directory=Path(directory)
    keys=sorted(dataset.group_keys)
    arrays={'X':dataset.X,'y':dataset.y,'groups':dataset.groups}
    arrays.update({f'group_{i}':dataset.group_keys[key] for i,key in enumerate(keys)})
    record={'schema':'bdmc.dataset.v1','feature_names':dataset.feature_names,
            'target_names':dataset.target_names,'task':asdict(dataset.task),
            'metadata':dataset.metadata,'group_keys':keys,'fingerprint':dataset.fingerprint,
            'data_sha256':'0'*64}
    manifest_bytes=len((json.dumps(record,indent=2,allow_nan=False)+'\n').encode())
    planned=_estimated_npz_bytes(arrays)+manifest_bytes
    if budget is not None:
        budget.require_output(directory,planned_bytes=planned)
    with budget.writer(directory,planned_bytes=planned) if budget is not None else nullcontext():
        directory.mkdir(parents=True,exist_ok=False)
        np.savez_compressed(directory/'data.npz',**arrays)
        record['data_sha256']=hashlib.sha256((directory/'data.npz').read_bytes()).hexdigest()
        (directory/'dataset.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    return directory


def load_dataset(directory):
    from .ml import Dataset,RegressionTask
    directory=Path(directory)
    record=json.loads((directory/'dataset.json').read_text())
    if record['schema']!='bdmc.dataset.v1':
        raise ValueError('unsupported dataset schema')
    if hashlib.sha256((directory/'data.npz').read_bytes()).hexdigest()!=record['data_sha256']:
        raise ValueError('dataset numeric checksum mismatch')
    keys=record['group_keys']
    if len(set(keys))!=len(keys):
        raise ValueError('duplicate dataset group key')
    task={**record['task'],'target_units':tuple(record['task']['target_units'])}
    with np.load(directory/'data.npz',allow_pickle=False) as arrays:
        dataset=Dataset(arrays['X'],arrays['y'],arrays['groups'],
            feature_names=tuple(record['feature_names']),target_names=tuple(record['target_names']),
            task=RegressionTask(**task),metadata=record['metadata'],
            group_keys={key:arrays[f'group_{i}'] for i,key in enumerate(keys)})
    if dataset.fingerprint!=record['fingerprint']:
        raise ValueError('dataset metadata/fingerprint mismatch')
    return dataset


_CHECKPOINT_ARRAYS = {
    'particle_ids':('state','particle_ids'),
    'positions':('state','positions'),
    'unwrapped_positions':('state','unwrapped_positions'),
    'fixed_bonds':('topology','fixed_bonds'),
    'dynamic_bonds':('topology','dynamic_bonds'),
}


def save_checkpoint(directory, simulation, *, budget=None, lineage_id=None):
    """Save physical/RNG state, distinct from an observation frame.

    The native core rejects unsupported custom components. This file format
    contains JSON metadata and numeric arrays only, never executable pickle.
    """
    directory=Path(directory)
    if directory.exists():
        raise FileExistsError(directory)
    lineage_id=str(uuid.uuid4()) if lineage_id is None else lineage_id
    if not isinstance(lineage_id,str) or not lineage_id.strip():
        raise ValueError('checkpoint lineage_id must be a nonempty string')
    snapshot=simulation.checkpoint()
    # A lineage label accompanies actual restart state; it cannot recreate that
    # state by itself or turn continuation frames into independent replicas.
    arrays={name:np.array(snapshot[group][key],copy=True)
            for name,(group,key) in _CHECKPOINT_ARRAYS.items()}
    # Keep the native schema intact while moving its five arrays into NPZ.
    metadata={key:value for key,value in snapshot.items() if key not in ('state','topology')}
    metadata['state']={key:value for key,value in snapshot['state'].items()
                       if ('state',key) not in _CHECKPOINT_ARRAYS.values()}
    metadata['topology']={key:value for key,value in snapshot['topology'].items()
                          if ('topology',key) not in _CHECKPOINT_ARRAYS.values()}
    from . import _core
    record={'schema':'bdmc.checkpoint-file.v1','snapshot':metadata,'lineage_id':lineage_id,
            'core_binary_sha256':hashlib.sha256(Path(_core.__file__).read_bytes()).hexdigest(),
            'compatibility':'exact replay is tested for the same backend/version; cross-build replay requires verification',
            'state_sha256':'0'*64}
    planned=_estimated_npz_bytes(arrays)+len((json.dumps(record,indent=2,allow_nan=False)+'\n').encode())
    if budget is not None:
        budget.require_output(directory,planned_bytes=planned)
    with budget.writer(directory,planned_bytes=planned) if budget is not None else nullcontext():
        directory.mkdir(parents=True,exist_ok=False)
        np.savez_compressed(directory/'state.npz',**arrays)
        record['state_sha256']=hashlib.sha256((directory/'state.npz').read_bytes()).hexdigest()
        (directory/'checkpoint.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    return directory


def read_checkpoint(directory):
    """Read the numeric snapshot; native restore validates coupled invariants."""
    directory=Path(directory)
    record=json.loads((directory/'checkpoint.json').read_text())
    if record['schema']!='bdmc.checkpoint-file.v1':
        raise ValueError('unsupported checkpoint file schema')
    if not isinstance(record.get('lineage_id'),str) or not record['lineage_id'].strip():
        raise ValueError('checkpoint lineage_id must be a nonempty string')
    if hashlib.sha256((directory/'state.npz').read_bytes()).hexdigest()!=record['state_sha256']:
        raise ValueError('checkpoint numeric checksum mismatch')
    snapshot=record['snapshot']
    with np.load(directory/'state.npz',allow_pickle=False) as arrays:
        if set(arrays.files)!=set(_CHECKPOINT_ARRAYS):
            raise ValueError('checkpoint array schema mismatch')
        for name,(group,key) in _CHECKPOINT_ARRAYS.items():
            value=np.array(arrays[name],copy=True)
            if value.dtype.kind not in 'ifu' or not np.isfinite(value).all():
                raise ValueError('checkpoint arrays must be finite numeric values')
            snapshot[group][key]=value
    return snapshot


def load_checkpoint(directory):
    # Reconstruct the native schema and let the core validate coupled clock,
    # topology, coordinate and RNG invariants. NPZ integrity alone is weaker.
    from . import Simulation
    return Simulation.from_checkpoint(read_checkpoint(directory))
