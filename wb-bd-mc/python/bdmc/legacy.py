"""Recorded subprocess compatibility for an explicitly selected legacy binary.

This adapter is intentionally separate from the direct modern-core binding.
Legacy energy rows mix recording phases and are exposed without relabelling.
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
from typing import Mapping

import numpy as np


def read_input(path: str | Path) -> dict[str, str]:
    values = {}
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if line and not line.startswith('#') and '=' in line:
            key, value = line.split('=', 1)
            if key.strip() in values:
                raise ValueError(f'duplicate legacy key: {key.strip()}')
            values[key.strip()] = value.strip()
    return values


def _hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@dataclass(frozen=True)
class LegacyResult:
    directory: Path
    parameters: dict[str, str]
    metadata: dict
    times: np.ndarray
    positions: np.ndarray
    auxiliary: np.ndarray
    energy_raw: np.ndarray

    @property
    def bridge_counts(self) -> np.ndarray:
        if self.parameters['Method'] != 'NAPLE_ASSOCIATION' or self.parameters['Step'] == 'EQUILIBRATION':
            raise ValueError('bridge counts require the transient-network backend')
        return self.energy_raw[:, 4].copy()


class LegacyRunner:
    """Run an explicitly selected executable and input template in fresh folders.

    Version one supports 3D, quiescent, fresh starts and one topology thread.
    These bounds describe inspected adapter coverage, not every legacy feature.
    """
    def __init__(self, executable: str | Path, template: str | Path):
        self.executable = Path(executable).resolve()
        self.template = Path(template).resolve()
        if not self.executable.is_file():
            raise FileNotFoundError(self.executable)
        self.defaults = read_input(self.template)

    def run(self, output_root: str | Path, overrides: Mapping[str, object] | None = None,
            *, timeout: float = 60) -> LegacyResult:
        values = dict(self.defaults)
        values.setdefault('probability_area', 'FALSE')
        values.setdefault('friction_junction', 'DEFAULT')
        overrides = dict(overrides or {})
        if 'output_path' in overrides or 'filename_base' in overrides:
            raise ValueError('adapter owns output_path and filename_base')
        unknown = set(overrides) - set(values) - {'dt/tauB'}
        if unknown:
            raise ValueError(f'unknown legacy keys: {sorted(unknown)}')
        for key, value in overrides.items():
            if isinstance(value, bool):
                value = 'TRUE' if value else 'FALSE'
            value = str(value)
            if '\n' in value or '\r' in value:
                raise ValueError('legacy values must be single lines')
            values[key] = value
        if values['Method'] not in {'BROWNIAN', 'DUMBBELL', 'REPULSIVE_BROWNIAN', 'NAPLE_ASSOCIATION'}:
            raise ValueError('unknown Method')
        if values['Integrator'] != 'Euler':
            raise ValueError('legacy backend supports only Euler')
        if values['Step'] not in {'EQUILIBRATION','NO_EQUILIBRATION'}:
            raise ValueError('unknown legacy Step')
        if int(values['N_dimension']) != 3:
            raise ValueError('legacy adapter coverage is 3D')
        if values['MECHANICAL_PERTURBATION'] != 'FALSE':
            raise ValueError('shear is outside verified adapter coverage')
        if values['friction_junction'] != 'DEFAULT':
            raise ValueError('nondefault friction requires separate adapter acceptance for the selected executable')
        if int(values['N_THREADS_SS']) != 1:
            raise ValueError('parallel topology is outside verified adapter coverage')
        if any(values[k] != 'FALSE' for k in ('CONTINUATION_TRAJ','CONTINUATION_CONNECTION','CONTINUATION_CHAIN')):
            raise ValueError('continuation is outside verified adapter coverage')
        for k in ('Np','Nt','N_skip_file','N_skip_ener','N_THREADS_BD'):
            if int(values[k]) <= 0:
                raise ValueError(f'{k} must be positive')
        if int(values['Nt']) < 2:
            raise ValueError('Nt must be >=2')
        if values['N_skip_file'] != values['N_skip_ener']:
            raise ValueError('aligned energy/trajectory rows require equal skip values')
        if values['Method'] == 'DUMBBELL' and (int(values['Np']) % 2 or values['cell_list'] != 'FALSE'):
            raise ValueError('legacy dumbbells require even Np and cell_list=FALSE')
        dt_key = 'dt/tauB' if values['Method'] in {'BROWNIAN','DUMBBELL'} else 'dt/tauR'
        dt = float(values.get(dt_key, 'nan'))
        if not np.isfinite(dt) or dt <= 0:
            raise ValueError(f'positive finite {dt_key} required')
        root = Path(output_root).resolve()
        root.mkdir(parents=True, exist_ok=True)
        run = Path(tempfile.mkdtemp(prefix='legacy-', dir=root))
        values.update(output_path=str(run), filename_base='case')
        inp = run/'effective.inp'
        inp.write_text(''.join(f'{k}={v}\n' for k,v in values.items()))
        metadata = {'schema':'bdmc.legacy-run.v1', 'created_at':datetime.now(timezone.utc).isoformat(),
                    'executable':str(self.executable), 'executable_sha256':_hash(self.executable),
                    'template':str(self.template), 'template_sha256':_hash(self.template),
                    'effective_input_sha256':_hash(inp), 'parameters':values,
                    'backend':'legacy_subprocess', 'noise':'legacy_uniform_variance_normalised',
                    'recording':'trajectory is old state; raw observable columns retain the selected executable semantics',
                    'validation':'bounded engineering checks, historical regression deferred'}
        env = {**os.environ, 'OMP_DYNAMIC':'FALSE', 'GSL_RNG_TYPE':'mt19937'}
        try:
            with (run/'stdout.log').open('w') as log:
                result = subprocess.run([str(self.executable),str(inp)],cwd=run,env=env,
                                        stdout=log,stderr=subprocess.STDOUT,timeout=timeout)
            metadata['returncode'] = result.returncode
            logtext = (run/'stdout.log').read_text()
            if result.returncode or 'ERR:' in logtext or 'Bad condition' in logtext:
                raise RuntimeError(f'legacy execution failed; inspect {run}/stdout.log')
            traj = np.loadtxt(run/'case.traj', ndmin=2)
            energy = np.loadtxt(run/'case.ener', ndmin=2)
            n, dim = int(values['Np']), int(values['N_dimension'])
            frames = len(range(0,int(values['Nt'])-1,int(values['N_skip_file'])))
            association = values['Method']=='NAPLE_ASSOCIATION' and values['Step']!='EQUILIBRATION'
            columns = 31 if association else 18 if values['Method']=='BROWNIAN' else 24
            if traj.shape != (frames,1+2*n*dim) or energy.shape != (frames,columns):
                raise ValueError('legacy output schema/record count mismatch')
            if not np.isfinite(traj).all() or not np.isfinite(energy).all():
                raise ValueError('nonfinite legacy output')
            expected = np.arange(frames)*int(values['N_skip_file'])*dt
            if not np.allclose(traj[:,0],expected,atol=1e-8,rtol=1e-7):
                raise ValueError('legacy trajectory timestamps mismatch')
            if np.allclose(energy[:,0],traj[:,0],atol=1e-8,rtol=1e-7):
                metadata['energy_clock']='aligned pre-step force/frame timestamps observed'
            elif values['Method']=='BROWNIAN' and np.all(energy[:,0]==0):
                metadata['energy_clock']='historical zero Brownian energy clock observed'
                metadata['energy_time_warning']='selected executable has zero Brownian energy time; use trajectory times'
            else:
                raise ValueError('legacy energy timestamps mismatch')
            particle = traj[:,1:].reshape(frames,n,2*dim)
            # Preserve the selected recorder's auxiliary/raw columns. Renaming
            # them as modern velocity/stress would imply physical semantics that
            # the compatibility adapter has not established.
            metadata.update(status='passed',files={p.name:_hash(p) for p in run.glob('case.*')})
            return LegacyResult(run,values,metadata,traj[:,0].copy(),particle[:,:,:dim].copy(),
                                particle[:,:,dim:].copy(),energy)
        except Exception as exc:
            metadata.update(status='failed',error=str(exc))
            raise
        finally:
            (run/'run.json').write_text(json.dumps(metadata,indent=2)+'\n')
