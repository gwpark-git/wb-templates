#!/usr/bin/env python3
"""Build/test an explicit backend in fresh budgeted directories.

Never installs a toolkit, changes a branch, selects a silent fallback, or resets
the data allowance. --probe is read-only. Apple presets require an Apple
Silicon macOS host; a Linux ARM host is not an Apple toolchain.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import importlib.util
import json
import os
from pathlib import Path
import platform
import shutil
import sys
import tempfile

CODE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CODE / 'python'))
from bdmc.storage import DataBudget

PRESETS = ('cpu','openmp','mpi','cuda','apple-accelerate','apple-metal','apple-mpi','legacy-cxx17')


def probe():
    apple = sys.platform == 'darwin' and platform.machine().lower() in ('arm64','aarch64')
    tools = {name: shutil.which(name) for name in ('cmake','mpicxx','mpiexec','nvcc','xcrun','clang++')}
    local_cmake = Path(sys.executable).absolute().parent/'cmake'
    if local_cmake.is_file(): tools['cmake'] = str(local_cmake)
    return {'schema':'bdmc.backend-toolchains.v1','platform':sys.platform,
            'machine':platform.machine(),'apple_silicon_macos':apple,'tools':tools,
            'toolchain_candidates':{name:(
                apple and bool(tools['xcrun']) and (name!='apple-mpi' or bool(tools['mpicxx']))
                if name.startswith('apple-') else bool(tools['nvcc']) if name=='cuda' else
                bool(tools['mpicxx']) if name=='mpi' else True) for name in PRESETS},
            'limits':'Tool discovery is not successful compilation, device availability, execution, or performance evidence.'}


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--backend',choices=PRESETS,default='cpu')
    parser.add_argument('--probe',action='store_true')
    parser.add_argument('--test',action='store_true')
    parser.add_argument('--native-only',action='store_true',help='build without Python extensions')
    parser.add_argument('--data-budget',type=Path,default=CODE/'DATA_BUDGET.json')
    parser.add_argument('--compiler',default=None)
    args=parser.parse_args(argv)
    info=probe()
    if args.probe:
        print(json.dumps(info,indent=2));return 0
    if not info['tools']['cmake'] or not info['toolchain_candidates'][args.backend]:
        print(json.dumps({'status':'unavailable','requested_backend':args.backend,'probe':info},indent=2))
        return 2
    budget=DataBudget(args.data_budget)
    run_root=CODE/'runs/validation'
    budget.require_output(run_root,512_000_000)
    with budget.writer(run_root,512_000_000):
        run_root.mkdir(parents=True,exist_ok=True)
        run=Path(tempfile.mkdtemp(prefix='backend-'+args.backend+'-',dir=run_root))
        build=CODE/'build/validation'/run.name
        budget.require_output(build,512_000_000)
    # Reuse bounded command/log/cleanup and source-hash guards without importing
    # or mutating an installed numerical extension.
    spec=importlib.util.spec_from_file_location('_backend_validator',CODE/'tools/validate.py')
    validator=importlib.util.module_from_spec(spec);spec.loader.exec_module(validator)
    report={'schema':'bdmc.backend-build.v1','created_at':datetime.now(timezone.utc).isoformat(),
            'requested_backend':args.backend,'status':'running','probe':info,
            'run_directory':str(run.resolve()),'build_directory':str(build.resolve()),
            'source_before':validator.source_identity(),'steps':[]}
    env={**os.environ,'CMAKE_BUILD_PARALLEL_LEVEL':'2','OMP_NUM_THREADS':'1',
         'OPENBLAS_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1',
         'BDMC_OUTER_BUDGET_MANAGED':'1'}
    cmake=info['tools']['cmake']
    command=[cmake,'--preset',args.backend,'-B',str(build.resolve()),
             '-DPython_EXECUTABLE='+str(Path(sys.executable).absolute())]
    if args.compiler: command+=['-DCMAKE_CXX_COMPILER='+args.compiler]
    if args.native_only: command+=['-DBDMC_BUILD_PYTHON=OFF']
    commands=[('configure',command),('build',[cmake,'--build',str(build.resolve()),'-j','2'])]
    if args.test:
        ctest=str(Path(cmake).parent/'ctest') if Path(cmake).is_absolute() else 'ctest'
        commands.append(('native_tests',[ctest,'--test-dir',str(build.resolve()),'--output-on-failure','--no-tests=error']))
    try:
        with budget.writer(build,512_000_000):
            for name,command in commands:
                result=validator.execute_step(command,run/(name+'.log'),cwd=CODE/'modern',env=env,timeout=600)
                report['steps'].append(dict(result,name=name))
                if result['status']!='passed': raise RuntimeError(name+': '+result['error'])
        report['source_after']=validator.source_identity()
        if report['source_after']!=report['source_before']:
            raise RuntimeError('sources/tests/config changed during build; verify a stable edition')
        report['status']='passed'
    except Exception as exc:
        report.update(status='failed',error=f'{type(exc).__name__}: {exc}')
    report['completed_at']=datetime.now(timezone.utc).isoformat()
    report['storage_after']=budget.check()
    with budget.writer(run/'report.json',4_000_000):
        (run/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(run/'report.json',report['status'])
    return 0 if report['status']=='passed' else 1


if __name__=='__main__':
    raise SystemExit(main())
