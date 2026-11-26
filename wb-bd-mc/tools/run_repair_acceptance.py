#!/usr/bin/env python3
"""Queue the repair acceptance stages through their existing budget guards.

This coordinator never holds the shared writer while launching a child that
owns that writer. It stops at the first failure and preserves all evidence.
"""
from datetime import datetime, timezone
import argparse
import importlib.util
import json
import os
from pathlib import Path
import shutil
import sys
import time

CODE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(CODE/'python'))
from bdmc.storage import DataBudget


def validate_unit_report(data, validator, compiler):
    """Admit only this exact stable source/build edition and required matrix."""
    if data.get('status') != 'passed' or data.get('profile') != 'unit':
        raise ValueError('acceptance requires a passed unit profile')
    if Path(data.get('compiler', '')).resolve() != Path(compiler).resolve():
        raise ValueError('acceptance requires the selected GCC compiler')
    source = validator.source_identity()
    if data.get('source_before') != source or data.get('source_after') != source:
        raise ValueError('unit source edition does not match current source')
    stages = {step.get('name'): step for step in data.get('steps', [])}
    required = ('legacy_build_and_units', 'modern_configure', 'modern_build',
                'modern_native', 'python_import', 'python_tests')
    if any(stages.get(name, {}).get('status') != 'passed' for name in required):
        raise ValueError('unit matrix is incomplete')
    configure = stages['modern_configure'].get('command', [])
    if not all(flag in configure for flag in ('-DBDMC_ENABLE_OPENMP=ON',
                   '-DBDMC_ENABLE_MPI=ON', '-DBDMC_BUILD_LEGACY=ON')):
        raise ValueError('unit matrix lacks OpenMP/MPI/legacy CMake coverage')
    build = Path(data['build_directory']).resolve()
    executable = build/'legacy/stochastic_simulation'
    if validator.digest(executable) != data.get('legacy_executable_sha256'):
        raise ValueError('unit legacy executable changed')
    package = build/'import/bdmc'
    if Path(data.get('staged_package', '')).resolve() != package:
        raise ValueError('unit staged package path changed')
    for name in ('_core', '_reference', '_network', '_analysis', '_hydrodynamics'):
        entry = data.get('native_extensions', {}).get(name, {})
        path = Path(entry.get('path', '')).resolve()
        if path.parent != package or validator.digest(path) != entry.get('sha256'):
            raise ValueError('unit extension changed: '+name)



def main():
    raise SystemExit("Historical workspace coordinator is not a standalone workflow. Use tools/validate.py with an explicit profile and data budget; frozen references are never recreated by verification.")


if __name__ == "__main__":
    main()
