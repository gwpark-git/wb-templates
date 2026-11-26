#!/usr/bin/env python3
"""Isolated native/Python/frozen-reference checks; never creates expected data.

Run with the project's Python environment. All optional ML skips are reported.
The sanitizer profile runs native and legacy-reference checks, without loading
an ASan module into Python. Logs, scratch outputs and compiler builds are kept.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import xml.etree.ElementTree as ET

CODE = Path(__file__).resolve().parents[1]

DIAGNOSTIC = re.compile(
    r'runtime error:|(?:ERROR|SUMMARY): (?:AddressSanitizer|UndefinedBehaviorSanitizer|LeakSanitizer)'
    r'|AddressSanitizer:DEADLYSIGNAL|\bBad condition\b|\bERR:', re.I)


def digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b''):
            value.update(chunk)
    return value.hexdigest()


def source_identity():
    paths = []
    for folder in ('src', 'lib', 'modern/src', 'modern/include', 'python/bdmc',
                   'tests', 'modern/tests', 'tools'):
        paths.extend(p for p in (CODE / folder).rglob('*')
                     if p.is_file() and p.suffix in ('.py', '.cpp', '.h', '.hpp', '.cu', '.mm', '.metal'))
    paths += [CODE / 'Makefile', CODE / 'modern/CMakeLists.txt', CODE / 'short_test.inp',
              CODE / 'pyproject.toml', CODE / 'modern/CMakePresets.json']
    paths.extend((CODE / 'python/bdmc').glob('*.json'))
    # Imported HI is a frozen, self-contained source dependency. Include its
    # numerical tables and license/provenance, never the other live worktree.
    paths.extend(p for p in (CODE / 'third_party/hi_native').rglob('*') if p.is_file())
    return {str(p.relative_to(CODE)): digest(p) for p in sorted(paths)}


def stop_process(process):
    if os.name == 'posix':
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
    elif process.poll() is None:
        process.kill()
    process.wait()


def execute_step(command, log, *, cwd=CODE, env=None, timeout=300, max_log_bytes=16_000_000,
                 expected_diagnostics=()):
    """Retain a bounded log; zero exit with sanitizer diagnostics still fails."""
    start = time.monotonic()
    result = {'command': [str(x) for x in command], 'cwd': str(cwd),
              'log': str(log), 'status': 'running'}
    process = None
    try:
        with Path(log).open('wb') as out:
            process = subprocess.Popen(result['command'], cwd=cwd, env=env,
                                       stdout=out, stderr=subprocess.STDOUT,
                                       start_new_session=(os.name == 'posix'))
            while process.poll() is None:
                if time.monotonic() - start > timeout:
                    result['error'] = 'command timeout'
                    stop_process(process)
                    break
                if Path(log).stat().st_size > max_log_bytes:
                    result['error'] = 'command log exceeded byte limit'
                    stop_process(process)
                    break
                time.sleep(0.05)
            # Reap any descendants left behind by a completed command as well.
            stop_process(process)
        result['returncode'] = process.returncode
        # Read at most the bounded prefix if a fast producer exceeded the cap.
        with Path(log).open('rb') as stream:
            text = stream.read(max_log_bytes + 1).decode('utf-8', errors='replace')
        if Path(log).stat().st_size > max_log_bytes:
            result.setdefault('error', 'command log exceeded byte limit')
        expected = Counter(expected_diagnostics)
        scanned_lines, acknowledged = [], []
        for line in text.splitlines():
            if expected[line] > 0:
                expected[line] -= 1
                acknowledged.append(line)
            else:
                scanned_lines.append(line)
        result['expected_rejection_diagnostics'] = acknowledged
        if any(expected.values()):
            result.setdefault('error', 'declared negative-test rejection diagnostics were not observed')
        if DIAGNOSTIC.search('\n'.join(scanned_lines)):
            result.setdefault('error', 'runtime/sanitizer diagnostic in log')
        if process.returncode != 0:
            result.setdefault('error', f'exit {process.returncode}')
        result['known_limitations'] = re.findall(r'^XFAIL(?::|\s).*$', text, re.M)
        result['status'] = 'failed' if 'error' in result else 'passed'
    except BaseException as exc:
        if process is not None:
            stop_process(process)
        if not isinstance(exc, Exception):
            raise
        result.update(status='failed', error=f'{type(exc).__name__}: {exc}')
    result['elapsed_seconds'] = time.monotonic() - start
    return result


def executable(name, python):
    local = python.parent / name
    found = str(local) if local.is_file() else shutil.which(name)
    if not found:
        raise RuntimeError(f'{name} is required; select the prepared Python/toolchain environment')
    return found


def junit_summary(path):
    root = ET.parse(path).getroot()
    suites = [root] if root.tag == 'testsuite' else list(root.findall('testsuite'))
    return {field: sum(int(suite.get(field, 0)) for suite in suites)
            for field in ('tests', 'failures', 'errors', 'skipped')}


def import_check_script(package, module, reference_module=None, network_module=None,
                        additional_modules=None):
    # Explicit failures remain active under PYTHONOPTIMIZE/-O.
    script = ('import pathlib,bdmc,bdmc._core\n'
            f'if pathlib.Path(bdmc.__file__).parent != pathlib.Path({str(package)!r}):\n'
            '    raise RuntimeError("wrong Python package imported")\n'
            f'if pathlib.Path(bdmc._core.__file__) != pathlib.Path({str(module)!r}):\n'
            '    raise RuntimeError("wrong native extension imported")\n'
            'print(bdmc._core.__file__)\n')
    if reference_module is not None:
        script += ('import bdmc._reference\n'
                   f'if pathlib.Path(bdmc._reference.__file__) != pathlib.Path({str(reference_module)!r}):\n'
                   '    raise RuntimeError("wrong reference extension imported")\n'
                   'print(bdmc._reference.__file__)\n')
    if network_module is not None:
        script += ('import bdmc._network\n'
                   f'if pathlib.Path(bdmc._network.__file__) != pathlib.Path({str(network_module)!r}):\n'
                   '    raise RuntimeError("wrong network extension imported")\n'
                   'print(bdmc._network.__file__)\n')
    for name, path in (additional_modules or {}).items():
        if name not in ('_analysis', '_hydrodynamics'):
            raise ValueError('unsupported additional native module')
        script += (f'import bdmc.{name}\n'
                   f'if pathlib.Path(bdmc.{name}.__file__) != pathlib.Path({str(path)!r}):\n'
                   f'    raise RuntimeError("wrong {name} extension imported")\n'
                   f'print(bdmc.{name}.__file__)\n')
    return script


def main(argv=None):
    # Lazy import keeps helper imports from changing a test process's package path.
    sys.path.append(str(CODE / 'python'))
    from bdmc.storage import DataBudget
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profile', choices=('unit', 'regression', 'full', 'sanitizer'), default='full')
    parser.add_argument('--compiler', default='g++')
    parser.add_argument('--python', type=Path, default=Path(sys.executable))
    parser.add_argument('--jobs', type=int, choices=range(1, 5), default=2)
    parser.add_argument('--enable-openmp', action='store_true', help='exercise parallel bulk CPU kernels')
    parser.add_argument('--enable-mpi', action='store_true', help='build and execute actual MPI integration checks')
    parser.add_argument('--legacy-cmake', action='store_true', help='also build/test the C++17 original simulator through CMake')
    parser.add_argument('--timeout', type=int, default=300)
    parser.add_argument('--run-dir', type=Path, default=CODE / 'runs/validation')
    parser.add_argument('--build-root', type=Path, default=CODE / 'build/validation')
    parser.add_argument('--reference', type=Path, default=CODE / 'tests/reference/v2')
    parser.add_argument('--budget', type=Path, default=CODE / 'DATA_BUDGET.json')
    args = parser.parse_args(argv)
    if not 1 <= args.timeout <= 600:
        parser.error('--timeout must be in 1..600 seconds per step')
    budget = DataBudget(args.budget)
    # Conservative allowance includes compiler/debug files, logs and all scratch.
    budget.require_output(args.run_dir, 768_000_000)
    budget.require_output(args.build_root, 768_000_000)
    python = args.python.absolute()  # Preserve venv invocation; resolving symlinks loses its prefix.
    compiler = shutil.which(args.compiler)
    if not compiler or not python.is_file():
        parser.error('compiler and selected Python must exist')
    with budget.writer(args.run_dir, 768_000_000):
        args.run_dir.mkdir(parents=True, exist_ok=True)
        run = Path(tempfile.mkdtemp(prefix=args.profile + '-', dir=args.run_dir.resolve()))
        build = args.build_root.resolve() / run.name
        budget.require_output(build, 768_000_000)
        build.mkdir(parents=True)
    report = {'schema': 'bdmc.validation.v1', 'created_at': datetime.now(timezone.utc).isoformat(),
              'profile': args.profile, 'status': 'running', 'run_directory': str(run),
              'build_directory': str(build), 'compiler': compiler, 'python': str(python),
              'source_before': source_identity(), 'steps': [],
              'scope': 'Engineering/unit/current-port reference checks; historical SSD/Intel regression deferred.'}
    report_path = run / 'report.json'
    env = {**os.environ, 'OMP_DYNAMIC': 'FALSE', 'OMP_THREAD_LIMIT': '8',
           'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1',
           'ASAN_OPTIONS': 'detect_leaks=0:halt_on_error=1',
           'UBSAN_OPTIONS': 'halt_on_error=1:print_stacktrace=1'}
    env['BDMC_OUTER_BUDGET_MANAGED']='1'
    # Ambient user pytest selection must not silently shrink the validation set.
    for key in ('PYTEST_ADDOPTS', 'PYTEST_PLUGINS'):
        env.pop(key, None)
    env['PYTEST_DISABLE_PLUGIN_AUTOLOAD'] = '1'

    def save():
        with budget.writer(report_path, 1_000_000):
            temporary = report_path.with_suffix('.tmp')
            temporary.write_text(json.dumps(report, indent=2) + '\n')
            temporary.replace(report_path)

    def step(name, command, *, custom_env=None, managed_budget=False, expected_diagnostics=()):
        print(f'RUN: {name}', flush=True)
        log = run / (name + '.log')
        if managed_budget:
            budget.require_output(log, 128_000_000)
            result = execute_step(command, log, env=custom_env or env, timeout=args.timeout,
                                  expected_diagnostics=expected_diagnostics)
        else:
            with budget.writer(log, 768_000_000):
                result = execute_step(command, log, env=custom_env or env, timeout=args.timeout,
                                      expected_diagnostics=expected_diagnostics)
        result['name'] = name
        report['steps'].append(result)
        save()
        if result['status'] != 'passed':
            raise RuntimeError(f'{name}: {result["error"]}; see {log}')

    try:
        # One bounded build/native/Python batch under the existing writer.
        # Nested saves remain checked; standalone reference children own their
        # writer after this block, so no cross-process recursive lock occurs.
        with budget.writer(build, 768_000_000):
            save()
            legacy = build / 'legacy'
            flags = '-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer' if args.profile == 'sanitizer' else '-O2 -g'
            step('legacy_build_and_units', ['make', '-j', str(args.jobs), 'all', 'test-unit',
                  f'CXX={compiler}', f'BUILD_DIR={legacy}', f'CXXFLAGS={flags}',
                  *(['LDFLAGS=-fsanitize=address,undefined'] if args.profile == 'sanitizer' else [])],
                 expected_diagnostics=(
                     'ERR: Dumbbells require positive even Np and two HASH slots (self and partner)',) * 2)
            exe = legacy / 'stochastic_simulation'
            report['legacy_executable_sha256'] = digest(exe)
            if args.profile == 'full':
                step('legacy_smoke', [python, CODE / 'tests/regression.py', '--executable', exe,
                                     '--run-dir', run / 'legacy-smoke'])
            if args.profile != 'regression':
                cmake = executable('cmake', python)
                ctest = executable('ctest', python)
                modern = build / 'modern'
                with_python = args.profile != 'sanitizer'
                step('modern_configure', [cmake, '-S', CODE / 'modern', '-B', modern,
                     f'-DCMAKE_CXX_COMPILER={compiler}', '-DCMAKE_BUILD_TYPE=Release',
                     f'-DPython_EXECUTABLE={python}',
                     f'-DBDMC_BUILD_PYTHON={"ON" if with_python else "OFF"}',
                     f'-DBDMC_ENABLE_OPENMP={"ON" if args.enable_openmp else "OFF"}',
                     f'-DBDMC_ENABLE_MPI={"ON" if args.enable_mpi else "OFF"}',
                     f'-DBDMC_BUILD_LEGACY={"ON" if args.legacy_cmake else "OFF"}',
                     '-DBDMC_NETWORK_LEGACY_ORACLE=ON',
                     f'-DBDMC_ENABLE_SANITIZERS={"ON" if args.profile == "sanitizer" else "OFF"}'])
                step('modern_build', [cmake, '--build', modern, '-j', str(args.jobs)])
                native_env = {**env, 'ASAN_OPTIONS': 'detect_leaks=1:halt_on_error=1'}
                step('modern_native', [ctest, '--test-dir', modern, '--output-on-failure',
                                      '--no-tests=error'], custom_env=native_env)
                if with_python:
                    staged = build / 'import'
                    package = staged / 'bdmc'
                    with budget.writer(package, 128_000_000):
                        shutil.copytree(CODE / 'python/bdmc', package,
                                        ignore=shutil.ignore_patterns('__pycache__', '*.pyc', '*.so'))
                        staged_modules = {}
                        for name in ('_core', '_reference', '_network', '_analysis', '_hydrodynamics'):
                            modules = list(modern.glob(name + '*.so')) + list(modern.glob(name + '*.pyd'))
                            if len(modules) != 1:
                                raise RuntimeError(f'expected exactly one freshly built {name} extension')
                            destination = package / modules[0].name
                            shutil.copy2(modules[0], destination)
                            staged_modules[name] = destination
                        module = staged_modules['_core']
                        # Editable-install meta finders can precede PYTHONPATH. A
                        # local finder pins only this package to the staged build,
                        # including in spawn workers and CLI subprocesses.
                        (staged / 'sitecustomize.py').write_text(
                            'import importlib.machinery, sys\n'
                            'class _ValidationPackageFinder:\n'
                            '    def find_spec(self, fullname, path=None, target=None):\n'
                            '        if fullname == "bdmc":\n'
                            f'            return importlib.machinery.PathFinder.find_spec(fullname, [{str(staged)!r}])\n'
                            '        if fullname.startswith("bdmc."):\n'
                            '            return importlib.machinery.PathFinder.find_spec(fullname, path)\n'
                            'sys.meta_path.insert(0, _ValidationPackageFinder())\n')
                    pyenv = {**env, 'PYTHONPATH': str(staged),
                             'BDMC_LEGACY_EXECUTABLE': str(exe), 'BDMC_ANALYSIS_TOOL_DIR':str(legacy),
                             'PYTHONDONTWRITEBYTECODE': '1'}
                    if args.enable_mpi:
                        pyenv['BDMC_MPI_EXECUTABLE']=str(modern / 'bdmc_mpi_ensemble')
                    report['staged_package'] = str(package)
                    report['core_extension_sha256'] = digest(module)
                    report['native_extensions'] = {
                        name: {'path': str(path), 'sha256': digest(path)}
                        for name, path in staged_modules.items()}
                    step('python_import', [python, '-c', import_check_script(
                         package, module, staged_modules['_reference'], staged_modules['_network'],
                         {name:staged_modules[name] for name in ('_analysis','_hydrodynamics')})], custom_env=pyenv)
                    step('python_tests', [python, '-m', 'pytest', '-q', '-ra', '-c', CODE / 'pyproject.toml',
                         '-o', 'cache_dir=' + str(run / 'pytest-cache'),
                         '--basetemp=' + str(run / 'pytest-temp'),
                         '--junitxml=' + str(run / 'python-junit.xml'),
                         CODE / 'tests/python', CODE / 'modern/tests'], custom_env=pyenv)
                    report['python_tests'] = junit_summary(run / 'python-junit.xml')
                    for name, path in staged_modules.items():
                        if digest(path) != report['native_extensions'][name]['sha256']:
                            raise RuntimeError(f'staged {name} changed during the verification window')
        if args.legacy_cmake and args.profile in ('full','sanitizer'):
            step('cmake_legacy_reference', [python, CODE / 'tests/reference_suite.py', 'verify',
                 '--executable', modern / 'bdmc_legacy', '--dataset', args.reference.resolve(),
                 '--run-dir', run / 'cmake-legacy-reference', '--budget', args.budget.resolve(),
                 '--compiler', compiler, '--multithread-invariants'], managed_budget=True,
                 custom_env={**env, 'BDMC_REFERENCE_MANAGED_GROUP':'1'})
        if args.profile in ('regression', 'full', 'sanitizer'):
            step('frozen_reference', [python, CODE / 'tests/reference_suite.py', 'verify',
                 '--executable', exe, '--dataset', args.reference.resolve(),
                 '--run-dir', run / 'reference', '--budget', args.budget.resolve(),
                 '--compiler', compiler,
                 '--multithread-invariants'], managed_budget=True,
                 custom_env={**env, 'BDMC_REFERENCE_MANAGED_GROUP': '1'})
        report['source_after'] = source_identity()
        if report['source_after'] != report['source_before']:
            raise RuntimeError('production source identity changed during validation; rerun a stable edition')
        report['status'] = 'passed'
    except Exception as exc:
        report.update(status='failed', error=f'{type(exc).__name__}: {exc}')
        print(report['error'], file=sys.stderr, flush=True)
    finally:
        report['known_limitations'] = [item for entry in report['steps']
                                       for item in entry.get('known_limitations', [])]
        report['completed_at'] = datetime.now(timezone.utc).isoformat()
        report['storage_after'] = budget.check()
        save()
        print('Report:', report_path, flush=True)
    return 0 if report['status'] == 'passed' else 1


if __name__ == '__main__':
    raise SystemExit(main())
