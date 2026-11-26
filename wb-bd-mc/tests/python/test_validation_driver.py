"""Failures in the test runner must not turn into successful validation reports."""
import importlib.util
import json
import os
from pathlib import Path
import sys
import subprocess
import time

import pytest

_source = Path(__file__).resolve().parents[2] / 'tools/validate.py'
_spec = importlib.util.spec_from_file_location('validation_driver', _source)
driver = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(driver)


def test_zero_exit_with_sanitizer_diagnostic_is_failure(tmp_path):
    result = driver.execute_step(
        [sys.executable, '-c', 'print("runtime error: division by zero")'],
        tmp_path / 'diagnostic.log')
    assert result['returncode'] == 0
    assert result['status'] == 'failed'
    assert 'diagnostic' in result['error']


def test_only_declared_negative_rejection_lines_are_acknowledged(tmp_path):
    expected = 'ERR: deliberately rejected test input'
    for printed, passed in (([expected, expected], True), ([expected], False),
                            ([expected] * 3, False), ([expected, 'ERR: unexpected'], False)):
        result = driver.execute_step([sys.executable, '-c', f'print({chr(10).join(printed)!r})'],
                                      tmp_path / 'expected.log', expected_diagnostics=(expected,) * 2)
        assert (result['status'] == 'passed') == passed


def test_nonzero_exit_and_successful_known_issue_reporting(tmp_path):
    failure = driver.execute_step([sys.executable, '-c', 'raise SystemExit(7)'], tmp_path / 'bad.log')
    assert failure['status'] == 'failed' and failure['returncode'] == 7
    success = driver.execute_step([sys.executable, '-c',
                                  'print("XFAIL: F09 known boundary"); '
                                  'print("XFAIL multi_box_gradient_crossing [F09]")'],
                                  tmp_path / 'ok.log')
    assert success['status'] == 'passed'
    assert success['known_limitations'] == ['XFAIL: F09 known boundary',
                                            'XFAIL multi_box_gradient_crossing [F09]']


@pytest.mark.skipif(os.name != 'posix', reason='POSIX process-group cleanup contract')
def test_timeout_terminates_descendant_before_it_can_write(tmp_path):
    sentinel = tmp_path / 'orphan-write'
    child = f'import time,pathlib; time.sleep(.5); pathlib.Path({str(sentinel)!r}).touch()'
    parent = f'import subprocess,sys,time; subprocess.Popen([sys.executable,"-c",{child!r}]); time.sleep(20)'
    start = time.monotonic()
    result = driver.execute_step([sys.executable, '-c', parent], tmp_path / 'timeout.log', timeout=.15)
    assert result['status'] == 'failed' and 'timeout' in result['error']
    assert time.monotonic() - start < 3
    time.sleep(.6)
    assert not sentinel.exists()


def test_excess_log_is_failure_even_if_process_exits_quickly(tmp_path):
    result = driver.execute_step([sys.executable, '-c', 'print("a"*20000)'],
                                  tmp_path / 'long.log', max_log_bytes=1000)
    assert result['status'] == 'failed' and 'byte limit' in result['error']


def test_junit_summary_keeps_missing_optional_dependencies_visible(tmp_path):
    p = tmp_path / 'checks.xml'
    p.write_text('<testsuites><testsuite tests="12" failures="1" errors="0" skipped="2"/>\n'
                 '<testsuite tests="3" failures="0" errors="1" skipped="0"/></testsuites>')
    assert driver.junit_summary(p) == {'tests': 15, 'failures': 1, 'errors': 1, 'skipped': 2}


def test_import_freshness_guard_survives_python_optimization():
    setup = ('import sys,types\n'
             'p=types.ModuleType("bdmc"); p.__path__=[]; p.__file__="/stale/bdmc/__init__.py"\n'
             'm=types.ModuleType("bdmc._core"); m.__file__="/stale/bdmc/_core.so"\n'
             'p._core=m; sys.modules["bdmc"]=p; sys.modules["bdmc._core"]=m\n')
    process = subprocess.run([sys.executable, '-O', '-c',
                              setup + driver.import_check_script('/expected/bdmc', '/expected/bdmc/_core.so')],
                             capture_output=True, text=True, timeout=5)
    assert process.returncode != 0 and 'wrong Python package imported' in process.stderr


def test_reference_import_guard_rejects_stale_extension_under_optimization():
    setup = ('import sys,types\n'
             'p=types.ModuleType("bdmc"); p.__path__=[]; p.__file__="/expected/bdmc/__init__.py"\n'
             'm=types.ModuleType("bdmc._core"); m.__file__="/expected/bdmc/_core.so"\n'
             'r=types.ModuleType("bdmc._reference"); r.__file__="/stale/_reference.so"\n'
             'p._core=m; p._reference=r\n'
             'sys.modules.update({"bdmc":p,"bdmc._core":m,"bdmc._reference":r})\n')
    process = subprocess.run([sys.executable, '-O', '-c', setup + driver.import_check_script(
        '/expected/bdmc', '/expected/bdmc/_core.so', '/expected/bdmc/_reference.so')],
        capture_output=True, text=True, timeout=5)
    assert process.returncode != 0 and 'wrong reference extension imported' in process.stderr


def test_network_import_guard_rejects_stale_extension_under_optimization():
    setup = ('import sys,types\n'
             'p=types.ModuleType("bdmc"); p.__path__=[]; p.__file__="/expected/bdmc/__init__.py"\n'
             'm=types.ModuleType("bdmc._core"); m.__file__="/expected/bdmc/_core.so"\n'
             'r=types.ModuleType("bdmc._reference"); r.__file__="/expected/bdmc/_reference.so"\n'
             'n=types.ModuleType("bdmc._network"); n.__file__="/stale/_network.so"\n'
             'p._core=m; p._reference=r; p._network=n\n'
             'sys.modules.update({"bdmc":p,"bdmc._core":m,"bdmc._reference":r,"bdmc._network":n})\n')
    process = subprocess.run([sys.executable, '-O', '-c', setup + driver.import_check_script(
        '/expected/bdmc', '/expected/bdmc/_core.so', '/expected/bdmc/_reference.so',
        '/expected/bdmc/_network.so')], capture_output=True, text=True, timeout=5)
    assert process.returncode != 0 and 'wrong network extension imported' in process.stderr


@pytest.mark.parametrize('name', ['_analysis', '_hydrodynamics'])
def test_new_module_import_guards_reject_stale_extensions_under_optimization(name):
    setup = ('import sys,types\n'
             'p=types.ModuleType("bdmc"); p.__path__=[]; p.__file__="/expected/bdmc/__init__.py"\n'
             'm=types.ModuleType("bdmc._core"); m.__file__="/expected/bdmc/_core.so"\n'
             f'n=types.ModuleType("bdmc.{name}"); n.__file__="/stale/{name}.so"\n'
             f'p._core=m; p.{name}=n\n'
             f'sys.modules.update({{"bdmc":p,"bdmc._core":m,"bdmc.{name}":n}})\n')
    process = subprocess.run([sys.executable, '-O', '-c', setup + driver.import_check_script(
        '/expected/bdmc', '/expected/bdmc/_core.so',
        additional_modules={name: '/expected/bdmc/' + name + '.so'})],
        capture_output=True, text=True, timeout=5)
    assert process.returncode != 0 and f'wrong {name} extension imported' in process.stderr


@pytest.mark.parametrize('fault', [None, 'source', 'compiler', 'stage', 'options',
                                  'executable', 'extension', 'missing_module'])
def test_repair_acceptance_admits_only_verified_current_build(tmp_path, fault):
    from types import SimpleNamespace
    source = _source.with_name('run_repair_acceptance.py')
    spec = importlib.util.spec_from_file_location('repair_acceptance_driver', source)
    acceptance = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(acceptance)
    build = tmp_path/'build'
    executable = build/'legacy/stochastic_simulation'
    executable.parent.mkdir(parents=True)
    executable.write_bytes(b'verified executable')
    package = build/'import/bdmc'
    package.mkdir(parents=True)
    entries = {}
    for name in ('_core', '_reference', '_network', '_analysis', '_hydrodynamics'):
        path = package/(name+'.so')
        path.write_bytes(name.encode())
        entries[name] = {'path':str(path), 'sha256':driver.digest(path)}
    stages = [{'name':name, 'status':'passed'} for name in
              ('legacy_build_and_units','modern_configure','modern_build',
               'modern_native','python_import','python_tests')]
    stages[1]['command'] = ['cmake','-DBDMC_ENABLE_OPENMP=ON',
                            '-DBDMC_ENABLE_MPI=ON','-DBDMC_BUILD_LEGACY=ON']
    hashes = {'physical_source.cpp':'frozen source hash'}
    record = {'status':'passed', 'profile':'unit', 'compiler':'/toolchain/g++',
              'source_before':hashes.copy(), 'source_after':hashes.copy(),
              'steps':stages, 'build_directory':str(build),
              'legacy_executable_sha256':driver.digest(executable),
              'staged_package':str(package), 'native_extensions':entries}
    validator = SimpleNamespace(source_identity=lambda:hashes.copy(), digest=driver.digest)
    if fault == 'source': record['source_after'] = {'physical_source.cpp':'old source hash'}
    elif fault == 'compiler': record['compiler'] = '/toolchain/clang++'
    elif fault == 'stage': stages[-1]['status'] = 'failed'
    elif fault == 'options': stages[1]['command'].remove('-DBDMC_ENABLE_MPI=ON')
    elif fault == 'executable': executable.write_bytes(b'changed after verification')
    elif fault == 'extension': Path(entries['_analysis']['path']).write_bytes(b'stale module')
    elif fault == 'missing_module': del entries['_hydrodynamics']
    if fault is None:
        acceptance.validate_unit_report(record,validator,'/toolchain/g++')
    else:
        with pytest.raises(ValueError):
            acceptance.validate_unit_report(record,validator,'/toolchain/g++')
