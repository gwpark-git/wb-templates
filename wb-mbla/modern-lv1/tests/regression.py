#!/usr/bin/env python3
"""Bounded stdlib integration checks; never regenerate the frozen fixtures."""
import argparse
import json
import math
from pathlib import Path
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]


def read_table(path):
    rows = [[float(v) for v in line.split()] for line in path.read_text().splitlines() if line.strip()]
    if not rows or any(len(row) != 10 or not all(map(math.isfinite, row)) for row in rows):
        raise AssertionError(f"invalid ten-column finite output: {path}")
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--exe', required=True, type=Path)
    p.add_argument('--work', required=True, type=Path)
    a = p.parse_args()
    a.exe = a.exe.resolve()
    a.work = a.work.resolve()
    a.work.mkdir(parents=True, exist_ok=True)
    results = []

    def run(name, config, expected=0):
        directory = a.work / name
        directory.mkdir(exist_ok=True)
        (directory / 'case.inp').write_text(config)
        start = time.monotonic()
        result = subprocess.run([str(a.exe), 'case.inp'], cwd=directory,
                                capture_output=True, text=True, timeout=35)
        (directory / 'stdout.log').write_text(result.stdout + result.stderr)
        assert result.returncode == expected, (name, result.returncode, result.stdout[-500:], result.stderr)
        return directory, time.monotonic() - start

    for name in ('noncake', 'cake', 'cake-mild'):
        config = (ROOT / 'examples' / f'{name}.inp').read_text()
        directory, seconds = run(name, config)
        actual = read_table(directory / 'result.dat')
        expected = read_table(ROOT / 'tests' / 'fixtures' / f'{name}-strict.dat')
        assert len(actual) == len(expected) == 25
        max_scaled_error = 0.0
        exact = actual == expected
        for row, reference in zip(actual, expected):
            for value, target in zip(row, reference):
                # The inherited printer stores six decimal digits after the
                # decimal point in scientific notation. Allow two last-digit
                # rounding units across toolchains, never a historical tolerance.
                allowance = 1e-12 + 2e-6 * abs(target)
                max_scaled_error = max(max_scaled_error, abs(value-target) / allowance)
                assert abs(value-target) <= allowance, (name, value, target)
        assert abs(actual[0][0]) < 1e-15 and abs(actual[-1][0]-0.5) < 1e-12
        assert all(actual[i][0] < actual[i+1][0] for i in range(24))
        assert actual[0][1] == 0.001
        assert all(row[7] >= 0 for row in actual)
        if name == 'noncake':
            assert all(row[7] == 0 for row in actual), 'noncake fixture must remain cake-free'
            assert max(row[1] for row in actual) < 0.494
        else:
            assert actual[-1][7] > 0, 'cake fixture must exercise resistance'
            assert abs(actual[-1][1]-0.494) < 1e-12
        if name == 'cake-mild':
            assert max(abs(row[6]/actual[0][6]-1) for row in actual) < 0.01
        assert 'MBLA_STATUS=converged' in (directory / 'stdout.log').read_text()
        results.append({'case': name, 'status': 'passed', 'seconds': seconds,
                        'printed_values_exact': exact, 'max_scaled_error': max_scaled_error,
                        'wall_phi_outlet': actual[-1][1],
                        'cake_resistance_ratio_outlet': actual[-1][7] / 0.0005,
                        'particle_flux_relative_drift': max(abs(r[6]/actual[0][6]-1) for r in actual)})

    # Constant table values must reproduce analytic CT choices; this exercises
    # valid spline allocation, units and repeated ownership/destruction.
    base_ct = (ROOT / 'examples' / 'noncake.inp').read_text().replace('ref_DTP=5000','ref_DTP=1000')
    base_ct = base_ct.replace('fcn_Dc_div_D0=SPHS','fcn_Dc_div_D0=CT').replace('fcn_eta_div_eta0=SPHS','fcn_eta_div_eta0=CT')
    inline_ct, _ = run('constant_analytic', base_ct)
    table_dir = a.work / 'constant_tables'
    table_dir.mkdir(exist_ok=True)
    (table_dir / 'unity.dat').write_text('0 1\n0.25 1\n0.5 1\n')
    table_ct, _ = run('constant_tables', base_ct.replace('fcn_Dc_div_D0=CT','fcn_Dc_div_D0=unity.dat').replace('fcn_eta_div_eta0=CT','fcn_eta_div_eta0=unity.dat'))
    assert read_table(inline_ct/'result.dat') == read_table(table_ct/'result.dat')
    details, _ = run('detailed_report', base_ct.replace('DETAILED_REPORT=FALSE','DETAILED_REPORT=TRUE'))
    assert len((details/'result.dat.Rout').read_text().splitlines()) == 120
    assert all(len(line.split()) == 6 for line in (details/'result.dat.Rout').read_text().splitlines())
    assert len((details/'result.dat.Fval').read_text().splitlines()) == 25
    # Inspect actual printed data through a separate formula/interface.
    checker = ROOT / 'tools' / 'inspect_output.py'
    diagnostic = subprocess.check_output([__import__('sys').executable, str(checker),
        str(a.work/'cake-mild/result.dat'),'--radius','0.0005','--permeability','6.7e-10',
        '--csv',str(a.work/'cake-mild/result.csv')], text=True)
    diagnostic = json.loads(diagnostic)
    assert diagnostic['darcy_starling_max_residual_normalized_by_max_wall_velocity'] < 2e-6
    assert diagnostic['max_recorded_particle_flux_relative_drift'] < 0.01

    # Each counterexample targets a plausible former failure: atoi("bad")==0,
    # silent duplicate override, undersized stencils, a missing file/parent,
    # malformed continuation data or GSL's default abort on invalid splines.
    base = (ROOT / 'examples' / 'noncake.inp').read_text()
    failures = {
        'missing_input': base.replace('Nz=25\n', ''),
        'zero_diffusion_scale': base.replace('gamma=1.0', 'gamma=0'),
        'nonfinite': base.replace('weight=0.1', 'weight=nan'),
        'bad_integer': base.replace('Nz=25', 'Nz=25oops'),
        'small_stencil': base.replace('Nz=25', 'Nz=2'),
        'duplicate': base + 'Nz=40\n',
        'typo': base + 'inlet_BCC=velocity\n',
        'unavailable_law': base.replace('fcn_Dc_div_D0=SPHS', 'fcn_Dc_div_D0=absent.dat'),
        'unwritable_output': base.replace('fout=result.dat', 'fout=absent/result.dat'),
        'bad_continuation': base.replace('CONT=FALSE', 'CONT=TRUE\nFN_CONT=short.dat'),
        'bad_property_table': base.replace('fcn_Dc_div_D0=SPHS', 'fcn_Dc_div_D0=bad.dat'),
    }
    for name, config in failures.items():
        directory = a.work / name
        directory.mkdir(exist_ok=True)
        if name == 'bad_continuation':
            (directory / 'short.dat').write_text('0 0.001\n1 0.001\n')
        if name == 'bad_property_table':
            (directory / 'bad.dat').write_text('0 1\n0.2 1\n0.1 1\n0.5 1\n')
        run(name, config, expected=2)
    _, seconds = run('iteration_limit', base.replace('N_iter=2400', 'N_iter=1'), expected=3)
    # The accepted config grammar includes spaces and comments, unlike the old parser.
    _, seconds = run('whitespace', base.replace('N_iter=2400', ' N_iter = 2400 # maximum iterations'))
    summary = {'cases': results, 'invalid_input_cases_rejected': len(failures),
               'iteration_limit_status': 'passed', 'whitespace_status': 'passed',
               'constant_table_equivalence': 'passed', 'detailed_output_schema': 'passed',
               'independent_printed_darcy_diagnostic': diagnostic,
               'limits': 'Current-source fixed regression, not a historical tolerance or a general physical/cake convergence claim.'}
    (a.work / 'RESULT.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
