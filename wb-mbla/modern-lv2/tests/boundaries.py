#!/usr/bin/env python3
"""Behavioral input/output contracts, including successful-open/failed-write cases.

This stdlib test never changes frozen data. POSIX RLIMIT_FSIZE is used only in
isolated child processes; unavailable POSIX checks are explicitly NOT_RUN.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]
# Preserve the signal disposition in the actual solver child. Python otherwise
# resets SIGXFSZ, which would test signal termination instead of EFBIG handling.
LIMITED = '''import resource,signal,subprocess,sys
resource.setrlimit(resource.RLIMIT_FSIZE,(1024,1024))
signal.signal(signal.SIGXFSZ,signal.SIG_IGN)
sys.exit(subprocess.run(sys.argv[1:],restore_signals=False).returncode)
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--work', type=Path, required=True)
    parser.add_argument('--section', choices=['input', 'output', 'all'], default='all')
    args = parser.parse_args()
    args.exe = args.exe.resolve()
    args.work = args.work.resolve()
    args.work.mkdir(parents=True, exist_ok=True)
    results = []
    config = (ROOT/'examples/noncake.inp').read_text()
    # Constant closures and mild CP give a fast supported input; no production
    # output is used as a physical expected solution in these interface tests.
    config = config.replace('ref_DTP=5000', 'ref_DTP=100')
    config = config.replace('phiw_init=0.01', 'phiw_init=0.001')
    config = config.replace('fcn_Dc_div_D0=SPHS', 'fcn_Dc_div_D0=CT')
    config = config.replace('fcn_eta_div_eta0=SPHS', 'fcn_eta_div_eta0=CT')

    def run(name, text, limit=False):
        directory = args.work/name
        directory.mkdir(exist_ok=True)
        (directory/'case.inp').write_text(text)
        argv = [str(args.exe), 'case.inp']
        if limit:
            argv = [sys.executable, '-c', LIMITED, *argv]
        start = time.monotonic()
        proc = subprocess.run(argv, cwd=directory, capture_output=True, text=True, timeout=35)
        (directory/'stdout.log').write_text(proc.stdout)
        (directory/'stderr.log').write_text(proc.stderr)
        (directory/'command.json').write_text(json.dumps({
            'argv': argv, 'cwd': str(directory), 'returncode': proc.returncode,
            'seconds': time.monotonic()-start, 'limit_bytes': 1024 if limit else None,
            'exe_sha256': hashlib.sha256(args.exe.read_bytes()).hexdigest()}, indent=2)+'\n')
        return directory, proc

    if args.section in ('input', 'all'):
        baseline, proc = run('integer-decimal', config)
        baseline_ok = proc.returncode == 0 and 'MBLA_STATUS=converged' in proc.stdout
        results.append({'behavior': 'typed-count-baseline', 'status': 'PASS' if baseline_ok else 'FAIL'})
        for name, spelling in [('scientific', '2.4e3'), ('hexadecimal', '0x960')]:
            directory, proc = run('integer-'+name, config.replace('N_iter=2400', 'N_iter='+spelling))
            consumed_line = next((x for x in proc.stdout.splitlines() if 'N_iter, weight,' in x), '')
            consumed = consumed_line.split(':', 1)[-1].split()[0] if consumed_line else None
            equal = baseline_ok and (directory/'result.dat').exists() and (
                (directory/'result.dat').read_bytes() == (baseline/'result.dat').read_bytes())
            ok = proc.returncode == 0 and consumed == '2400' and equal
            results.append({'behavior': 'typed-count-'+name, 'input': spelling,
                            'consumed_count': consumed, 'same_printed_output_as_decimal': equal,
                            'returncode': proc.returncode, 'status': 'PASS' if ok else 'FAIL'})

        directory = Path(tempfile.mkdtemp(prefix='ordered-batch-', dir=args.work))
        first = config.replace('fout=result.dat', 'fout=first.dat')
        second = config.replace('fout=result.dat', 'fout=second.dat')
        second = second.replace('CONT=FALSE', 'CONT=TRUE\nFN_CONT=first.dat')
        (directory / 'first.inp').write_text(first)
        (directory / 'second.inp').write_text(second)
        batch = subprocess.run([str(args.exe), 'first.inp', 'second.inp'], cwd=directory,
                               capture_output=True, text=True, timeout=35)
        ordered = (batch.returncode == 0 and batch.stdout.count('MBLA_STATUS=converged') == 2
                   and all((directory / f'{name}.dat').is_file() for name in ('first', 'second')))
        results.append({'behavior': 'ordered-batch-continuation', 'status': 'PASS' if ordered else 'FAIL'})

        # Case three retains its diagnostic table. Case four must never run.
        (directory / 'third.inp').write_text(config.replace('N_iter=2400', 'N_iter=1')
                                             .replace('fout=result.dat', 'fout=third.dat'))
        (directory / 'fourth.inp').write_text(config.replace('fout=result.dat', 'fout=fourth.dat'))
        stopped = subprocess.run([str(args.exe), 'first.inp', 'second.inp', 'third.inp', 'fourth.inp'],
                                 cwd=directory, capture_output=True, text=True, timeout=35)
        limited = (stopped.returncode == 3 and stopped.stdout.count('MBLA_STATUS=converged') == 2
                   and 'MBLA_STATUS=iteration_limit' in stopped.stdout
                   and (directory / 'third.dat').is_file() and not (directory / 'fourth.dat').exists())
        results.append({'behavior': 'batch-stop-on-iteration-limit', 'status': 'PASS' if limited else 'FAIL'})

    if args.section in ('output', 'all'):
        try:
            import resource
            import signal
            available = hasattr(resource, 'RLIMIT_FSIZE') and hasattr(signal, 'SIGXFSZ')
        except ImportError:
            available = False
        if not available:
            results.append({'behavior': 'late-write-failure', 'status': 'NOT_RUN',
                            'reason': 'POSIX RLIMIT_FSIZE/SIGXFSZ unavailable'})
        else:
            for name, target, text in [
                ('primary', 'result.dat', config),
                ('outlet-sidecar', 'result.dat.Rout', config.replace('DETAILED_REPORT=FALSE', 'DETAILED_REPORT=TRUE')),
                # Twelve transverse rows keep Rout below1024 bytes; the25-row
                # Fval exceeds the same bound, exercising its own write path.
                ('integral-sidecar', 'result.dat.Fval', config.replace('Nr=120', 'Nr=12').replace('DETAILED_REPORT=FALSE', 'DETAILED_REPORT=TRUE')),
            ]:
                baseline, normal = run(name+'-baseline', text)
                baseline_ok = normal.returncode == 0 and (baseline/target).exists() and (baseline/target).stat().st_size > 1024
                directory, limited = run(name+'-limited', text, limit=True)
                ok = baseline_ok and limited.returncode == 2 and target in limited.stderr and 'MBLA_STATUS=converged' not in limited.stdout
                results.append({'behavior': 'late-write-'+name, 'target': target,
                    'baseline_ok': baseline_ok, 'baseline_bytes': (baseline/target).stat().st_size if (baseline/target).exists() else None,
                    'limited_bytes': (directory/target).stat().st_size if (directory/target).exists() else None,
                    'returncode': limited.returncode, 'diagnostic': limited.stderr.strip(),
                    'status': 'PASS' if ok else 'FAIL'})
    (args.work/'RESULT.json').write_text(json.dumps({'cases': results}, indent=2)+'\n')
    print(json.dumps(results, indent=2))
    return 1 if any(r['status'] == 'FAIL' for r in results) else 0


if __name__ == '__main__':
    sys.exit(main())
