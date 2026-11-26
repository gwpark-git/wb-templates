#!/usr/bin/env python3
"""Bounded, verified timing of legacy repairs, CPU paths, observers and HI.

No historical-physics or GPU speed claim follows from these shared-host runs.
Keep raw samples, inputs, outputs and immutable executable/module identities.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import re
import statistics
import sys
import tempfile
import time

CODE = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def summary(values):
    return dict(samples_seconds=values, median_seconds=statistics.median(values),
                min_seconds=min(values), max_seconds=max(values))


def clocked(operation):
    start = time.perf_counter()
    value = operation()
    return value, time.perf_counter() - start


def measure(operation, repeats):
    operation()  # one unmeasured warmup; no retry on error
    return summary([clocked(operation)[1] for _ in range(repeats)])


def load_helper(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--before', type=Path, required=True)
    parser.add_argument('--after', type=Path, required=True)
    parser.add_argument('--package', type=Path, required=True, help='fresh staged bdmc package directory')
    parser.add_argument('--repeats', type=int, choices=range(3, 10), default=5)
    parser.add_argument('--threads', type=int, choices=range(2, 9), default=4)
    parser.add_argument('--run-dir', type=Path, default=CODE / 'runs/repair-analysis-hi-20261111/performance')
    parser.add_argument('--budget', type=Path, default=CODE / 'DATA_BUDGET.json')
    args = parser.parse_args(argv)
    args.before = args.before.resolve(strict=True)
    args.after = args.after.resolve(strict=True)
    args.package = args.package.resolve(strict=True)
    # Editable install finders must not redirect this benchmark to stale code.
    sys.path.insert(0, str(args.package.parent))
    import importlib.machinery
    class Finder:
        def find_spec(self, fullname, path=None, target=None):
            if fullname == 'bdmc':
                return importlib.machinery.PathFinder.find_spec(fullname, [str(args.package.parent)])
            if fullname.startswith('bdmc.'):
                return importlib.machinery.PathFinder.find_spec(fullname, path)
    sys.meta_path.insert(0, Finder())
    import numpy as np
    import bdmc
    from bdmc import _reference, _analysis, _hydrodynamics
    from bdmc.reference import RepulsiveSpec, RepulsiveSimulation
    from bdmc.analysis import pair_correlation, directional_pair_correlation, coordination_number, autocorrelation
    from bdmc.hydrodynamics import HydrodynamicSpec, HydrodynamicSimulation, StokesianSpec, solve_stokesian
    from bdmc.storage import DataBudget
    if Path(bdmc.__file__).resolve().parent != args.package:
        raise RuntimeError('benchmark imported the wrong package')
    modules = {name: Path(module.__file__).resolve() for name, module in
               (('_reference', _reference), ('_analysis', _analysis), ('_hydrodynamics', _hydrodynamics))}
    if any(path.parent != args.package for path in modules.values()):
        raise RuntimeError('benchmark imported a stale native extension')
    suite = load_helper('_timing_reference_suite', CODE / 'tests/reference_suite.py')
    validator = load_helper('_timing_validator', CODE / 'tools/validate.py')
    budget = DataBudget(args.budget)
    with budget.writer(args.run_dir, 192_000_000):
        args.run_dir.mkdir(parents=True, exist_ok=True)
        run = Path(tempfile.mkdtemp(prefix='run-', dir=args.run_dir))
        report = dict(schema='bdmc.workflow-performance.v1', status='running',
            created_at=datetime.now(timezone.utc).isoformat(),
            source_before=validator.source_identity(),
            executable_identity={name: dict(path=str(path), sha256=digest(path))
                                 for name, path in (('before',args.before),('after',args.after))},
            module_identity={name: dict(path=str(path), sha256=digest(path)) for name,path in modules.items()},
            host=dict(platform=platform.platform(), cpu_count=os.cpu_count(),
                      affinity=sorted(os.sched_getaffinity(0)) if hasattr(os,'sched_getaffinity') else None,
                      cpu=Path('/proc/cpuinfo').read_text().split('model name',1)[-1].splitlines()[0].strip()
                          if Path('/proc/cpuinfo').exists() else platform.processor(),
                      load_before=os.getloadavg() if hasattr(os,'getloadavg') else None),
            protocol=dict(repeats=args.repeats, warmups=1, threads=args.threads,
                          legacy_seeds=[100,150], modern_seed=100, analysis_seed=150, hi_seed=810729,
                          note='Fixed workloads/seeds. No retries, reseeding, clipping or automatic backend fallback.'),
            legacy={}, modern_cpu={}, analysis={}, hydrodynamics={},
            limitations=['Shared host; concurrent HI session may contend for CPU.',
                         'Legacy repairs intentionally alter buggy physics/scheduling; cost comparison is not equivalence.',
                         'HI and free-draining engines represent different dynamics and RNGs; cost only.',
                         'No CUDA/Apple device timings; no complete modern network controller timing.'])
        def save():
            (run / 'report.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
        try:
            save()
            # Equal input, N, dt, seed, threads and recording within each before/after pair.
            cases = {name: {**overrides} for name,overrides in suite.CASES.items()
                     if name != 'equilibration'}
            cases['network_record_dense'] = {**suite.CASES['network'], 'N_skip_file':'5', 'N_skip_ener':'5'}
            cases['network_bd4'] = {**suite.CASES['network'], 'N_THREADS_BD':str(args.threads)}
            cases['network_ss4'] = {**suite.CASES['network'], 'N_THREADS_SS':str(args.threads)}
            # Before SS>1 is known unsafe; retain cost/observations, never bless its topology.
            cases.pop('network_ss4')
            for name,overrides in cases.items():
                effective = {'Np':'128','Nt':'301','N_skip_file':'30','N_skip_ener':'30',
                             'N_chains_per_particle':'3', **overrides}
                # All benchmark cases use the common bounded size/time, including brownian.
                effective.update(Np='128',Nt='301')
                text = suite.input_text(effective)
                entry = dict(input_sha256=hashlib.sha256(text.encode()).hexdigest(),
                             parameters=suite.parse_input(text), before={}, after={})
                samples = {key:[] for key in ('before','after')}
                compute_samples = {key:[] for key in samples}
                paths = {key:[] for key in samples}
                for repetition in range(-1,args.repeats):
                    # Alternate pair order to reduce a simple warm-cache/order bias.
                    order = ('before','after') if repetition%2 else ('after','before')
                    for label in order:
                        exe = getattr(args,label)
                        case_name = f'{name}-{label}-{repetition+1}'
                        (directory,validated),elapsed = clocked(lambda: suite.run_case(
                            exe,run,case_name,text,budget,
                            'legacy-v1' if label=='before' else 'repaired-v2'))
                        if repetition>=0:
                            samples[label].append(elapsed)
                            log = (directory/'stdout.log').read_text()
                            native = re.findall(r'Total simulation time\s*=\s*([\deE+.-]+)', log)
                            if len(native)!=1: raise RuntimeError('missing native simulation timer')
                            compute_samples[label].append(float(native[0]))
                            paths[label].append(directory)
                        entry[label]['validation']=validated
                for label in samples:
                    entry[label].update(summary(samples[label]))
                    entry[label]['native_simulation']=summary(compute_samples[label])
                    entry[label]['workflow_overhead']=summary([
                        full-native for full,native in zip(samples[label],compute_samples[label])])
                    entry[label]['timer_scope']='End-to-end includes process/setup, I/O and independent Python invariant validation; native_simulation is the rounded internal loop timer.'
                    log = (paths[label][-1]/'stdout.log').read_text()
                    phases = re.findall(r'time consuming: LV = ([\deE+.-]+).*?SS = ([\deE+.-]+).*?AN = ([\deE+.-]+).*?FILE = ([\deE+.-]+).*?DIST = ([\deE+.-]+)',log)
                    entry[label]['last_printed_cumulative_phases_seconds'] = (
                        dict(zip(('BD','topology','analysis','file','distance'),map(float,phases[-1]))) if phases else None)
                    entry[label]['phase_caveat']='Rounded cumulative instrumentation at final recorded frame; last unrecorded tail excluded.'
                    entry[label]['run_directories']=[str(p) for p in paths[label]]
                try:
                    entry['before_after_agreement']=suite.compare_runs(paths['before'][0],paths['after'][0])
                except AssertionError as exc:
                    entry['before_after_agreement']=dict(status='intentionally_not_identical', first_difference=str(exc))
                entry['after_over_before_time']=entry['after']['median_seconds']/entry['before']['median_seconds']
                entry['after_over_before_native_time']=(entry['after']['native_simulation']['median_seconds']/
                                                        entry['before']['native_simulation']['median_seconds'])
                report['legacy'][name]=entry
                print('TIMED legacy',name,entry['after_over_before_time'],flush=True)
                save()

            for n in (128,512,1024):
                length=(n/.4)**(1/3)
                q=np.random.default_rng(100).uniform(0,length,(n,3))
                steps=20
                def simulation(backend,search,threads):
                    return RepulsiveSimulation(RepulsiveSpec(n_particles=n,box_length=length,dt=1e-5,
                        diffusion=.1,epsilon=1.,sigma=1.25,seed=100,
                        force_search=search,execution_backend=backend,backend_threads=threads),q)
                settings={'serial_cell':('serial','cell_list',1),
                          'serial_brute':('serial','brute_force',1),
                          'openmp_brute':('openmp','brute_force',args.threads)}
                frames={}
                for label,setting in settings.items():
                    sim=simulation(*setting);sim.advance(steps);frames[label]=sim.observe()
                errors={}
                for label,frame in frames.items():
                    err=float(np.max(np.abs(frame['positions']-frames['serial_cell']['positions'])))
                    np.testing.assert_allclose(frame['positions'],frames['serial_cell']['positions'],rtol=2e-12,atol=2e-12)
                    np.testing.assert_allclose(frame['forces'],frames['serial_cell']['forces'],rtol=2e-11,atol=2e-11)
                    errors[label]=err
                entry=dict(n_particles=n,steps=steps,box_length=length,max_position_errors=errors,
                           interaction_check=float(np.linalg.norm(frames['serial_cell']['forces'])))
                if entry['interaction_check']==0: raise RuntimeError('timing fixture has no active pair forces')
                for label,setting in settings.items():
                    operation=lambda setting=setting: (lambda sim:(sim.advance(steps),sim.observe()))(simulation(*setting))
                    entry[label]=measure(operation,args.repeats)
                entry['scope']='Construction + C++ advance + copied Python observation; same noise/model/initial state.'
                report['modern_cpu'][str(n)]=entry
                print('TIMED modern',n,flush=True);save()

            for n in (256,2048):
                length=(n/.4)**(1/3)
                q=np.random.default_rng(150).uniform(0,length,(3,n,3))
                edges=np.linspace(0,1.25,21)
                operations={backend:(lambda backend=backend:pair_correlation(q,[length]*3,edges,backend=backend))
                            for backend in ('cell_list','brute_force')}
                values={label:op() for label,op in operations.items()}
                np.testing.assert_array_equal(values['cell_list']['counts'],values['brute_force']['counts'])
                np.testing.assert_allclose(coordination_number(values['cell_list']),
                                           coordination_number(values['cell_list'],source='g'),rtol=1e-14,atol=1e-14)
                directional=directional_pair_correlation(q,[length]*3,edges,[-1,0,1],np.linspace(0,2*np.pi,5))
                np.testing.assert_array_equal(np.sum(directional['counts'],axis=(1,2)),values['cell_list']['counts'])
                report['analysis'][str(n)]=dict(n_frames=3,n_particles=n,ordered_counts=int(values['cell_list']['counts'].sum()),
                    candidate_pairs={label:int(v['candidate_pairs']) for label,v in values.items()},
                    timings={label:measure(op,args.repeats) for label,op in operations.items()},
                    directional=measure(lambda:directional_pair_correlation(q,[length]*3,edges,[-1,0,1],np.linspace(0,2*np.pi,5)),args.repeats),
                    checks='Exact brute/cell counts; directional marginal; direct/g coordination integral.')
                save()
            series=np.random.default_rng(150).normal(size=(4096,6))
            report['analysis']['autocorrelation']=measure(lambda:autocorrelation(series,max_lag=128,center=True),args.repeats)

            pi=np.pi
            for n in (2,8,32,64):
                side=int(np.ceil(n**(1/3)))
                q=np.array([[i%side,(i//side)%side,i//(side*side)] for i in range(n)],dtype=float)*1.1
                s=HydrodynamicSpec(n_particles=n,radius=.2,viscosity=1/(6*pi*.2),dt=1e-5,
                                  thermal_energy=.1,seed=810729,mobility='rpy')
                def hi():
                    sim=HydrodynamicSimulation(s,q);sim.advance(10);return sim.observe()
                first=hi();second=hi()
                np.testing.assert_array_equal(first['positions'],second['positions'])
                # Analytic self-mobility path is an explicit cost reference,
                # not a substitute for (or output-equivalent to) the HI model.
                noise=np.random.default_rng(810729).normal(size=(10,n,3))
                def fd():
                    x=q.copy()
                    for z in noise: x+=np.sqrt(2*.1*1e-5)*z
                    return x
                entry=dict(n_particles=n,steps=10,hi=measure(hi,args.repeats),
                           free_draining_numpy=measure(fd,args.repeats),
                           note='HI assembles/factorizes mobility and validates candidate PD every step; FD uses precomputed independent NumPy noise. Different work/physics/RNG: cost comparison only.')
                if n<=32:
                    forces=np.zeros((n,3));forces[:,0]=.1
                    spec=StokesianSpec(radius=.5,viscosity=1.,lubrication=True)
                    solve=lambda:solve_stokesian(q,spec,forces=forces)
                    solved=solve()
                    if solved['relative_residual']>1e-8: raise RuntimeError('SD residual exceeds acceptance')
                    entry['stokesian']=measure(solve,args.repeats)
                    entry['stokesian_residual']=solved['relative_residual']
                    entry['stokesian_spec']=dict(radius=.5,viscosity=1.,lubrication=True,
                                                nearest_distance=1.1,gap_over_radius=.2,
                                                note='Nearest pairs inside active lubrication cutoff 3.99a.')
                report['hydrodynamics'][str(n)]=entry
                print('TIMED HI',n,flush=True);save()
            report['source_after']=validator.source_identity()
            if report['source_before']!=report['source_after']: raise RuntimeError('source edition changed during benchmark')
            for label,path in (('before',args.before),('after',args.after)):
                if digest(path)!=report['executable_identity'][label]['sha256']: raise RuntimeError('timed executable changed')
            for name,path in modules.items():
                if digest(path)!=report['module_identity'][name]['sha256']: raise RuntimeError('timed module changed')
            report['status']='passed'
        except Exception as exc:
            report.update(status='failed',error=f'{type(exc).__name__}: {exc}')
            raise
        finally:
            report['completed_at']=datetime.now(timezone.utc).isoformat()
            report['host']['load_after']=os.getloadavg() if hasattr(os,'getloadavg') else None
            report['storage_after']=budget.check()
            save()
            print('Report:',run/'report.json',flush=True)
    return 0


if __name__=='__main__':
    raise SystemExit(main())
