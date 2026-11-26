#!/usr/bin/env python3
"""Bounded independent theory comparisons; never refreezes legacy references.

Acceptance is fixed before generation: each Gaussian moment/independent-replica
estimate must lie within six estimated standard errors of its declared oracle.
Both time steps must pass; this is a finite pilot, not asymptotic convergence.
Failures are retained, without automatic reseeding or dropping trajectories.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.integrate import quad, solve_ivp

from bdmc import _reference
from bdmc.reference import DumbbellSimulation, DumbbellSpec, RepulsiveSimulation, RepulsiveSpec
from bdmc.storage import DataBudget

CODE = Path(__file__).resolve().parents[1]


def compare(value, expected, sem, **labels):
    value, expected, sem = map(float, (value, expected, sem))
    limit = 6 * sem + 1e-12
    return dict(labels, observed=value, expected=expected, standard_error=sem,
                acceptance_absolute_error=limit, passed=bool(
                    np.isfinite([value, expected, sem]).all() and abs(value - expected) <= limit))


def moment_reference(protocol, dt, steps):
    """Independent moment propagation; continuum solved as a separate ODE."""
    m = np.diag([.25, .5, .75]) if protocol == 'step_shear' else np.zeros((3, 3))
    if protocol == 'step_shear':
        f = np.eye(3); f[0, 1] = .7
        m = f @ m @ f.T
    initial = m.copy()

    def rate(t):
        return {'simple_shear_positive': 1., 'simple_shear_negative': -1.,
                'step_shear': 0., 'oscillatory_shear': .6 * 2.3 * np.cos(2.3 * t)}[protocol]

    for n in range(steps):
        a = (1 - 2 * dt) * np.eye(3); a[0, 1] = dt * rate(n * dt)
        m = a @ m @ a.T + 2 * dt * np.eye(3)  # mu=H=1,D=.5

    def rhs(t, flat):
        tensor = flat.reshape(3, 3)
        k = np.zeros((3, 3)); k[0, 1] = rate(t)
        return (k @ tensor + tensor @ k.T - 4 * tensor + 2 * np.eye(3)).ravel()

    continuum = solve_ivp(rhs, (0., steps * dt), initial.ravel(),
                          rtol=1e-11, atol=1e-13).y[:, -1].reshape(3, 3)
    return m, continuum


def hookean_checks(seeds):
    results = []
    for protocol in ('simple_shear_positive', 'simple_shear_negative', 'step_shear', 'oscillatory_shear'):
        for dt in (.005, .0025):
            steps = round(1.5 / dt)
            expected, continuum = moment_reference(protocol, dt, steps)
            for seed in seeds:
                initial = None
                if protocol == 'step_shear':
                    # Exact zero mean and prescribed raw second moment; no random
                    # initialization variance is hidden in the reference.
                    initial = np.tile(np.vstack([np.diag(np.sqrt([.75, 1.5, 2.25])),
                                                  -np.diag(np.sqrt([.75, 1.5, 2.25]))]), (1024, 1))
                spec = DumbbellSpec(n_dumbbells=6144, dt=dt, mobility=1., diffusion=.5,
                                    spring_constant=1., protocol=(
                                        'simple_shear' if protocol.startswith('simple_shear') else protocol),
                                    shear_rate=(-1. if protocol.endswith('negative') else 1.)
                                    if protocol.startswith('simple_shear') else 0.,
                                    step_strain=.7 if protocol == 'step_shear' else 0.,
                                    strain_amplitude=.6 if protocol == 'oscillatory_shear' else 0.,
                                    angular_frequency=2.3, seed=seed)
                sim = DumbbellSimulation(spec) if initial is None else DumbbellSimulation(spec, initial)
                sim.advance(steps)
                frame = sim.observe(); q = frame['connectors']
                for i, j in ((0, 0), (0, 1), (1, 1), (2, 2)):
                    samples = q[:, i] * q[:, j]
                    results.append(compare(samples.mean(), expected[i, j],
                                           samples.std(ddof=1) / np.sqrt(len(samples)),
                                           model='hookean', protocol=protocol, dt=dt, seed=seed,
                                           specification=spec.to_dict(),
                                           component=f'{i}{j}', continuum=float(continuum[i, j]),
                                           discretization_bias=float(expected[i, j] - continuum[i, j])))
    return results


def fene_checks(seeds):
    results = []
    # b=H Q0²/kT=18; d=3, exact equilibrium <q²>=27/23.
    expected_q2 = 27 / 23
    for dt in (.001, .0005):
        for seed in seeds:
            spec = DumbbellSpec(n_dumbbells=4096, dt=dt, mobility=1., diffusion=.5,
                                spring_constant=1., spring_law='fene', q0=3., seed=seed)
            sim = DumbbellSimulation(spec); sim.advance(round(3. / dt))
            frame = sim.observe(); q = frame['connectors']; f = frame['spring_forces']
            q2 = (q * q).sum(axis=1)
            results.append(compare(q2.mean(), expected_q2, q2.std(ddof=1) / np.sqrt(len(q2)),
                                   model='fene', dt=dt, seed=seed, component='mean_q_squared',
                                   specification=spec.to_dict(),
                                   horizon=3., max_extension=float(np.sqrt(q2.max())), q0=3.,
                                   oracle='continuum equilibrium; finite-time/dt pilot'))
            for i in range(3):
                samples = q[:, i] * f[:, i]
                results.append(compare(samples.mean(), .5, samples.std(ddof=1) / np.sqrt(len(q)),
                                       model='fene', dt=dt, seed=seed, component=f'qF_{i}{i}',
                                       oracle='equilibrium integration-by-parts <q_i F_i>=kT'))
    return results


def repulsive_checks(seeds):
    results = []
    sigma, epsilon, kt, radius, box = 1., 5., 1., 1.8, 4.

    def weight(r):
        x = r / sigma
        u = epsilon * (2 - 3 * x + x**3) / 3 if x < 1 else 0.
        return r * r * np.exp(-u / kt)

    norm = quad(weight, 0., radius, points=[sigma], epsabs=1e-12)[0]
    expected = [quad(weight, 0., sigma, epsabs=1e-12)[0] / norm,
                quad(lambda r: r*r*weight(r), 0., radius, points=[sigma], epsabs=1e-12)[0] / norm]
    for dt in (.002, .001):
        replicates = []
        for seed in seeds:
            spec = RepulsiveSpec(n_particles=2, box_length=box, dt=dt, mobility=1., diffusion=1.,
                                 epsilon=epsilon, sigma=sigma, seed=seed, force_search='brute_force')
            sim = RepulsiveSimulation(spec, np.array([[.5, .5, .5], [2., 2., 2.]]))
            sim.advance(round(20. / dt))
            batch = sim.sample(round(120. / dt), stride=round(.04 / dt), include_initial=False)
            delta = batch['positions'][:, 1] - batch['positions'][:, 0]
            delta -= box * np.floor(delta / box + .5)
            r2 = (delta * delta).sum(axis=1); inside = r2 < radius**2
            replicates.append([np.mean(r2 < sigma**2), np.mean(r2 * inside), np.mean(inside)])
        replicas = np.array(replicates)
        for index, label in enumerate(('contact_probability_given_r_lt_R', 'mean_r_squared_given_r_lt_R')):
            ratio = replicas[:, index].mean() / replicas[:, 2].mean()
            # Independent trajectories are the statistical units. Time frames
            # inside each trajectory are deliberately not treated as iid.
            sem = (replicas[:, index] - ratio * replicas[:, 2]).std(ddof=1)
            sem /= np.sqrt(len(seeds)) * replicas[:, 2].mean()
            results.append(compare(ratio, expected[index], sem, model='repulsive_pair', dt=dt,
                                   component=label, radius=radius, epsilon=epsilon, thermal_energy=kt,
                                   replica_joint_means=replicates, burn_time=20., sample_time=120.,
                                   oracle='3D conditional radial Boltzmann density; R<L/2'))
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--budget', type=Path, default=CODE / 'DATA_BUDGET.json')
    args = parser.parse_args()
    budget = DataBudget(args.budget)
    start = time.monotonic()
    manifest_path = CODE / 'examples/two_reference_physics_seeds.json'
    manifest = json.loads(manifest_path.read_text())
    report = {'schema': 'bdmc.two-reference-physics.v1',
              'created_at': datetime.now(timezone.utc).isoformat(), 'status': 'running',
              'seed_manifest': manifest, 'checks': [],
              'seed_manifest_sha256': hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
              'native_extension': str(_reference.__file__),
              'native_extension_sha256': hashlib.sha256(Path(_reference.__file__).read_bytes()).hexdigest(),
              'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'acceptance': '6 estimated standard errors +1e-12 per preregistered check; both dt values',
              'limits': 'Finite independent ensembles/replicas; not proof of equilibration, asymptotic convergence or many-particle rheology. No direct Bird page verification.'}
    with budget.writer(args.output, 4_000_000):
        args.output.mkdir(parents=True, exist_ok=False)
        try:
            for key, function in [('hookean', hookean_checks), ('fene', fene_checks), ('repulsive', repulsive_checks)]:
                report['active_family'] = key
                report['checks'].extend(function(manifest[key]))
            report['status'] = 'passed' if all(c['passed'] for c in report['checks']) else 'failed'
            report.pop('active_family', None)
        except Exception as exc:
            report.update(status='failed', error=f'{type(exc).__name__}: {exc}')
        report['elapsed_seconds'] = time.monotonic() - start
        report['storage_after'] = budget.check()
        (args.output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(args.output / 'report.json', report['status'], flush=True)
    return 0 if report['status'] == 'passed' else 1


if __name__ == '__main__':
    raise SystemExit(main())
