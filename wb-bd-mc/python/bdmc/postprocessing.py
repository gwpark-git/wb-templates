"""Physically defined replacements for the recovered experimental processors.

Sparse topology, bond multiplicity, density normalization, time origins and
units are explicit. Numerical kernels require only NumPy. Positional RDFs use
the separately tested :mod:`bdmc.analysis`; recorded distance PDFs are not RDFs.
"""
from __future__ import annotations

from collections import Counter
from contextlib import ExitStack
from dataclasses import dataclass
from pathlib import Path
import re
import numpy as np


def _finite(value, name):
    a = np.asarray(value, dtype=float)
    if not np.isfinite(a).all():
        raise ValueError(f'{name} must be finite')
    return a


def _positive(value, name):
    a = _finite(value, name)
    if a.ndim != 0 or a <= 0:
        raise ValueError(f'{name} must be a positive scalar')
    return float(a)


def _nonnegative(value, name):
    a = _finite(value, name)
    if a.ndim != 0 or a < 0:
        raise ValueError(f'{name} must be a nonnegative scalar')
    return float(a)


def _integer(value, name, minimum=0):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value < minimum:
        raise ValueError(f'{name} must be an integer >= {minimum}')
    return int(value)


def _box(box, dimension):
    box = _finite(box, 'box lengths')
    if box.shape != (dimension,) or np.any(box <= 0):
        raise ValueError('one positive box length per dimension required')
    return box


def minimum_image(displacements, box_lengths, *, shear_displacement=0., flow_axis=0, gradient_axis=1):
    """Euclidean nearest lattice image, including rectangular Lees–Edwards cells.

    The gradient lattice vector contains +shear_displacement along flow.
    The search bounds follow from the best candidate norm, so skew cells and
    multiple lattice crossings do not rely on a fixed 27-image assumption.
    """
    d = _finite(displacements, 'displacements')
    if d.ndim < 1 or not 1 <= d.shape[-1] <= 3:
        raise ValueError('last dimension must be 1, 2 or 3')
    box = _box(box_lengths, d.shape[-1])
    shift = float(_finite(shear_displacement, 'shear displacement'))
    out = d.copy()
    if shift == 0:
        return out - box * np.floor(out / box + .5)
    f = _integer(flow_axis, 'flow axis'); g = _integer(gradient_axis, 'gradient axis')
    if f == g or max(f, g) >= len(box):
        raise ValueError('distinct valid shear axes required')
    shift = (shift + box[f] / 2) % box[f] - box[f] / 2
    turns = np.floor(out[..., g] / box[g] + .5)
    out[..., g] -= turns * box[g]
    out[..., f] -= turns * shift
    out[..., f] -= box[f] * np.floor(out[..., f] / box[f] + .5)
    for axis in range(len(box)):
        if axis not in (f, g):
            out[..., axis] -= box[axis] * np.floor(out[..., axis] / box[axis] + .5)
    best = out.copy()
    norm2 = np.sum(best * best, axis=-1)
    radius = int(np.ceil(np.sqrt(norm2.max(initial=0)) / box[g] + .5))
    for turn in range(-radius, radius + 1):
        if turn == 0:
            continue
        candidate = out.copy()
        candidate[..., g] -= turn * box[g]
        candidate[..., f] -= turn * shift
        candidate[..., f] -= box[f] * np.floor(candidate[..., f] / box[f] + .5)
        trial = np.sum(candidate * candidate, axis=-1)
        better = trial < norm2
        best = np.where(np.asarray(better)[..., None], candidate, best)
        norm2 = np.minimum(norm2, trial)
    return best


def legacy_geometry(parameters, time):
    """Current HEUR input convention, traced to its actual integrator.

    Despite its name, Wi_tau_C is assigned directly to Wi_tau_R in this source.
    With the recorded clock t, shift=Wi_tau_C*L_gradient*t; no extra C factor.
    """
    d = int(parameters['N_dimension'])
    if d != 3:
        raise ValueError('this legacy network schema requires three dimensions')
    box = _box([float(parameters.get('LB'+axis, parameters['box_dimension'])) for axis in 'xyz'], d)
    f = int(parameters.get('shear_axis', 0)); g = int(parameters.get('shear_grad_axis', 1))
    mode = parameters.get('MECHANICAL_PERTURBATION', 'FALSE')
    if mode not in ('FALSE', 'SIMPLE_SHEAR', 'STEP_SHEAR'):
        raise ValueError('unknown legacy perturbation mode')
    if mode != 'FALSE' and (f, g) != (0, 1):
        raise ValueError('current legacy interacting shear supports xy only')
    if mode == 'SIMPLE_SHEAR':
        shift = float(parameters['Wi_tau_C'])*box[g]*float(_finite(time, 'time'))
    elif mode == 'STEP_SHEAR':
        shift = float(parameters['gamma_0'])*box[g]
    else:
        shift = 0.
    return {'box_lengths': box, 'shear_displacement': shift,
            'flow_axis': f, 'gradient_axis': g}


def validate_adjacency(adjacency):
    """Compact rows (indices, weights), self first; loop weight counts two ends."""
    rows = []
    n = len(adjacency)
    if n < 1:
        raise ValueError('at least one particle required')
    for i, (indices, weights) in enumerate(adjacency):
        h = np.asarray(indices); w = np.asarray(weights)
        if (h.ndim != 1 or w.shape != h.shape or not len(h)
                or h.dtype.kind not in 'iu' or w.dtype.kind not in 'iu'
                or h[0] != i or len(set(h.tolist())) != len(h)
                or np.any(h < 0) or np.any(h >= n) or np.any(w < 0)
                or np.any(w[1:] <= 0) or w[0] % 2):
            raise ValueError(f'invalid compact adjacency row {i}')
        rows.append(dict(zip(h.tolist(), w.tolist())))
    for i, row in enumerate(rows):
        for j, weight in row.items():
            if j != i and rows[j].get(i) != weight:
                raise ValueError('adjacency must be reciprocal with equal weights')
    return rows


def topology_statistics(adjacency):
    rows = validate_adjacency(adjacency)
    ends = np.array([sum(row.values()) for row in rows], dtype=np.int64)
    loops = sum(row[i] for i, row in enumerate(rows)) // 2
    bridges = sum(weight for i, row in enumerate(rows) for j, weight in row.items() if j > i)
    distinct = sum(j > i for i, row in enumerate(rows) for j in row)
    unseen = set(range(len(rows))); sizes = []
    while unseen:
        pending = [unseen.pop()]; count = 0
        while pending:
            i = pending.pop(); count += 1
            for j in rows[i]:
                if j != i and j in unseen:
                    unseen.remove(j); pending.append(j)
        sizes.append(count)
    return {'attached_ends': ends, 'loop_chains': loops, 'bridge_chains': bridges,
            'total_chains': loops + bridges, 'distinct_edges': distinct,
            'mean_distinct_degree': 2 * distinct / len(rows),
            'bridge_fraction': bridges / (loops + bridges) if loops + bridges else None,
            'cluster_sizes': np.array(sorted(sizes, reverse=True)),
            'cluster_definition': 'connected particle components, isolated particles included; not spatial percolation'}


@dataclass
class LegacyFrame:
    index: int
    segment: int
    time: float
    positions: np.ndarray
    adjacency: tuple
    chain_ends: np.ndarray | None = None


def iter_legacy_frames(prefix, n_particles, dimension=3, *, start=0, stop=None, stride=1,
                       with_chains=False, allow_time_resets=False):
    """Read synchronized current compact .traj/.hash/.weight[/chain] streams.

    Consume full frames before selection, including skipped frames. Clean EOF,
    partial frames, extra topology and malformed records are distinct. The
    unused auxiliary trajectory columns are never interpreted as velocities.
    """
    n = _integer(n_particles, 'particle count', 1); d = _integer(dimension, 'dimension', 1)
    if d > 3:
        raise ValueError('dimension must be <= 3')
    start = _integer(start, 'start'); stride = _integer(stride, 'stride', 1)
    if stop is not None:
        stop = _integer(stop, 'stop')
        if stop < start:
            raise ValueError('stop must be >= start')
    prefix = str(prefix)
    def ints(line):
        return np.array([int(token) for token in line.split()], dtype=np.int64)
    with ExitStack() as stack:
        traj = stack.enter_context(open(prefix + '.traj'))
        hashes = stack.enter_context(open(prefix + '.hash'))
        weights = stack.enter_context(open(prefix + '.weight'))
        chains = stack.enter_context(open(prefix + '.chain')) if with_chains else None
        index = segment = 0; previous = None
        while stop is None or index < stop:
            line = traj.readline()
            if not line:
                if hashes.read().strip() or weights.read().strip() or (chains and chains.read().strip()):
                    raise ValueError('topology has extra frames beyond trajectory EOF')
                return
            values = _finite([float(token) for token in line.split()], 'trajectory frame')
            if values.shape != (1 + 2 * d * n,):
                raise ValueError('trajectory frame width mismatch')
            time = float(values[0])
            if previous is not None and time <= previous:
                if not allow_time_resets or time >= previous:
                    raise ValueError('nonincreasing clock; split appended runs explicitly')
                segment += 1
            previous = time
            adjacency = []
            for _ in range(n):
                # Sparse rows have variable widths but one row per particle.
                # Read them even for skipped output frames, or trajectory and
                # topology would silently refer to different physical times.
                h = hashes.readline(); w = weights.readline()
                if not h or not w:
                    raise ValueError('partial topology frame')
                adjacency.append((ints(h), ints(w)))
            adjacency = tuple(adjacency)
            rows = validate_adjacency(adjacency)
            ends = None
            if chains:
                raw = chains.readline()
                if not raw:
                    raise ValueError('partial labeled-chain frame')
                ends = ints(raw)
                if len(ends) % 2 or np.any(ends < 0) or np.any(ends >= n):
                    raise ValueError('invalid labeled chain endpoints')
                ends = ends.reshape(2, -1)
                expected = [Counter() for _ in range(n)]
                for a, b in ends.T:
                    expected[a][int(b)] += 1; expected[b][int(a)] += 1
                if any(dict(expected[i]) != {j: w for j, w in row.items() if w}
                       for i, row in enumerate(rows)):
                    raise ValueError('labeled chains disagree with adjacency')
            if index >= start and (index - start) % stride == 0:
                yield LegacyFrame(index, segment, time, values[1:].reshape(n, 2*d)[:, :d].copy(), adjacency, ends)
            index += 1


def bond_observables(positions, adjacency, box_lengths, *, measure='chains', **geometry):
    """Unordered bridge vectors; chain-weighted or distinct-edge measure.

    RR_per_particle matches the historical mathematical purpose. RR_per_bridge
    normalizes by the selected weight. Loops have zero connector and are counted
    separately; they are not quietly mixed into the bridge-only PDF.
    """
    p = _finite(positions, 'positions'); rows = validate_adjacency(adjacency)
    if p.ndim != 2 or len(p) != len(rows) or not 1 <= p.shape[1] <= 3:
        raise ValueError('positions must have shape (particles, dimension)')
    if measure not in ('chains', 'distinct_edges'):
        raise ValueError('measure must be chains or distinct_edges')
    pairs = [(i, j) for i, row in enumerate(rows) for j in row if j > i]
    vectors = np.array([p[j] - p[i] for i, j in pairs]).reshape(-1, p.shape[1])
    vectors = minimum_image(vectors, box_lengths, **geometry)
    weights = np.array([rows[i][j] if measure == 'chains' else 1 for i, j in pairs], dtype=float)
    # Multiple chains sharing one edge affect chain-weighted conformation, but
    # not distinct-neighbor geometry. Both measures are valid for their stated
    # question; do not collapse them into a single unnamed "bond average".
    rr = np.einsum('n,ni,nj->ij', weights, vectors, vectors)
    total = float(weights.sum())
    stats = topology_statistics(adjacency)
    return {'vectors': vectors, 'lengths': np.linalg.norm(vectors, axis=1), 'weights': weights,
            'pairs': np.array(pairs, dtype=int).reshape(-1, 2), 'measure': measure,
            'bridge_weight': total, 'RR_sum': rr, 'RR_per_particle': rr / len(p),
            'RR_per_bridge': rr / total if total else None,
            'loop_chains_excluded': stats['loop_chains']}


def chain_midpoints(frame, box_lengths, **geometry):
    """Instantaneous geometric midpoint proxy; not the polymer's resolved COM."""
    if frame.chain_ends is None:
        raise ValueError('labeled chain endpoints required')
    a, b = frame.chain_ends
    q = minimum_image(frame.positions[b] - frame.positions[a], box_lengths, **geometry)
    return {'midpoints': frame.positions[a] + q / 2, 'connectors': q,
            'definition': 'labeled endpoint midpoint; image/transfer jumps are not COM diffusion'}


def histogram_density(samples, edges, *, weights=None, measure='recorded observations'):
    """Count/weight-normalized PDF with exact cell widths and explicit outside mass.

    A radial PDF divides by dr, a Cartesian PDF by dx*dy*dz. Neither is RDF g.
    Empty selections have an undefined PDF (None), rather than a fabricated zero.
    """
    a = _finite(samples, 'samples')
    if a.ndim == 1:
        a = a[:, None]; edges = [edges]
    if a.ndim != 2:
        raise ValueError('samples must be a vector or a table')
    if len(edges) != a.shape[1]:
        raise ValueError('one bin-edge array per coordinate required')
    bins = []
    for edge in edges:
        edge = _finite(edge, 'bin edges')
        if edge.ndim != 1 or len(edge) < 2 or np.any(np.diff(edge) <= 0):
            raise ValueError('strictly increasing bin edges required')
        bins.append(edge)
    w = np.ones(len(a)) if weights is None else _finite(weights, 'weights')
    if w.shape != (len(a),) or np.any(w < 0):
        raise ValueError('one nonnegative weight per sample required')
    counts, _ = np.histogramdd(a, bins=bins, weights=w)
    volumes = np.ones(counts.shape)
    for axis, edge in enumerate(bins):
        shape = [1] * len(bins); shape[axis] = len(edge) - 1
        volumes *= np.diff(edge).reshape(shape)
    total = float(w.sum())
    return {'counts': counts, 'bin_measure': volumes, 'total_weight': total,
            'outside_weight': max(0., total - float(counts.sum())),
            'pdf': counts / (total * volumes) if total else None,
            'edges': bins, 'measure': measure}


def spherical_bond_density(vectors, radial_edges, mu_edges, phi_edges, *, weights=None):
    """P per physical volume sector: dV=(b^3-a^3)/3*dmu*dphi.

    mu=qz/r and phi=atan2(qy,qx) modulo 2pi. Zero vectors have no direction;
    their mass is explicit and stays in the denominator. No opposite-vector
    symmetrization is invented for labeled/oriented connectors.
    """
    q = _finite(vectors, 'vectors')
    if q.ndim != 2 or q.shape[1] != 3:
        raise ValueError('three-dimensional vectors required')
    w = np.ones(len(q)) if weights is None else _finite(weights, 'weights')
    if w.shape != (len(q),) or np.any(w < 0):
        raise ValueError('one nonnegative weight per vector required')
    r = np.linalg.norm(q, axis=1); nz = r > 0
    mu = q[nz, 2] / r[nz]; phi = np.arctan2(q[nz, 1], q[nz, 0]) % (2*np.pi)
    result = histogram_density(np.column_stack([r[nz], mu, phi]),
                               [radial_edges, mu_edges, phi_edges], weights=w[nz])
    re, me, pe = result['edges']
    if re[0] < 0 or me[0] < -1 or me[-1] > 1 or pe[0] < 0 or pe[-1] > 2*np.pi:
        raise ValueError('invalid spherical coordinate ranges')
    volumes = np.diff(re**3)[:, None, None] / 3 * np.diff(me)[None, :, None] * np.diff(pe)[None, None, :]
    total = float(w.sum())
    result.update(bin_measure=volumes, total_weight=total,
                  pdf=result['counts'] / (total * volumes) if total else None,
                  zero_vector_weight=float(w[~nz].sum()),
                  outside_weight=max(0., total - float(result['counts'].sum())),
                  measure='connector probability per physical volume; zero directions unresolved')
    return result


def _time_grid(times, n):
    t = np.arange(n, dtype=float) if times is None else _finite(times, 'times')
    if t.shape != (n,) or n < 2:
        raise ValueError('at least two matching time samples required')
    dt = np.diff(t)
    if np.any(dt <= 0) or not np.allclose(dt, dt[0], rtol=2e-6, atol=1e-12):
        raise ValueError('uniform increasing time grid required; split appended runs')
    return t, float(dt[0])


def correlation(x, y=None, *, times=None, max_lag=None, center=False, window='all_origins', method='fft'):
    """Linear <x(t)y(t+lag)> per component, with actual origin counts.

    The FFT is zero-padded and conjugated; even/odd lengths use the same lag
    definition. fixed-window means floor(N/2)-1 origins, as declared in legacy.
    """
    x = _finite(x, 'x'); y = x.copy() if y is None else _finite(y, 'y')
    if x.ndim not in (1, 2) or y.shape != x.shape:
        raise ValueError('equal vectors or (time,components) tables required')
    n = len(x); _, dt = _time_grid(times, n)
    if window not in ('all_origins', 'legacy_fixed_window') or method not in ('fft', 'direct'):
        raise ValueError('unknown correlation window or method')
    limit = n if window == 'all_origins' else n // 2
    if window == 'legacy_fixed_window' and limit < 2:
        raise ValueError('legacy window requires at least four samples')
    lag = max(0, n//2 - 1) if max_lag is None else _integer(max_lag, 'maximum lag')
    if lag >= limit:
        raise ValueError('maximum lag exceeds origin window')
    if center:
        x = x - x.mean(axis=0); y = y - y.mean(axis=0)
    counts = n - np.arange(lag + 1) if window == 'all_origins' else np.full(lag + 1, limit - 1)
    if method == 'fft':
        # Zero padding prevents circular wrap-around of late data into early
        # lags. The per-lag origin denominator is part of the estimator, rather
        # than a plotting normalization that can be changed after the FFT.
        nfft = 1 << (2*n - 2).bit_length()
        first = x.copy()
        if window == 'legacy_fixed_window':
            first[limit - 1:] = 0
        products = np.fft.irfft(np.conj(np.fft.rfft(first, n=nfft, axis=0)) *
                               np.fft.rfft(y, n=nfft, axis=0), n=nfft, axis=0)[:lag + 1]
        values = products / counts.reshape((-1,) + (1,) * (x.ndim - 1))
    else:
        values = np.array([np.sum(x[:m] * y[k:k+m], axis=0) / m for k, m in enumerate(counts)])
    return {'lag_times': np.arange(lag + 1) * dt, 'values': values, 'origin_counts': counts,
            'centered': bool(center), 'window': window, 'method': method,
            'definition': '<x(t) y(t+lag)>; component products, no C(0) normalization'}


def stress_correlation(stress, other=None, *, projection='isotropic', center=True, **options):
    """Symmetric 3D stress columns xx,yy,zz,xy,xz,yz.

    Isotropic projection = sum(offdiag correlations)/5 +
    sum(normal-difference correlations)/30 = <dev(A):dev(B)>/10.
    It is a projection in a driven state, not a nonequilibrium Green–Kubo law.
    """
    a = _finite(stress, 'stress'); b = a if other is None else _finite(other, 'other stress')
    if a.ndim != 2 or a.shape[1] != 6 or b.shape != a.shape:
        raise ValueError('matching six-column stress tables required')
    if projection == 'xy':
        result = correlation(a[:, 3], b[:, 3], center=center, **options)
    elif projection == 'isotropic':
        def modes(s):
            return np.column_stack([s[:, 3:]/np.sqrt(5),
                                    (s[:, 0]-s[:, 1])/np.sqrt(30),
                                    (s[:, 0]-s[:, 2])/np.sqrt(30),
                                    (s[:, 1]-s[:, 2])/np.sqrt(30)])
        result = correlation(modes(a), modes(b), center=center, **options)
        result['values'] = result['values'].sum(axis=1)
    else:
        raise ValueError('projection must be isotropic or xy')
    result.update(projection=projection, units='stress squared', observable='configurational stress correlation')
    return result


def stress_components(energy, *, schema):
    """Explicit schema; recorded random-force RF is not configurational stress."""
    a = _finite(energy, 'energy table')
    widths = {'heur31': 31, 'repulsive24': 24, 'dumbbell24': 24, 'brownian18': 18}
    if schema not in widths or a.ndim != 2 or a.shape[1] != widths[schema]:
        raise ValueError('table width disagrees with explicit energy schema')
    result = {'time': a[:, 0].copy(), 'random_force_diagnostic': a[:, 12:18].copy()}
    if schema in ('heur31', 'repulsive24'):
        result['repulsive'] = a[:, 18:24].copy()
    if schema in ('heur31', 'dumbbell24'):
        result['connector'] = a[:, 24:30].copy() if schema == 'heur31' else a[:, 18:24].copy()
    result['configurational'] = sum((result[k] for k in ('repulsive', 'connector') if k in result),
                                    np.zeros((len(a), 6)))
    return result


def green_kubo(corr, *, volume, thermal_energy, equilibrium):
    """Finite configurational contribution G=V*C/kBT and its running integral.

    Caller declares equilibrium and kBT in the same energy units as stress*V.
    This does not add solvent/instantaneous terms or select a convergence window.
    """
    if equilibrium is not True or not corr.get('centered'):
        raise ValueError('equilibrium and connected (centered) correlation required')
    v = _positive(volume, 'volume'); theta = _positive(thermal_energy, 'thermal energy')
    c = _finite(corr['values'], 'correlation'); t = _finite(corr['lag_times'], 'lag times')
    if c.ndim != 1 or t.shape != c.shape or not len(t) or t[0] != 0 or np.any(np.diff(t) <= 0):
        raise ValueError('scalar correlation on increasing nonnegative lag grid required')
    g = v / theta * c
    eta = np.r_[0., np.cumsum((g[1:] + g[:-1]) * np.diff(t) / 2)]
    return {'lag_times': t.copy(), 'modulus_contribution': g, 'running_viscosity_contribution': eta,
            'definition': 'finite configurational equilibrium contribution; convergence not inferred'}


def unwrap_positions(positions, box_lengths, *, images=None, assume_small_steps=False,
                     shear_displacements=None, flow_axis=0, gradient_axis=1):
    """Unwrap using exact image counters, or declared small-step orthogonal sampling.

    Coarse snapshots alone do not determine winding. Evolving shear requires
    image counters; guessing boundary histories from two snapshots is refused.
    """
    p = _finite(positions, 'positions')
    if p.ndim != 3 or not len(p) or not p.shape[1] or not 1 <= p.shape[2] <= 3:
        raise ValueError('(frames,particles,dimension) positions required')
    box = _box(box_lengths, p.shape[2])
    if images is not None:
        # Exact image history is additional information carried by the recorder.
        # A nearest-image difference between coarse frames cannot reconstruct
        # an unknown number of crossings, especially under evolving shear.
        im = np.asarray(images)
        if im.shape != p.shape or im.dtype.kind not in 'iu':
            raise ValueError('matching integer image counters required')
        result = p + im * box
        if shear_displacements is not None:
            shifts = _finite(shear_displacements, 'shear displacements')
            f = _integer(flow_axis, 'flow axis'); g = _integer(gradient_axis, 'gradient axis')
            if shifts.shape != (len(p),) or f == g or max(f, g) >= p.shape[2]:
                raise ValueError('matching per-frame shear and valid axes required')
            result[:, :, f] += im[:, :, g] * shifts[:, None]
        return result
    if shear_displacements is not None:
        raise ValueError('sheared trajectory unwrapping requires image counters')
    if assume_small_steps is not True:
        raise ValueError('image counters or explicit assume_small_steps=True required')
    delta = np.diff(p, axis=0)
    reduced = minimum_image(delta, box)
    if np.any(np.isclose(np.abs(reduced), box / 2, rtol=0, atol=1e-12)):
        raise ValueError('half-box jump has ambiguous winding')
    return np.concatenate([p[:1], p[:1] + np.cumsum(reduced, axis=0)])


def mean_square_displacement(positions, *, times=None, max_lag=None, origin='all',
                             affine_displacement=None):
    """MSD/tensor of explicitly unwrapped positions; no implicit PBC inversion.

    Optional affine_displacement is the supplied integrated flow history at
    every particle/time. Single-origin output uses lag times from zero.
    """
    p = _finite(positions, 'unwrapped positions')
    if p.ndim != 3 or not p.shape[1] or not 1 <= p.shape[2] <= 3:
        raise ValueError('(frames,particles,dimension) unwrapped positions required')
    _, dt = _time_grid(times, len(p))
    if affine_displacement is not None:
        flow = _finite(affine_displacement, 'integrated affine displacement')
        if flow.shape != p.shape:
            raise ValueError('matching integrated flow history required')
        p = p - flow
    if origin not in ('all', 'first'):
        raise ValueError('origin must be all or first')
    lag = len(p)-1 if max_lag is None else _integer(max_lag, 'maximum lag')
    if lag >= len(p):
        raise ValueError('maximum lag exceeds trajectory')
    tensors = []; counts = []
    for k in range(lag + 1):
        m = len(p)-k if origin == 'all' else 1
        delta = p[k:k+m] - p[:m]
        tensors.append(np.einsum('tni,tnj->ij', delta, delta) / (m * p.shape[1])); counts.append(m)
    tensor = np.array(tensors)
    return {'lag_times': np.arange(lag + 1)*dt, 'tensor': tensor,
            'msd': np.trace(tensor, axis1=1, axis2=2), 'origin_counts': np.array(counts),
            'origin': origin, 'coordinates': 'explicitly unwrapped; supplied affine history removed' if affine_displacement is not None else 'explicitly unwrapped'}


def structure_factor(positions, wavevectors):
    """Direct S(q)=mean_frames |sum_j exp(i q.r_j)|^2/N, self term included.

    q=0 is the finite-system value N. Periodic interpretation requires reciprocal
    lattice wavevectors; arbitrary q still defines the stated finite point sum.
    """
    p = _finite(positions, 'positions'); q = _finite(wavevectors, 'wavevectors')
    if p.ndim == 2:
        p = p[None, ...]
    if (p.ndim != 3 or not len(p) or not p.shape[1] or q.ndim != 2
            or q.shape[1] != p.shape[2] or not 1 <= p.shape[2] <= 3):
        raise ValueError('matching position and wavevector dimensions required')
    values = []
    for vector in q:
        rho = np.exp(1j * (p @ vector)).sum(axis=1)
        values.append(np.mean(np.abs(rho)**2) / p.shape[1])
    return np.array(values)


def rdf_to_structure_factor(r, g, wave_numbers, *, number_density, dimension=3):
    """Continuum isotropic transform of supplied bulk g; finite cutoff is explicit.

    3D: 1+4pi*rho*integral r^2*(g-1)*sinc(qr/pi) dr.
    2D: 1+2pi*rho*integral r*(g-1)*J0(qr) dr (optional SciPy).
    q=0 is the integral limit, not the fixed-N direct forward-scattering peak.
    """
    r = _finite(r, 'r'); g = _finite(g, 'g'); q = _finite(wave_numbers, 'wave numbers')
    rho = _positive(number_density, 'number density')
    if (r.ndim != 1 or len(r) < 2 or g.shape != r.shape or r[0] < 0
            or np.any(np.diff(r) <= 0) or q.ndim != 1 or np.any(q < 0) or dimension not in (2, 3)):
        raise ValueError('valid radial samples, wave numbers and dimension required')
    if dimension == 3:
        kernel = 4*np.pi*r[None, :]**2 * np.sinc(q[:, None]*r/np.pi)
    else:
        try:
            from scipy.special import j0
        except ImportError as exc:
            raise ImportError('2D isotropic RDF transform requires scipy.special.j0') from exc
        kernel = 2*np.pi*r[None, :] * j0(q[:, None]*r)
    integrand = kernel*(g-1)
    integral = np.sum((integrand[:, 1:] + integrand[:, :-1])*np.diff(r)/2, axis=1)
    return {'values': 1 + rho*integral,
            'wave_numbers': q.copy(), 'radial_range': [float(r[0]), float(r[-1])],
            'definition': 'truncated continuum transform of supplied bulk g; no tail extrapolation'}


def conservative_virial(positions, adjacency, box_lengths, *, connector='modified_gaussian',
                        alpha=1., repulsion_coefficient=1., repulsion_coupling='pre_average',
                        repulsion_base=.25, fene_maximum=None, effective_distance=1., **geometry):
    """Reconstruct current HEUR conservative energy and stored RF/V convention.

    Uc/kBT=d*alpha^2*r^2/2 (or FENE); Ur/(C*kBT)=(2-3r+r^3)/3
    for r<1. All distinct particle pairs enter repulsion, not only linked pairs.
    Stored RF is +q*(dU/dr)*q/r/V; pressure virial has the opposite sign.
    Energy/stress units here are C*kBT and C*kBT/length^d respectively.
    """
    p = _finite(positions, 'positions')
    if p.ndim != 2 or not len(p) or not 1 <= p.shape[1] <= 3:
        raise ValueError('particle position table required')
    box = _box(box_lengths, p.shape[1]); volume = float(np.prod(box))
    c = _positive(repulsion_coefficient, 'repulsion coefficient'); a = _positive(alpha, 'alpha')
    if effective_distance != 1.:
        raise ValueError('legacy repulsive polynomial is consistent at effective_distance=1 only')
    if connector not in ('gaussian', 'modified_gaussian', 'fene'):
        raise ValueError('unknown connector law')
    if repulsion_coupling not in ('pre_average', 'geometric_mean'):
        raise ValueError('unknown repulsive coupling')
    base = _nonnegative(repulsion_base, 'repulsion base')
    bonds = bond_observables(p, adjacency, box, measure='chains', **geometry)
    r = bonds['lengths']; q = bonds['vectors']; w = bonds['weights']
    stiffness = p.shape[1] * (1. if connector == 'gaussian' else a*a)
    if connector == 'fene':
        maximum = _positive(fene_maximum, 'FENE maximum')
        if np.any(r >= maximum):
            raise ValueError('FENE connector reaches/exceeds maximum extension')
        factor = stiffness / (1 - (r/maximum)**2)
        u = -.5*stiffness*maximum**2*np.log1p(-(r/maximum)**2)
    else:
        factor = np.full_like(r, stiffness); u = .5*stiffness*r*r
    con = np.einsum('n,ni,nj->ij', w*factor/c, q, q) / volume
    ucon = float(w @ u / c)
    ii, jj = np.triu_indices(len(p), 1)
    delta = minimum_image(p[jj]-p[ii], box, **geometry); lengths = np.linalg.norm(delta, axis=1)
    take = lengths < 1.; rr = lengths[take]; dr = delta[take]
    if np.any(rr == 0):
        raise ValueError('coincident distinct particles have undefined radial repulsive force')
    coupling = np.ones(len(rr))
    if repulsion_coupling == 'geometric_mean':
        ends = topology_statistics(adjacency)['attached_ends']
        coupling = base/c * ends[ii[take]] * ends[jj[take]]
    rep = np.einsum('n,ni,nj->ij', coupling*(rr*rr-1)/rr, dr, dr) / volume
    urep = float(np.sum(coupling*(2-3*rr+rr**3)/3))
    return {'connector': con, 'repulsive': rep, 'configurational': con+rep,
            'connector_energy': ucon, 'repulsive_energy': urep,
            'stress_convention': 'stored RF/V = +q outer force_on_i / V; negative pressure virial',
            'thermal_energy_in_recorded_units': 1/c}


def bell_factor(lengths, *, dimension=3, alpha=1., capture_length, weights=None):
    """Bridge-conditioned mean exp(tension*l_cap), for a Gaussian connector.

    This is not an automatically inferred MC acceptance probability or rate.
    """
    r = _finite(lengths, 'lengths'); a = _positive(alpha, 'alpha')
    d = _integer(dimension, 'dimension', 1); cap = _nonnegative(capture_length, 'capture length')
    if r.ndim != 1 or np.any(r < 0):
        raise ValueError('nonnegative distance vector required')
    w = np.ones(len(r)) if weights is None else _finite(weights, 'weights')
    if w.shape != r.shape or np.any(w < 0) or not w.sum():
        raise ValueError('nonempty positive total weight required')
    exponent = d*a*a*r*cap
    if np.any(exponent > np.log(np.finfo(float).max)):
        raise ValueError('Bell factor overflows; retain a log-domain statistic')
    return float(np.average(np.exp(exponent), weights=w))


def connector_correlation(connectors, *, center=False, **options):
    """Labeled chain vector dot correlation, averaged over the same chain IDs.

    Chain transfer changes q for that chain and is retained. This is neither
    continuous-bond survival nor correlation of unlabeled edge lists.
    """
    q = _finite(connectors, 'labeled connectors')
    if q.ndim != 3 or not q.shape[1] or not 1 <= q.shape[2] <= 3:
        raise ValueError('(time,chains,dimension) labeled connectors required')
    result = correlation(q.reshape(len(q), -1), center=center, **options)
    per_chain = result['values'].reshape(-1, q.shape[1], q.shape[2]).sum(axis=2)
    result.update(values=per_chain.mean(axis=1), per_chain=per_chain,
                  definition='mean over labeled chains of q(t) dot q(t+lag)')
    return result


def empirical_cdf(samples, *, weights=None):
    """Right-continuous weighted empirical CDF; duplicate observations combined.

    A spline derivative is not used as a PDF because it need not be positive.
    """
    x = _finite(samples, 'samples')
    w = np.ones(len(x)) if weights is None else _finite(weights, 'weights')
    if x.ndim != 1 or w.shape != x.shape or np.any(w < 0) or not w.sum():
        raise ValueError('sample vector and positive total nonnegative weight required')
    unique, inverse = np.unique(x, return_inverse=True)
    mass = np.bincount(inverse, weights=w, minlength=len(unique))
    return {'x': unique, 'cdf': np.cumsum(mass)/w.sum(), 'mass': mass/w.sum(),
            'definition': 'P(recorded value <= x), selected measure'}


def diffusion_fit(msd, lag_times, *, dimension, fit_range):
    """Explicit-window Einstein slope/(2*d), with intercept retained.

    Assumes unwrapped equilibrium diffusive MSD. No automatic regime selection,
    uncertainty estimate, or interpretation of driven/ballistic MSD is supplied.
    """
    values = _finite(msd, 'MSD'); t = _finite(lag_times, 'lag times')
    d = _integer(dimension, 'dimension', 1); limits = _finite(fit_range, 'fit range')
    if (t.ndim != 1 or values.shape != t.shape or np.any(np.diff(t) <= 0)
            or np.any(values < 0) or limits.shape != (2,) or limits[0] >= limits[1]):
        raise ValueError('increasing times, matching nonnegative MSD and fit range required')
    select = (t >= limits[0]) & (t <= limits[1])
    if np.count_nonzero(select) < 2:
        raise ValueError('fit range needs at least two lag samples')
    design = np.column_stack([t[select], np.ones(np.count_nonzero(select))])
    slope, intercept = np.linalg.lstsq(design, values[select], rcond=None)[0]
    if slope < 0:
        raise ValueError('negative fitted slope cannot define an equilibrium diffusion coefficient')
    return {'diffusion': float(slope/(2*d)), 'slope': float(slope), 'intercept': float(intercept),
            'fit_range': limits.copy(), 'samples': int(select.sum()),
            'interpretation': 'Einstein estimate under declared diffusive equilibrium conditions'}


def distance_stream(path, *, legacy_scientific_printer=False):
    """Ragged selected-distance stream; optional explicit fixed-printer joins.

    This yields recorded samples only. It never invents zeros/frame boundaries.
    """
    pattern = re.compile(r'[+-]?\d\.\d{6}e[+-]\d{2}')
    with Path(path).open() as stream:
        for line_number, line in enumerate(stream, 1):
            if legacy_scientific_printer:
                end = 0; tokens = []
                for match in pattern.finditer(line):
                    if line[end:match.start()].strip():
                        raise ValueError(f'unrecognized distance at line {line_number}')
                    tokens.append(match.group()); end = match.end()
                if line[end:].strip():
                    raise ValueError(f'unrecognized distance suffix at line {line_number}')
            else:
                tokens = line.split()
            for token in tokens:
                value = float(token)
                if not np.isfinite(value) or value < 0:
                    raise ValueError(f'invalid distance at line {line_number}')
                yield value
