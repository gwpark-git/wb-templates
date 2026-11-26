"""Run the bounded radial thermal pair benchmark from free-space mobility grids."""

import argparse
import csv
from pathlib import Path

import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.special import exprel

MODELS = ("pure_BD_K1", "RPY_BD", "thermal_SD_pair")
SEEDS = (51367, 93481)
LAGS = np.linspace(0, 0.4, 51)


def write_rows(path, labels, rows):
    with path.open("x", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(labels)
        for row in rows:
            if not all(isinstance(value, str) or np.isfinite(value) for value in row):
                raise ValueError(f"nonfinite output for {path}")
            writer.writerow(value if isinstance(value, str) else f"{value:.16e}" for value in row)


def coefficients(table, model):
    if model == MODELS[0]:
        ones = np.ones(table.shape[0])
        return 2*ones, 2*ones, .5*ones, .5*ones
    start = 3 if model == MODELS[1] else 7
    return table[:, start:start+4].T


def radial_generator(radius, weight, diffusion_parallel):
    h = radius[2] - radius[0]
    center_weights = weight[1::2]
    conductance = weight[2:-1:2] * diffusion_parallel[2:-1:2] / h**2
    plus = conductance / center_weights[:-1]
    minus = conductance / center_weights[1:]
    diagonal = -np.r_[plus, 0.] - np.r_[0., minus]
    offdiag = conductance / np.sqrt(center_weights[:-1]*center_weights[1:])
    stationary = center_weights / center_weights.sum()
    np.testing.assert_allclose(
        diagonal*np.sqrt(stationary)
        + np.r_[offdiag*np.sqrt(stationary[1:]), 0.]
        + np.r_[0., offdiag*np.sqrt(stationary[:-1])],
        0., atol=1e-10,
    )
    return diagonal, offdiag, stationary


def thermal_samples(diagonal, offdiag, stationary, radius, values):
    """Sample endpoints of the reversible continuous-time dynamics at 2 τσ."""
    eigenvalues, eigenvectors = eigh_tridiagonal(diagonal, offdiag)
    assert abs(eigenvalues[-1]) < 1e-9
    equilibrium_values = stationary @ values
    rows = []
    for seed, start in zip(SEEDS, (1.15, 2.4)):
        source = int(np.argmin(abs(radius - start)))
        decay = np.exp(2. * eigenvalues)
        probability = eigenvectors @ (decay * eigenvectors[source])
        probability *= np.sqrt(stationary / stationary[source])
        assert probability.min() >= -1e-12 and abs(probability.sum() - 1.) < 1e-9
        probability = np.maximum(probability, 0.)
        probability /= probability.sum()
        total_variation = .5 * np.sum(abs(probability - stationary))
        assert total_variation < .003, (seed, total_variation)
        ids = np.random.default_rng(seed).choice(radius.size, size=8192, p=probability)
        measurements = values[ids]
        estimate = measurements.mean(axis=0)
        se = measurements.std(axis=0, ddof=1) / np.sqrt(len(ids))
        assert np.all(abs(estimate-equilibrium_values) <= 6*se + .01*abs(equilibrium_values) + 1e-12)
        rows.append([str(seed), radius[source], total_variation, *estimate, *se])
    return rows



def sampled_green_kubo(diagonal, offdiag, stationary, virial, angular_decay, prefactor):
    """Check the spectral result with 4096 equilibrium radial paths per seed."""
    eigenvalues, basis = eigh_tridiagonal(diagonal, offdiag)
    dt = LAGS[1] - LAGS[0]
    reversible = (basis*np.exp(dt*eigenvalues)) @ basis.T
    transition = reversible * np.sqrt(stationary[None, :] / stationary[:, None])
    assert transition.min() > -1e-12
    np.testing.assert_allclose(transition.sum(axis=1), 1., atol=1e-10)
    cdf = np.cumsum(np.maximum(0., transition), axis=1)
    cdf /= cdf[:, -1, None]
    rows = []
    for seed in SEEDS:
        rng = np.random.default_rng(seed+1)
        states = rng.choice(stationary.size, size=4096, p=stationary)
        start = virial[states]
        attenuation = np.ones(states.size)
        previous = start*start
        integral = np.zeros(states.size)
        for _ in LAGS[1:]:
            updated = np.sum(cdf[states] < rng.random(states.size)[:, None], axis=1)
            attenuation *= np.exp(-.5*dt*(angular_decay[states]+angular_decay[updated]))
            correlation = start*virial[updated]*attenuation
            integral += .5*dt*(previous+correlation)
            states, previous = updated, correlation
        values = prefactor*integral
        rows.append((str(seed+1), values.mean(), values.std(ddof=1)/np.sqrt(values.size)))
    return rows


def analyze(table, model, sample):
    radius, energy, virial = table[:, :3].T
    weight = radius**2 * np.exp(-energy)
    radial, tangential, cm_radial, cm_tangential = coefficients(table, model)
    diagonal, offdiag, stationary = radial_generator(radius, weight, radial)
    r = radius[1::2]
    volume = (4*np.pi/3) * (radius[-1]**3 - radius[0]**3)
    pair_h = 2*(cm_radial[1::2] + 2*cm_tangential[1::2])/3
    internal_pressure = virial[1::2]/(3*volume)
    # The outer thermodynamic load differentiates the relative-coordinate
    # partition function at fixed inner wall. It includes no COM confinement.
    outer_pressure = np.exp(-energy[-1]) / (4*np.pi*weight[1::2].sum()*(radius[2]-radius[0]))

    # Given a radial path, the l=2 angular correlation contributes exp(-∫6 D_t/r²).
    # The corresponding killed radial generator integrates that factor directly.
    angular_decay = 6*tangential[1::2]/r**2
    decay, basis = eigh_tridiagonal(diagonal - angular_decay, offdiag)
    assert decay[-1] < 0.
    spectral_weight = (basis.T @ (np.sqrt(stationary)*virial[1::2]))**2
    prefactor = 3*np.pi / (15*volume)  # Green–Kubo tensor contraction and solvent scale.
    running = prefactor * (LAGS * (spectral_weight @ exprel(decay[:, None]*LAGS)))
    long_limit = prefactor * np.sum(spectral_weight / -decay)
    high_frequency = (stationary @ table[1::2, -1]) / (2*(1/(3*np.pi))*volume)
    metrics = [stationary @ pair_h, stationary @ r, stationary @ internal_pressure,
               outer_pressure, running[-1], long_limit]
    samples = thermal_samples(
        diagonal, offdiag, stationary, r, np.column_stack((r, pair_h, internal_pressure)),
    ) if sample else []
    stress = sampled_green_kubo(
        diagonal, offdiag, stationary, virial[1::2], angular_decay, prefactor,
    ) if sample else []
    for _, estimate, se in stress:
        assert abs(estimate-running[-1]) <= 6*se + .03*running[-1]
    return metrics, running, high_frequency, samples, stress


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grids", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir()
    rows, samples, curves, stress_checks, spectra = [], [], [], [], {}
    for cells in (72, 144, 288):
        table = np.loadtxt(args.grids / f"mobility-{cells}.csv", delimiter=",", skiprows=1)
        assert table.shape == (2*cells+1, 12) and np.isfinite(table).all()
        for model in MODELS:
            metrics, running, high_frequency, trajectories, stress = analyze(table, model, cells == 144)
            rows.append([str(cells), model, *metrics,
                         high_frequency if model == MODELS[2] else ""])
            spectra[cells, model] = np.asarray(metrics)
            if cells == 144:
                curves.extend([model, lag, value] for lag, value in zip(LAGS, running))
                samples.extend([model, *row] for row in trajectories)
                stress_checks.extend([model, *row] for row in stress)
                print(model, " ".join(f"{x:.6e}" for x in metrics))
        expected = spectra[cells, MODELS[0]]
        np.testing.assert_allclose(expected[0], 1., atol=2e-15)
        # Across dynamics the Boltzmann density and pressure are the same.
        for model in MODELS[1:]:
            np.testing.assert_array_equal(spectra[cells, model][1:4], expected[1:4])
    for model in MODELS:
        np.testing.assert_allclose(spectra[144, model], spectra[288, model], rtol=.008)
    write_rows(args.output/"pair-summary.csv",
               ("cells", "model", "pair_cm_mobility_K1_normalized", "mean_r_over_sigma",
                "soft_virial_pressure_star", "outer_confinement_pressure_star",
                "eta_conf_over_eta_s_cutoff_0p4", "eta_conf_over_eta_s_pair_limit",
                "sd_eta_high_frequency_increment_over_eta_s"), rows)
    write_rows(args.output/"thermal-endpoints.csv",
               ("model", "seed", "initial_r_over_sigma", "tv_to_equilibrium_at_2tau",
                "mean_r_over_sigma", "pair_cm_mobility_K1_normalized", "soft_virial_pressure_star",
                "se_r", "se_mobility", "se_virial_pressure"), samples)
    write_rows(args.output/"sampled-green-kubo.csv",
               ("model", "seed", "eta_conf_over_eta_s_cutoff_0p4", "standard_error"), stress_checks)
    write_rows(args.output/"stress-relaxation.csv",
               ("model", "cutoff_over_tau_sigma", "eta_conf_over_eta_s"), curves)


if __name__ == "__main__":
    main()
