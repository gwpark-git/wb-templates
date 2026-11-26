"""Sample synthetic reduced properties outside the continuum solve."""

import argparse
import csv
import math
import sys


PHI_NODES = (0., 0.1, 0.2, 0.3, 0.4, 0.494)


def analytical_pair(phi):
    """Return the synthetic formula used as the inference oracle."""
    return 1. + 1.5 * phi, 1. + 2.5 * phi + 4. * phi * phi


def torch_pairs(nodes):
    """Evaluate a batch of states on CPU with explicit float64 inputs."""
    import torch

    class DemoTransport(torch.nn.Module):
        def forward(self, phi):
            diffusivity = 1. + 1.5 * phi
            viscosity = 1. + 2.5 * phi + 4. * phi * phi
            return torch.cat((diffusivity, viscosity), dim=1)

    model = DemoTransport().eval()
    with torch.inference_mode():
        features = torch.tensor(nodes, dtype=torch.float64, device="cpu").reshape(-1, 1)
        outputs = model(features).tolist()
    return outputs


def table(backend="python"):
    """Return aligned rows of volume fraction, D_c/D_0 and eta/eta_s."""
    if backend == "python":
        values = [analytical_pair(phi) for phi in PHI_NODES]
    elif backend == "torch":
        values = torch_pairs(PHI_NODES)
    else:
        raise ValueError(f"Unsupported example backend: {backend}")
    return [(phi, diffusivity, viscosity)
            for phi, (diffusivity, viscosity) in zip(PHI_NODES, values, strict=True)]


def linear_pair(phi, nodes, diffusivity, viscosity):
    """Interpolate validated ordered columns, with no automatic extrapolation."""
    if not math.isfinite(phi) or not nodes[0] <= phi <= nodes[-1]:
        raise ValueError("Volume fraction is outside the qualified property table")
    upper = 1
    while nodes[upper] < phi:
        upper += 1
    lower = upper - 1
    fraction = (phi - nodes[lower]) / (nodes[upper] - nodes[lower])
    return (diffusivity[lower] + fraction * (diffusivity[upper] - diffusivity[lower]),
            viscosity[lower] + fraction * (viscosity[upper] - viscosity[lower]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--backend", choices=("python", "torch"), default="python")
    args = parser.parse_args()
    try:
        rows = table(args.backend)
    except ModuleNotFoundError as error:
        parser.error(f"The torch backend requires a compatible PyTorch installation: {error}")
    writer = csv.writer(sys.stdout, lineterminator="\n")
    writer.writerow(("phi", "D_c_over_D0", "eta_over_eta_s"))
    writer.writerows(tuple(f"{value:.16e}" for value in row) for row in rows)


if __name__ == "__main__":
    main()
