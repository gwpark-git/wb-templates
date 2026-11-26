#!/usr/bin/env python3
"""Read the inherited ten-column SI table and expose its hydraulic/flux checks."""
import argparse
import csv
import json
import math
from pathlib import Path

COLUMNS = ['z_m', 'phi_wall', 'pressure_Pa', 'permeate_m_s', 'axial_factor_m_s',
           'osmotic_pressure_Pa', 'particle_flux_m3_s', 'cake_ratio_times_radius_m',
           'transmembrane_pressure_Pa', 'cake_thickness_over_radius']


def inspect(path, radius, permeability):
    rows = [[float(x) for x in line.split()] for line in path.read_text().splitlines() if line.strip()]
    if radius <= 0 or permeability <= 0:
        raise ValueError('radius and clean permeability must be positive')
    if not rows or any(len(r) != 10 or not all(map(math.isfinite, r)) for r in rows):
        raise ValueError('expected a nonempty finite ten-column table')
    if any(rows[i+1][0] <= rows[i][0] for i in range(len(rows)-1)):
        raise ValueError('axial coordinates must increase strictly')
    if rows[0][6] == 0:
        raise ValueError('zero inlet particle flux; relative diagnostic undefined')
    cake = [r[7] / radius for r in rows]
    if min(cake) < 0:
        raise ValueError('negative cake resistance')
    predicted = [permeability*(r[8]-r[5])/(1+resistance) for r, resistance in zip(rows, cake)]
    velocity_scale = max(abs(r[3]) for r in rows)
    pperm = [r[2]-r[8] for r in rows]
    return rows, {
        'rows': len(rows), 'columns': COLUMNS,
        'phi_wall_outlet': rows[-1][1], 'cake_resistance_ratio_outlet': cake[-1],
        'max_recorded_particle_flux_relative_drift': max(abs(r[6]/rows[0][6]-1) for r in rows),
        'darcy_starling_max_residual_normalized_by_max_wall_velocity':
            max(abs(v-r[3]) for v,r in zip(predicted,rows))/velocity_scale if velocity_scale else None,
        'inferred_permeate_pressure_range_Pa': [min(pperm),max(pperm)],
        'limits': 'Printed-value diagnostics only; iterative stopping and conservation/grid acceptance are separate.'}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('table', type=Path)
    p.add_argument('--radius', type=float, required=True, help='channel radius in metres')
    p.add_argument('--permeability', type=float, required=True, help='clean Lp in m/(Pa s)')
    p.add_argument('--csv', type=Path, help='explicit optional CSV destination; overwritten')
    a = p.parse_args()
    rows, result = inspect(a.table,a.radius,a.permeability)
    if a.csv:
        with a.csv.open('w', newline='') as out:
            writer=csv.writer(out);writer.writerow(COLUMNS);writer.writerows(rows)
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
