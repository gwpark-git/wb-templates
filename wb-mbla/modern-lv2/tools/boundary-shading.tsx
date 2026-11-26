import { useMemo } from 'react';

export type BoundaryProfiles = { wall: number[][]; outlet: number[][] };

// Use the site's log colour scale only on the two recorded boundaries.
export function BoundaryShading({ profiles, x, y }: {
  profiles: BoundaryProfiles;
  x: (z: number) => number;
  y: (distance: number) => number;
}) {
  const cells = useMemo(() => {
    const logPhi = (value: number) => Math.log(Math.max(value, 1e-15));
    const logs = [...profiles.wall, ...profiles.outlet].map(row => logPhi(row[1]));
    const lower = Math.min(...logs);
    const span = Math.max(...logs) - lower;
    const color = (phi: number) => {
      const f = span ? (logPhi(phi) - lower) / span : 0;
      return `rgb(${Math.round(244 - 100 * f)},${Math.round(242 - 164 * f)},${Math.round(233 - 192 * f)})`;
    };
    // The rectangles are display cells centered on the exact stored coordinates.
    const cellsFor = (points: number[][]) => points.map(([at, phi], index) => ({
      from: index ? (points[index - 1][0] + at) / 2 : at,
      to: index + 1 < points.length ? (points[index + 1][0] + at) / 2 : at,
      color: color(phi),
    }));
    return { wall: cellsFor(profiles.wall), outlet: cellsFor(profiles.outlet) };
  }, [profiles]);
  return <g data-boundary-shading="true" shapeRendering="crispEdges">
    {cells.wall.map((cell, index) => <rect key={`wall-${index}`} x={x(cell.from)}
      y={y(0) - 11} width={x(cell.to) - x(cell.from)} height="11" fill={cell.color} />)}
    {cells.outlet.map((cell, index) => <rect key={`outlet-${index}`} x={x(1) - 11}
      y={y(cell.to)} width="11" height={y(cell.from) - y(cell.to)} fill={cell.color} />)}
  </g>;
}
