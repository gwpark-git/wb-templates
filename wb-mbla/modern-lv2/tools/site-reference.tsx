import { useState } from 'react';
import { FiltrationSchematic } from './schematic';
import { LinePlot } from './plots';
import { cpReferenceScale } from './cp-scale';
import { RadialCoordinates } from './radial-coordinates';
import { axialProfileAtY, profileAtZ } from './sampling';
import { screenToWall, wallTicks, wallToScreen, type RadialScale } from './coordinates';

type Case = NonNullable<Parameters<typeof FiltrationSchematic>[0]['result']> &
  Parameters<typeof profileAtZ>[0] & Parameters<typeof cpReferenceScale>[0] & {
  source: { id: string; sha256: string; commit: string; nominalTmp: number; radiusNm: number; bulkFraction: number };
};

export function SiteReference() {
  const [{ data, cp }] = useState(() => {
    const data = JSON.parse(document.getElementById('site-reference-data')!.textContent!) as Case;
    return { data, cp: cpReferenceScale(data) };
  });
  const [z, setZ] = useState(0.5);
  const [y, setY] = useState(1 - Math.sqrt(2 / 3));
  const [scale, setScale] = useState<RadialScale>('wall-log');
  const [velocity, setVelocity] = useState<'u_m_s' | 'v_m_s'>('v_m_s');
  const velocityLabel = velocity === 'v_m_s' ? 'Transverse velocity v (µm/s)' : 'Axial velocity u (mm/s)';
  const multiplier = velocity === 'v_m_s' ? 1e6 : 1e3;
  const radialScale = {
    forward: (value: number) => wallToScreen(value, 1, scale),
    inverse: (fraction: number) => screenToWall(fraction, 1, scale),
  };
  const radialTicks = wallTicks(1, scale).map(value => ({ value, label: String(value) }));
  const selected = (position: number, distance: number) => { setZ(position); setY(distance); };
  return <section className="site-reference" aria-label="Site hollow-fibre profiles">
    <section className="card" data-site-cutaway="true">
      <div className="card-heading">
        <h2>Site hollow-fibre reference</h2>
        <label>Radial spacing <select aria-label="Site radial spacing" value={scale}
          onChange={event => setScale(event.target.value as RadialScale)}>
          <option value="wall-log">Wall-log</option><option value="linear">Linear</option>
        </select></label>
      </div>
      <p className="scope-hint">Precomputed Python preview at {data.source.nominalTmp} kPa,
        particle radius {data.source.radiusNm} nm and φ<sub>b</sub> = {data.source.bulkFraction}.
        This preview is separate from the C++ regression.</p>
      <FiltrationSchematic geometry="HF" result={data} position={z} wallDistance={y} radialScale={scale}
        cpRatio={cp?.ratio} onSelect={selected} />
      <label className="slider">Axial position z/L <input aria-label="Site axial position" type="range"
        min="0" max="1" step="0.001" value={z} onChange={event => setZ(+event.target.value)} /></label>
      <label className="slider">Wall distance y/R <input aria-label="Site wall distance" type="range"
        min="0" max="1" step="0.001" value={y} onChange={event => setY(+event.target.value)} /></label>
      <p className="coordinate" aria-live="polite">z/L = {z.toPrecision(5)} / <RadialCoordinates y={y} ratio={cp?.ratio} /></p>
      {cp && <p className="scope-hint">Fixed reference δ<sub>CP</sub> = {cp.delta.toExponential(6)} m.
        The dashed line marks one reference CP length from the membrane.</p>}
    </section>
    <section className="card">
      <div className="card-heading">
        <h2>Site axial and radial profiles</h2>
        <label>Velocity component <select aria-label="Site velocity component" value={velocity}
          onChange={event => setVelocity(event.target.value as typeof velocity)}>
          <option value="v_m_s">Transverse v</option><option value="u_m_s">Axial u</option>
        </select></label>
      </div>
      <div className="profile-grid" data-site-profiles="true">
        <LinePlot points={axialProfileAtY(data, 'phi', y, data.length_m)} label="Axial concentration φ"
          xLabel="Axial position z/L" selectedX={z} onSelectX={setZ} height={255} />
        <LinePlot points={profileAtZ(data, 'phi', z, data.length_m)} label="Radial concentration φ"
          xLabel="Distance from membrane y/R" xScale={radialScale} xTicks={radialTicks}
          selectedX={y} onSelectX={setY} height={255} />
        <LinePlot points={axialProfileAtY(data, velocity, y, data.length_m).map(([at, value]) => [at, value * multiplier])}
          label={`Axial profile of ${velocityLabel.toLowerCase()}`} xLabel="Axial position z/L"
          selectedX={z} onSelectX={setZ} height={255} />
        <LinePlot points={profileAtZ(data, velocity, z, data.length_m).map(([at, value]) => [at, value * multiplier])}
          label={`Radial profile of ${velocityLabel.toLowerCase()}`} xLabel="Distance from membrane y/R"
          xScale={radialScale} xTicks={radialTicks} selectedX={y} onSelectX={setY} height={255} />
      </div>
      <p className="scope-hint">Cutaway colour uses the site's 24 × 80 adaptive concentration cells.
        Charts sample its saved 24 × 320 display profiles in physical coordinates. No solver is run.</p>
      <details><summary>Site preset provenance</summary><p>{data.source.id} / Python source revision{' '}
        <code>{data.source.commit}</code>. Preset SHA-256 <code>{data.source.sha256}</code>.</p></details>
    </section>
  </section>;
}
