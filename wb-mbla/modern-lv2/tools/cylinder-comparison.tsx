// Adapted from the 2026-11-25 mBLA cylinder regression renderer under Gunwoo Park’s direction.
import { useMemo, useState } from 'react';
import { createRoot } from 'react-dom/client';
import { FiltrationSchematic } from './schematic';
import { LinePlot } from './plots';
import { SiteReference } from './site-reference';
import { RadialCoordinates } from './radial-coordinates';
import { screenToWall, wallTicks, wallToScreen, type RadialScale } from './coordinates';

type Case = {
  id: string; label: string; coverage: number; upper: number;
  meanTmp?: number; meanPermeate: number; checked: boolean;
  archiveMeanPermeate?: number; level2MeanPermeate?: number;
  lengthMm?: number; archive?: number[][]; recomputed?: number[][];
  source?: Record<string, string>;
  outletArchive?: number[][]; outletRecomputed?: number[][];
};
type SearchSample = {
  search: number; trials: number; resistance: {cpp: number; python: number; numba: number};
};
type CakeWorkload = {
  case: string; radial_nodes: number; searches: number; trial_count: number; cpp_seconds: number;
  python: {median_seconds: number; max_residual_abs_difference: number};
  numba: {median_seconds: number; max_residual_abs_difference: number};
  samples: SearchSample[];
};
type Payload = {
  schema: string; initial: string; quantities: string[]; cases: Case[];
  cakeSearch?: {cases: CakeWorkload[]; numba_first_call_seconds: number};
};
const schema = 'wb.mbla.cylinder-comparison.v3';
const displayLabel = (label: string) => label.replaceAll(' \u00b7 ', ' / ');
const format = (value: number) => String(value);

type Decimal = { integer: bigint; power: number };
const zero: Decimal = { integer: 0n, power: 0 };

// v3 stores the source's printed decimals. Align their decimal scales before
// subtraction so equally sized changes tie at the first stored station.
function decimal(value: number): Decimal {
  const [mantissa, exponent = '0'] = String(value).split('e');
  const fraction = mantissa.split('.')[1]?.length ?? 0;
  return { integer: BigInt(mantissa.replace('.', '')), power: Number(exponent) - fraction };
}

function subtract(left: Decimal, right: Decimal): Decimal {
  const power = Math.min(left.power, right.power);
  const integer = left.integer * 10n ** BigInt(left.power - power) -
    right.integer * 10n ** BigInt(right.power - power);
  return { integer, power };
}

function storedDelta(archive: number, replay: number) {
  return archive === replay ? zero : subtract(decimal(replay), decimal(archive));
}

function formatDelta({ integer, power }: Decimal) {
  if (!integer) return '0';
  const negative = integer < 0n;
  const digits = String(negative ? -integer : integer);
  const exponent = digits.length + power - 1;
  const significant = digits.replace(/0+$/, '');
  const sign = negative ? '-' : '';
  if (exponent < -4 || exponent >= 16) {
    const fraction = significant.slice(1);
    return sign + significant[0] + (fraction && '.' + fraction) + 'e' + exponent;
  }
  const point = exponent + 1;
  if (point <= 0) return sign + '0.' + '0'.repeat(-point) + significant;
  if (point >= significant.length) return sign + significant.padEnd(point, '0');
  return sign + significant.slice(0, point) + '.' + significant.slice(point);
}

function nearest(rows: number[][], value: number) {
  return rows.reduce((best, row, i) =>
    Math.abs(row[0] - value) < Math.abs(rows[best][0] - value) ? i : best, 0);
}

function storedDifferences(archive: number[][], replay: number[][], quantities: string[]) {
  return quantities.slice(1).flatMap((quantity, offset) => {
    const column = offset + 1;
    let count = 0;
    let maximum = zero;
    let station = 0;
    for (const [index, row] of archive.entries()) {
      const delta = storedDelta(row[column], replay[index][column]);
      if (!delta.integer) continue;
      count++;
      const magnitude = { ...delta, integer: delta.integer < 0n ? -delta.integer : delta.integer };
      if (subtract(magnitude, maximum).integer > 0n) {
        maximum = magnitude;
        station = row[0];
      }
    }
    return count ? [{ quantity, count, maximum, station }] : [];
  });
}

function StoredValues({ rows, label }: {
  rows: (readonly [string, number, number])[];
  label: string;
}) {
  return <div className="table-wrap"><table aria-label={label}>
    <thead><tr><th>Observable</th><th>Archived</th><th>Level 1</th><th>Level 2</th><th>Level 2 − archive</th></tr></thead>
    <tbody>{rows.map(([name, before, after]) => <tr key={name}>
      <th scope="row">{name}</th><td>{format(before)}</td><td>{format(after)}</td><td>{format(after)}</td>
      <td>{formatDelta(storedDelta(before, after))}</td>
    </tr>)}</tbody>
  </table></div>;
}

function plotSeries(archive: number[][], replay: number[][], column: number) {
  // The original receipt identifies identical Level 1 and Level 2 tables.
  const points = (rows: number[][]) => rows.map(row => [row[0], row[column]]);
  const replayPoints = points(replay);
  return [
    {label: 'Level 1', points: replayPoints, className: 'recomputed'},
    {label: 'Archived', points: points(archive), className: 'archive'},
    {label: 'Level 2', points: replayPoints, className: 'level2'},
  ];
}

type CPScale = { delta_m: number; ratio: number };

function Comparison({ item, quantities, cp }: { item: Case; quantities: string[]; cp?: CPScale }) {
  const [column, setColumn] = useState(1);
  const [index, setIndex] = useState(0);
  const [radialIndex, setRadialIndex] = useState(0);
  const [scale, setScale] = useState<RadialScale>('wall-log');
  const archive = item.archive!;
  const recomputed = item.recomputed!;
  const outlet = item.outletArchive!;
  const outletNew = item.outletRecomputed!;
  const station = recomputed[index][0];
  const distance = outletNew[radialIndex][0];
  const boundaryProfiles = useMemo(() => ({
    wall: recomputed.map(row => [row[0], row[1]]), outlet: outletNew,
  }), [recomputed, outletNew]);
  const radialTransform = {
    forward: (value: number) => wallToScreen(value, 1, scale),
    inverse: (fraction: number) => screenToWall(fraction, 1, scale),
  };
  const selectAxial = (z: number) => setIndex(nearest(recomputed, z));
  const selectRadial = (y: number) => setRadialIndex(nearest(outletNew, y));
  const axialDifferences = useMemo(() => storedDifferences(archive, recomputed, quantities),
    [archive, recomputed, quantities]);

  return <div className="workspace" data-regression-profiles="true">
    <section className="card profile-cutaway">
      <div className="card-heading"><h2>Hollow-fibre selector</h2>
        <label>Radial spacing <select aria-label="Radial spacing" value={scale} onChange={event => setScale(event.target.value as RadialScale)}>
          <option value="wall-log">Wall-log</option><option value="linear">Linear</option>
        </select></label></div>
      <FiltrationSchematic geometry="HF" position={station} wallDistance={distance} radialScale={scale}
        cpRatio={cp?.ratio} boundaryProfiles={boundaryProfiles} onSelect={(z, y) => { selectAxial(z); selectRadial(y); }} />
      <p className="scope-hint">Log colours show stored replay concentration along the wall and at the outlet.
        The radial probe always reads the outlet section at z/L = 1.</p>
      {cp && <p className="scope-hint">Fixed reference δ<sub>CP</sub> = {cp.delta_m.toExponential(6)} m.
        The dashed line marks one reference CP length from the membrane.</p>}
    </section>
    <div className="profile-grid">
      <section className="card" data-axial-profile="true">
        <div className="card-heading"><h2>Axial profile</h2>
          <label>Quantity <select aria-label="Comparison quantity" value={column}
            onChange={event => setColumn(+event.target.value)}>
            {quantities.slice(1).map((name, i) => <option key={name} value={i + 1}>{name}</option>)}
          </select></label>
        </div>
        <LinePlot series={plotSeries(archive, recomputed, column)} label={quantities[column]}
          xLabel="Axial position z/L" selectedX={station} onSelectX={selectAxial} height={255} />
        <label className="slider">Stored axial station
          <input type="range" min="0" max={recomputed.length - 1} value={index}
            aria-label="Stored axial station"
            aria-valuetext={`Station ${index + 1} of ${recomputed.length}, z/L = ${format(station)}`}
            onChange={event => setIndex(+event.target.value)} />
        </label>
        <p className="coordinate" aria-live="polite">Station {index + 1} of {recomputed.length} /{' '}
          <output aria-label="Selected axial location">z/L = {format(station)} / z = {(station * item.lengthMm!).toFixed(2)} mm</output>
        </p>
        <StoredValues label="Stored axial values" rows={quantities.slice(1).map((name, i) =>
          [name, archive[index][i + 1], recomputed[index][i + 1]])} />
        <details>
          <summary>Whole-axis regression differences</summary>
          {axialDifferences.length ? <div className="table-wrap"><table aria-label="Whole-axis regression differences">
            <thead><tr><th>Observable</th><th>Max |replay − archive|</th><th>Changed stations</th><th>First maximum z/L</th></tr></thead>
            <tbody>{axialDifferences.map(({ quantity, count, maximum, station: first }) => <tr key={quantity}>
              <th scope="row">{quantity}</th><td>{formatDelta(maximum)}</td><td>{count}</td><td>{format(first)}</td>
            </tr>)}</tbody>
          </table></div> : <p className="scope-hint">All shared axial columns agree at the stored precision.</p>}
        </details>
      </section>
      <section className="card" data-outlet-profile="true">
        <div className="card-heading"><h2>Outlet radial profile</h2><span>z/L = 1</span></div>
        <LinePlot series={plotSeries(outlet, outletNew, 1)} label="Outlet volume fraction φ"
          xLabel="Distance from membrane y/R" xScale={radialTransform}
          xTicks={wallTicks(1, scale).map(value => ({value, label: String(value)}))}
          selectedX={distance} onSelectX={selectRadial} height={255} />
        <label className="slider">Stored radial sample
          <input type="range" min="0" max={outletNew.length - 1} value={radialIndex}
            aria-label="Stored radial sample"
            aria-valuetext={`Sample ${radialIndex + 1} of ${outletNew.length}, y/R = ${format(distance)}`}
            onChange={event => setRadialIndex(+event.target.value)} />
        </label>
        <p className="coordinate" aria-live="polite">Sample {radialIndex + 1} of {outletNew.length} /{' '}
          <output aria-label="Selected outlet location">z/L = 1 / <RadialCoordinates y={distance} ratio={cp?.ratio} exact /></output>
        </p>
        <StoredValues label="Stored outlet values" rows={[
          ['Outlet volume fraction φ', outlet[radialIndex][1], outletNew[radialIndex][1]],
        ]} />
        <p className="scope-hint">The outlet is the stored last-iteration sidecar before final axial analysis.
          The saved Level 1 and Level 2 digests identify the same axial and outlet tables.</p>
      </section>
    </div>
  </div>;
}

function CaseList({ cases, onSelect }: { cases: Payload; onSelect: (id: string) => void }) {
  return <section className="card result-list" id="case-list" aria-labelledby="case-list-title">
    <h2 id="case-list-title">Case results</h2>
    <p className="scope-hint">Mean permeate velocity in µm/s. Select a case to inspect its full axial and outlet results.</p>
    <div className="table-wrap"><table aria-label="Mean permeate velocity by case">
      <thead><tr><th scope="col">Case / nominal TMP</th><th scope="col">Mean TMP (kPa)</th>
        <th scope="col">Archived</th><th scope="col">Level 1</th><th scope="col">Level 2</th>
        <th scope="col">Level 1–2 axial match</th></tr></thead>
      <tbody>{cases.cases.filter(item => item.checked).map(item => <tr key={item.id}>
        <th scope="row"><a href="#case-detail" onClick={() => onSelect(item.id)}>{displayLabel(item.label)}</a></th>
        <td>{item.meanTmp!.toFixed(6)}</td><td>{item.archiveMeanPermeate!.toFixed(9)}</td>
        <td>{item.meanPermeate.toFixed(9)}</td><td>{item.level2MeanPermeate!.toFixed(9)}</td>
        <td>{item.recomputed!.length} / {item.recomputed!.length}</td>
      </tr>)}</tbody>
    </table></div>
  </section>;
}

const workloadLabels: Record<string, string> = {
  mild: 'Mild cake', demanding: 'Demanding cake', 'wide-mild': 'Mild cake',
};

function CakeSearch({ search }: { search: NonNullable<Payload['cakeSearch']> }) {
  const [id, setId] = useState('demanding');
  const workload = search.cases.find(row => row.case === id)!;
  const matchingSearches = search.cases.reduce((count, row) => count + row.searches, 0);
  const trials = search.cases.reduce((count, row) => count + row.trial_count, 0);
  const ms = (seconds: number) => (seconds * 1000).toFixed(3);
  const name = (row: CakeWorkload) => `${workloadLabels[row.case]} / ${row.radial_nodes} radial nodes`;
  return <section className="card result-list" id="cake-search" aria-labelledby="cake-search-title">
    <h2 id="cake-search-title">Numba inner cake-resistance benchmark</h2>
    <p className="scope-hint">Level 2 C++, Python and Numba matched all {matchingSearches} captured searches and {trials.toLocaleString('en-US')} trials. These measurements cover the inner search. The axial comparison above uses the full C++ solver.</p>
    <div className="table-wrap"><table aria-label="Inner search benchmark timings in milliseconds">
      <thead><tr><th scope="col">Workload</th><th scope="col">Searches</th><th scope="col">Trials</th>
        <th scope="col">Level 2 C++ (ms)</th><th scope="col">Python (ms)</th><th scope="col">Numba (ms)</th>
        <th scope="col">Max |Numba − C++ residual|</th></tr></thead>
      <tbody>{search.cases.map(row => <tr key={row.case}>
        <th scope="row">{name(row)}</th><td>{row.searches}</td><td>{row.trial_count}</td>
        <td>{ms(row.cpp_seconds)}</td><td>{ms(row.python.median_seconds)}</td><td>{ms(row.numba.median_seconds)}</td>
        <td>{row.numba.max_residual_abs_difference.toExponential(2)}</td>
      </tr>)}</tbody>
    </table></div>
    <p className="scope-hint">Median batch times after warm-up. First Numba call including compilation took {search.numba_first_call_seconds.toFixed(2)} s.</p>
    <details>
      <summary>Compare returned resistance Rc/Rm by search</summary>
      <label>Workload <select aria-label="Cake search workload" value={id} onChange={event => setId(event.target.value)}>
        {search.cases.map(row => <option key={row.case} value={row.case}>{name(row)}</option>)}
      </select></label>
      <div className="table-wrap"><table aria-label="Returned cake resistance Rc/Rm">
        <thead><tr><th scope="col">Captured search</th><th scope="col">Trials</th><th scope="col">Level 2 C++</th>
          <th scope="col">Python</th><th scope="col">Numba</th></tr></thead>
        <tbody>{workload.samples.map(sample => <tr key={sample.search}>
          <th scope="row">{sample.search}</th><td>{sample.trials}</td><td>{format(sample.resistance.cpp)}</td>
          <td>{format(sample.resistance.python)}</td><td>{format(sample.resistance.numba)}</td>
        </tr>)}</tbody>
      </table></div>
    </details>
  </section>;
}

function App({ cases, cpScales }: { cases: Payload; cpScales: Record<string, CPScale> }) {
  const [id, setId] = useState(cases.initial);
  const [view, setView] = useState<'regression' | 'site'>('regression');
  const item = cases.cases.find(entry => entry.id === id)!;
  return <>
    <nav className="view-switch" aria-label="Profile data view">
      <button type="button" aria-pressed={view === 'regression'} onClick={() => setView('regression')}>Validation comparison</button>
      <button type="button" aria-pressed={view === 'site'} onClick={() => setView('site')}>Site hollow-fibre reference</button>
    </nav>
    {view === 'site' ? <SiteReference /> : <>
    <p className="scope-hint">Stored archive and replay comparison. The inlet participates in numerical regression; physical entrance development is outside this model.</p>
    <section className="card overview" id="case-detail">
      <label>Case <select aria-label="Case" value={id} onChange={event => setId(event.target.value)}>
        <optgroup label="Full axial comparisons">{cases.cases.filter(entry => entry.checked).map(entry =>
          <option value={entry.id} key={entry.id}>{displayLabel(entry.label)}</option>)}</optgroup>
        <optgroup label="Summary only / over half covered">{cases.cases.filter(entry => !entry.checked).map(entry =>
          <option value={entry.id} key={entry.id}>{displayLabel(entry.label)}</option>)}</optgroup>
      </select></label>
      <dl>
        {item.checked && <div><dt>Mean TMP</dt><dd>{item.meanTmp!.toFixed(6)} kPa</dd></div>}
        <div><dt>Mean permeate velocity</dt><dd>{item.meanPermeate.toFixed(9)} µm/s</dd></div>
        <div><dt>Estimated cake coverage</dt><dd>{item.coverage.toFixed(2)}%{item.upper !== item.coverage && ('–' + item.upper.toFixed(2) + '%')}</dd></div>
      </dl>
    </section>
    {item.checked ? <Comparison key={item.id} item={item} quantities={cases.quantities} cp={cpScales[item.id]} /> :
      <section className="card" data-report-only="true"><h2>Summary only / extensive cake</h2>
        <p>The upper estimate of cake coverage exceeds half of the membrane. The influence of this coverage may be material.</p></section>}
    <CaseList cases={cases} onSelect={setId} />
    {cases.cakeSearch && <CakeSearch search={cases.cakeSearch} />}
    </>}
  </>;
}

const axialCount = 200;
const outletCount = 400;
const axialQuantities = [
  'Axial position z/L', 'Wall volume fraction φw', 'TMP (kPa)', 'Permeate velocity (µm/s)',
  'Axial particle flux (m³/s)', 'Cake resistance Rc/Rm', 'Bulk pressure (Pa)',
  'Axial velocity scale (m/s)', 'Wall osmotic pressure (Pa)',
];
const isText = (value: unknown): value is string =>
  typeof value === 'string' && value.trim().length > 0 && value.length <= 240;
const finite = (value: unknown): value is number => typeof value === 'number' && Number.isFinite(value);
const validDigest = (value: unknown): value is string =>
  typeof value === 'string' && /^[a-f0-9]{64}$/.test(value);

function validatePairedRows(original: unknown, replay: unknown, count: number, columns: number) {
  if (!Array.isArray(original) || !Array.isArray(replay) ||
      original.length !== count || replay.length !== count) return false;
  return original.every((row, i) =>
    Array.isArray(row) && row.length === columns && row.every(finite) &&
    Array.isArray(replay[i]) && replay[i].length === columns && replay[i].every(finite) &&
    row[0] === replay[i][0] && row[0] >= 0 && row[0] <= 1 &&
    (i === 0 ? row[0] === 0 : row[0] > original[i - 1][0]));
}

function validateBenchmark(search: Payload['cakeSearch']) {
  if (search === undefined) return;
  const nonnegative = (value: unknown) => finite(value) && value >= 0;
  if (!search || !nonnegative(search.numba_first_call_seconds) ||
      !Array.isArray(search.cases) || !search.cases.length || search.cases.length > 32 ||
      !search.cases.some(row => row?.case === 'demanding')) {
    throw new Error('Invalid saved cake-search benchmark.');
  }
  const seen = new Set();
  for (const row of search.cases) {
    if (!row || !Object.hasOwn(workloadLabels, row.case) || seen.has(row.case) ||
        !Number.isInteger(row.radial_nodes) || row.radial_nodes <= 0 ||
        !Number.isInteger(row.searches) || !Number.isInteger(row.trial_count) ||
        !nonnegative(row.cpp_seconds) || !nonnegative(row.python?.median_seconds) ||
        !nonnegative(row.numba?.median_seconds) || !nonnegative(row.numba?.max_residual_abs_difference) ||
        !Array.isArray(row.samples) || row.samples.length !== row.searches || row.samples.length > 1000 ||
        !row.samples.every((sample, index) => sample && Number.isInteger(sample.search) &&
          sample.search > (index ? row.samples[index - 1].search : 0) &&
          Number.isInteger(sample.trials) && sample.trials > 0 &&
          finite(sample.resistance?.cpp) && sample.resistance.cpp === sample.resistance.python &&
          sample.resistance.cpp === sample.resistance.numba) ||
        row.samples.reduce((total, sample) => total + sample.trials, 0) !== row.trial_count) {
      throw new Error('Invalid saved cake-search workload.');
    }
    seen.add(row.case);
  }
}

function validateComparison(data: Payload) {
  if (!data || data.schema !== schema || !isText(data.initial) ||
      !Array.isArray(data.quantities) || data.quantities.length !== axialQuantities.length ||
      !axialQuantities.every((quantity, index) => quantity === data.quantities[index]) ||
      !Array.isArray(data.cases) || !data.cases.length || data.cases.length > 32) {
    throw new Error('Expected a cylinder comparison with nine recorded columns and a case list.');
  }
  const ids = new Set();
  for (const item of data.cases) {
    if (!item || !isText(item.id) || !isText(item.label) || ids.has(item.id) ||
        typeof item.checked !== 'boolean' || !finite(item.coverage) || !finite(item.upper) ||
        item.coverage < 0 || item.upper < item.coverage || item.upper > 100 ||
        !finite(item.meanPermeate) || (item.meanTmp !== undefined && !finite(item.meanTmp))) {
      throw new Error('Invalid cylinder case summary.');
    }
    ids.add(item.id);
    if (!item.checked) {
      if (item.upper <= 50 || [item.archive, item.recomputed, item.outletArchive, item.outletRecomputed]
          .some(rows => rows !== undefined)) {
        throw new Error('An extensive-cake case must be a summary without plotted tables.');
      }
      continue;
    }
    if (item.upper > 50 || !finite(item.lengthMm) || item.lengthMm <= 0 || !finite(item.meanTmp) ||
        !finite(item.archiveMeanPermeate) || item.level2MeanPermeate !== item.meanPermeate ||
        !validatePairedRows(item.archive, item.recomputed, axialCount, axialQuantities.length) ||
        item.archive![axialCount - 1][0] !== 1 ||
        !validatePairedRows(item.outletArchive, item.outletRecomputed, outletCount, 2)) {
      throw new Error('Expected 200 paired axial rows from z/L = 0 to 1 and 400 paired outlet samples.');
    }
    const sources = item.source;
    if (!sources || ![sources.archiveSha256, sources.recomputedSha256, sources.level2Sha256,
                      sources.archiveOutletSha256, sources.recomputedOutletSha256,
                      sources.level2OutletSha256].every(validDigest) ||
        sources.level2Sha256 !== sources.recomputedSha256 ||
        sources.level2OutletSha256 !== sources.recomputedOutletSha256) {
      throw new Error('The saved Level 1 and Level 2 table digests must identify the same tables.');
    }
  }
  if (!ids.has(data.initial)) throw new Error('Initial cylinder case is missing from the data.');
  validateBenchmark(data.cakeSearch);
}

declare global {
  interface Window {
    wbMountCylinderComparison: (element: HTMLElement, data: Payload | null, savedReference?: boolean) => void;
  }
}

let root: ReturnType<typeof createRoot> | undefined;
let importSequence = 0;
window.wbMountCylinderComparison = (element, data, savedReference) => {
  if (data) {
    validateComparison(data);
    // A separately imported file can reuse case IDs without sharing the reference inputs.
    const cpScales = savedReference
      ? JSON.parse(document.getElementById('validation-cp-scales')!.textContent!)
      : {};
    root ??= createRoot(element);
    root.render(<App key={++importSequence} cases={data} cpScales={cpScales} />);
  } else {
    root?.render(null);
  }
};
