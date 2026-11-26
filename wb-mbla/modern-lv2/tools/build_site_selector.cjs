// Regenerate the standalone explorer from the existing research-site geometry.
const fs = require('node:fs');
const os = require('node:os');
const path = require('node:path');
const { createHash } = require('node:crypto');

const site = process.argv[2];
const esbuildModule = process.argv[3];
if (!site || !esbuildModule) {
  throw new Error('Usage: node build_site_selector.cjs RESEARCH_SITE ESBUILD_MODULE');
}
const siteRoot = path.resolve(site);
const source = path.join(siteRoot, 'app/interactive-tools/filtration');
const modules = path.join(siteRoot, 'node_modules');
const esbuild = require(path.resolve(esbuildModule));
const selectorFiles = ['schematic.tsx', 'coordinates.ts', 'plots.tsx', 'model.ts', 'sampling.ts', 'cp-scale.ts'];
const upstreamHashes = {};
const digest = (value) => createHash('sha256').update(value).digest('hex');
const begin = '<!-- research-site selector start -->';
const end = '<!-- research-site selector end -->';
const referenceTag = '<script id="profile-data" type="application/json">';

// Embed exactly the recorded regression JSON, including the original case order and initial case.
const report = path.join(__dirname, '../docs/validation-report');
const reference = fs.readFileSync(path.join(report, 'CUSTOM_DATA.json'), 'utf8');
const receipt = JSON.parse(fs.readFileSync(path.join(report, 'VIEWER_PROVENANCE.json'), 'utf8'));
if (digest(reference) !== receipt.embedded_data_sha256 || reference.includes('<')) {
  throw new Error('The recorded regression payload failed its digest or HTML embedding check.');
}

// Export only the public preset fields needed by the site plots and CP reference.
const presets = path.join(siteRoot, 'public/tools/mbla/presets');
const rawSite = fs.readFileSync(path.join(presets, 'HF.json'), 'utf8');
const listedSite = JSON.parse(fs.readFileSync(path.join(presets, 'manifest.json'), 'utf8')).cases
  .find(item => item.geometry === 'HF' && item.file === 'HF.json');
if (!listedSite || listedSite.sha256 !== digest(rawSite) || listedSite.bytes !== Buffer.byteLength(rawSite)) {
  throw new Error('The selected site HF preset differs from its manifest.');
}
const hf = JSON.parse(rawSite);
const siteData = {
  source: {
    id: hf.preset.id, sha256: listedSite.sha256, commit: hf.commit,
    nominalTmp: hf.parameters.tmp_kpa, radiusNm: hf.parameters.radius_nm, bulkFraction: hf.parameters.phi_bulk,
  },
  geometry: { code: hf.geometry.code }, length_m: hf.length_m, rows: hf.rows.map(row => [row[0]]),
  fields: { y_over_R: hf.fields.y_over_R, phi: hf.fields.phi, radius_m: hf.fields.radius_m },
  profiles: {
    y_over_R: hf.profiles.y_over_R, phi: hf.profiles.phi, u_m_s: hf.profiles.u_m_s, v_m_s: hf.profiles.v_m_s,
  },
  properties: { D0_m2_s: hf.properties.D0_m2_s },
  boundary_conditions: {
    solvent_outlet_pressure_pa: hf.boundary_conditions.solvent_outlet_pressure_pa,
    permeate_pressure_pa: hf.boundary_conditions.permeate_pressure_pa,
    reference_axial_pressure_drop_pa: hf.boundary_conditions.reference_axial_pressure_drop_pa,
  },
};
const numbers = value => Array.isArray(value) && value.length > 0 && value.every(Number.isFinite);
const axialCount = siteData.rows.length;
function validGrid(coordinates, fields) {
  return numbers(coordinates) && coordinates[0] === 0 && coordinates.at(-1) === 1 &&
    coordinates.every((value, index) => !index || value > coordinates[index - 1]) && fields.every(field =>
      Array.isArray(field) && field.length === axialCount &&
      field.every(row => numbers(row) && row.length === coordinates.length));
}
if (siteData.geometry.code !== 'HF' || !numbers(siteData.rows.map(([z]) => z)) ||
    !validGrid(siteData.fields.y_over_R, [siteData.fields.phi]) ||
    !validGrid(siteData.profiles.y_over_R, [siteData.profiles.phi, siteData.profiles.u_m_s, siteData.profiles.v_m_s])) {
  throw new Error('Unexpected field geometry in the saved site preview.');
}
const cpScales = fs.readFileSync(path.join(__dirname, 'validation-cp-scales.json'), 'utf8').trim();

const scratch = fs.mkdtempSync(path.join(os.tmpdir(), 'wb-filtration-selector-'));
try {
  for (const name of selectorFiles) {
    const bytes = fs.readFileSync(path.join(source, name));
    upstreamHashes[name] = digest(bytes);
    fs.writeFileSync(path.join(scratch, name), bytes);
  }
  for (const name of ['site-selector.tsx', 'cylinder-comparison.tsx', 'boundary-shading.tsx', 'site-reference.tsx', 'radial-coordinates.tsx']) {
    fs.copyFileSync(path.join(__dirname, name), path.join(scratch, name));
  }
  fs.writeFileSync(path.join(scratch, 'entry.tsx'),
    'import "./site-selector";\nimport "./cylinder-comparison";\n');

  // The v3 regression has boundary profiles. Keep the site's full-field renderer for its own case.
  const schematic = path.join(scratch, 'schematic.tsx');
  const replacements = [
    ["import type { PointerEvent } from 'react';", "import type { PointerEvent } from 'react';\nimport { BoundaryShading, type BoundaryProfiles } from './boundary-shading';"],
    ['  result?: Result;', "  result?: Pick<Result, 'rows' | 'length_m'> & { geometry: Pick<Result['geometry'], 'code'>; fields: Pick<Result['fields'], 'phi' | 'y_over_R'> };\n  boundaryProfiles?: BoundaryProfiles;"],
    ['result, cpRatio, onSelect }: SchematicProps', 'result, boundaryProfiles, cpRatio, onSelect }: SchematicProps'],
    ['          {!hasData && <g>', '          {boundaryProfiles && <BoundaryShading profiles={boundaryProfiles} x={x} y={y} />}\n          {!hasData && !boundaryProfiles && <g>'],
    [": 'The particle layer is an illustration while data are loading.'", ": boundaryProfiles ? 'Colour shows the stored wall and outlet concentration only.' : 'The particle layer is an illustration while data are loading.'"],
  ];
  let adapted = fs.readFileSync(schematic, 'utf8');
  for (const [from, to] of replacements) {
    if (!adapted.includes(from)) throw new Error('The site schematic changed. Check the adapter before rebuilding.');
    adapted = adapted.replace(from, to);
  }
  fs.writeFileSync(schematic, adapted.replaceAll(' \u00b7 ', ' - '));

  // Narrow the compile-time signatures to the fields the original implementations read.
  const cpHelper = path.join(scratch, 'cp-scale.ts');
  fs.writeFileSync(cpHelper, fs.readFileSync(cpHelper, 'utf8').replace(
    'result: Result',
    "result: { boundary_conditions: Pick<NonNullable<Result['boundary_conditions']>, 'solvent_outlet_pressure_pa' | 'permeate_pressure_pa' | 'reference_axial_pressure_drop_pa'>; properties: Pick<Result['properties'], 'D0_m2_s'>; fields: Pick<Result['fields'], 'radius_m'> }",
  ));
  const sampling = path.join(scratch, 'sampling.ts');
  fs.writeFileSync(sampling, fs.readFileSync(sampling, 'utf8')
    .replace("import type { Comparison, FieldKey, Result } from './model';", "type FieldKey = 'phi' | 'u_m_s' | 'v_m_s';")
    .replace("type Solution = Pick<Result | Comparison, 'rows' | 'profiles'>;",
      "type Solution = { rows: number[][]; profiles: { y_over_R: number[] } & Record<FieldKey, number[][]> };"));


  const result = esbuild.buildSync({
    entryPoints: [path.join(scratch, 'entry.tsx')],
    nodePaths: [modules], bundle: true, minify: true, write: false,
    format: 'iife', platform: 'browser', jsx: 'automatic',
    define: {'process.env.NODE_ENV': '"production"'}, legalComments: 'inline',
  });
  const javascript = result.outputFiles[0].text.replaceAll('</script', '<\\/script');
  const packages = ['react', 'react-dom', 'scheduler'];
  const licenses = [...new Set(packages.map(name =>
    fs.readFileSync(path.join(modules, name, 'LICENSE'), 'utf8')
  ))].join('\n');
  const explorer = path.join(__dirname, 'profile-explorer.html');
  const html = fs.readFileSync(explorer, 'utf8');
  const from = html.indexOf(begin), through = html.indexOf(end);
  if (from < 0 || through < from || html.indexOf(begin, from + 1) >= 0) {
    throw new Error('Expected a single research-site selector block.');
  }
  const block = [
    begin,
    '<!-- Research-site FiltrationSchematic and LinePlot. Original project direction: Gunwoo Park. -->',
    '<style>',
    fs.readFileSync(path.join(__dirname, 'cylinder-comparison.css'), 'utf8').trimEnd(),
    '</style>',
    `<script id="research-site-selector" data-schematic-sha256="${upstreamHashes['schematic.tsx']}">\n${javascript}</script>`,
    '<!-- Bundled React, React DOM and Scheduler license texts',
    licenses,
    '-->',
    end,
  ].join('\n');
  const withSelector = html.slice(0, from) + block + html.slice(through + end.length);
  const tag = withSelector.indexOf(referenceTag);
  const contentStart = tag + referenceTag.length;
  const contentEnd = withSelector.indexOf('</script>', contentStart);
  if (tag < 0 || contentEnd < 0 || withSelector.indexOf(referenceTag, contentStart) >= 0) {
    throw new Error('Expected a single reference data block.');
  }
  let output = withSelector.slice(0, contentStart) + reference + withSelector.slice(contentEnd);
  for (const [id, data] of [['site-reference-data', JSON.stringify(siteData)], ['validation-cp-scales', cpScales]]) {
    if (data.includes('<')) throw new Error(`The selected ${id} failed its HTML embedding check.`);
    const marker = `<script id="${id}" type="application/json">`;
    const existing = output.indexOf(marker);
    if (existing >= 0) {
      const close = output.indexOf('</script>', existing);
      output = output.slice(0, existing) + output.slice(close + '</script>'.length);
    }
    output = output.replace(referenceTag, `${marker}${data}</script>\n${referenceTag}`);
  }
  fs.writeFileSync(explorer, output);
  console.log(JSON.stringify({
    upstreamHashes, esbuild: esbuild.version, bundle: digest(javascript), reference: digest(reference),
    siteReferenceSource: listedSite.sha256, siteReferenceSlice: digest(JSON.stringify(siteData)), cpScales: digest(cpScales),
    dependencies: Object.fromEntries(packages.map(name => [name,
      JSON.parse(fs.readFileSync(path.join(modules, name, 'package.json'), 'utf8')).version
    ])),
  }, null, 2));
} finally {
  fs.rmSync(scratch, {recursive: true, force: true});
}
