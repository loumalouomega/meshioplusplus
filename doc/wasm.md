# WebAssembly / JavaScript

## Exodus series and pipeline Version 2

`createExodusTimeSeriesWriter(path)` returns a handle wrapper with `writePointsCells(mesh)`, `writeData(time, mesh)`, `flush`, `finalize`, `numSteps`, `finalized` and idempotent `close`. It uses MEMFS, fixes geometry/sets/attributes once and retains no step history. Shipped netCDF-enabled artifacts support it; custom builds without netCDF throw a named dependency error. Invalid/stale native handles fail rather than accessing freed memory. See [Exodus](./formats/exodus.md#stateful-series-writing).

`runPipeline` accepts [Version 2 spatial documents](./pipeline.md#version-2-spatial-multi-mesh-steps) with per-step file `Inputs` and terminal Split/Partition `Output.Pattern` (`{key}`/`{part}`), all in MEMFS. Version 1, `convertSurfaceOps`/`surfaceBuffersOps` and the browser viewer's per-mesh chip pipeline are unchanged; spatial branches and transient sequence schemas are separate.

## Pipeline reports

`runPipeline` already returns the shared structured report object (`steps`, `warnings`). The new C/Fortran/Julia/R caller-buffer report APIs use the same names and PascalCase counters; non-finite counters serialize as JSON `null`. See [pipeline reports](./pipeline.md#structured-reports-on-the-flat-abi).

For a sequential-only local artifact, run the complete WASM smoke suite with `MESHIOPLUSPLUS_WASM_VARIANT=seq node tests/wasm/smoke.mjs`; the default still exercises the threaded artifact and checks its OpenMP backend.

## MED named meshes

HDF5-enabled artifacts expose `medMeshNames(path)`, `readMedNamed(path, name, {timeStep: 0, lenient: false})` and `writeMedMulti(path, meshes, names, {version: '4.1.0'})`; paths refer to MEMFS. Named reads include the MED side channel as `mesh.info`. Names must be unique/nonempty, without `/` or `@`; ordinary nodal/element profiles expand with NaN fill. Artifacts without HDF5 throw a capability error. See [MED](./formats/med.md#native-named-meshes-and-profiles).

XDMF series `writePointsCells(mesh)` now stores the mesh's fixed point/cell/side regions once with the shared topology. Existing step, flush and append calls retain them without new arguments; `readMeshSelective` returns the regions at every selected step. See [shared named regions](xdmf_time_series.md#shared-named-regions).

Native serial VTK XML reads (`vtu`, `vtp`, `vts`, `vtr`, `vti`) support multiple pieces without welding, appended raw/base64 arrays, UInt32/UInt64 headers and either byte order, including selective reads and metadata summaries. VTP/VTS/VTR/VTI writers remain inline; optional codecs follow the WASM build. Legacy `vtk` also reads structured points, structured grids and rectilinear grids in ASCII and big-endian binary through the existing `readMesh` API. See [formats](formats.md).

The C++ core also compiles to WebAssembly and ships as an npm package, [`@meshioplusplus/wasm`](https://www.npmjs.com/package/@meshioplusplus/wasm), for reading and writing meshes in the browser or Node.js. (It is one of two "flat" bindings over the same core and shared format-dispatch registry — the other is the [C API](/c_api).)

## Install

```sh
npm install @meshioplusplus/wasm
```

## Lean and full builds

The package ships two pairs of artifacts behind one loader (v17.2.0, roadmap 6.2). The **lean** pair, `meshioplusplus_wasm` and its threaded flavour `meshioplusplus_wasm_mt`, is the default: MIT, every library in it permissively licensed. The **full** pair, `meshioplusplus_wasm_full` and `meshioplusplus_wasm_full_mt`, has the same core and the same JS API, linked with every optional library the Emscripten toolchain cross-compiles: [OpenCASCADE](https://github.com/Open-Cascade-SAS/OCCT) 7.6.3 (reading [STEP](./formats/step.md#reading) and [IGES](./formats/iges.md), and [`fitCad`](#fitcad)), [zstd](./codecs.md) and lz4 (VTK XML block compression) and bzip2 (libMesh `.xda.bz2`/`.xdr.bz2`).

```js
const lean = await loadMeshioPlusPlus();                    // the default, MIT
const full = await loadMeshioPlusPlus({}, { full: true });  // OpenCASCADE + codecs
full.buildFlavour();   // 'full'
full.hasOcct();        // true; hasZstd(), hasLz4(), hasBzip2() too
```

`full` takes the same `variant` (`'auto'`, `'mt'` or `'seq'`) as the lean pair, and only the artifact loaded is fetched: a page that never asks for `full` never downloads it. Five reporters tell which artifact loaded and what it carries, beside `hasCgnslib()`: `buildFlavour()` (`'lean'` or `'full'`), `hasOcct()`, `hasZstd()`, `hasLz4()` and `hasBzip2()`, all false in the lean pair and true in the full one, and `availableFormats()` lists `step` and `iges` among the readers only where OpenCASCADE is linked (`step` is a writer in both: the faceted solid needs no kernel). Reading a STEP or IGES file in the lean artifact, or calling `fitCad`, throws an `Error` that names `loadMeshioPlusPlus({}, {full: true})`. The loader checks `buildFlavour()` against what it was asked for, as it checks `parallelBackend()`, so a `locateFile` that hands the full glue the lean `.wasm` (or the reverse) fails at load with a `MeshioPlusPlusLoadError` (its `full` field set) rather than at the first CAD call.

| Artifact | `.wasm` | `.mjs` | Carries |
| --- | --- | --- | --- |
| `meshioplusplus_wasm` (lean) | 13.5 MB (4.5 MB gzipped) | 94 kB | HDF5, netCDF, CGNS, zlib |
| `meshioplusplus_wasm_mt` (lean, threaded) | 14.2 MB (4.7 MB gzipped) | 113 kB | the same, with OpenMP over pthreads |
| `meshioplusplus_wasm_full` | 26.2 MB (8.8 MB gzipped) | 98 kB | the lean libraries, plus OpenCASCADE, zstd, lz4 and bzip2 |
| `meshioplusplus_wasm_full_mt` (threaded) | 26.6 MB (9.0 MB gzipped) | 118 kB | the same, with OpenMP over pthreads |

Measured for v17.2.0 with Emscripten 6.0.3 (`gzip -9` in parentheses; a server's Brotli does better). OpenCASCADE's modelling, healing, meshing and STEP/IGES toolkits about double the binary, which is why the full pair is opt-in.

**Bundlers.** The lean pair's imports are literal, so Vite, webpack and the rest emit both chunks as before. The full pair's import is computed at run time on purpose: a bundler that saw it would put the full bytes into every lean bundle (the browser viewer, a single-file Vite build, is checked in CI to carry none of them). Under Node and native browser ESM it resolves next to the package's `src/index.mjs` with no setup; a bundled app copies `dist/meshioplusplus_wasm_full*.{mjs,wasm}` into its assets and passes `{ full: true, fullGlue: url }` with a `locateFile` returning the matching `.wasm`.

### Licence

The lean pair is MIT, its bundled libraries (zlib, HDF5, netCDF, CGNS, the vendored gidpost) all permissively licensed. The full pair also links **OpenCASCADE, licensed under the LGPL 2.1 with the OCCT exception**, statically, so it is not plain MIT; the package's licence expression says so (`MIT AND (LGPL-2.1-only WITH OCCT-exception-1.0)`). An application that ships the full artifact takes on the LGPL's obligations for OpenCASCADE (not for meshio++), the central one being that a recipient can relink against a modified OpenCASCADE: the artifact is rebuilt from source by `build/build-wasm-deps.sh --full` (OpenCASCADE at the pinned tag `V7_6_3`, commit `b079fb98…`, unpatched) and `build/configure-wasm.sh --full-only --build`. The package's `THIRD_PARTY_NOTICES.md` lists every bundled component with its version and licence; natively the same kernel stays bring-your-own and dynamically linked ([fit CAD](./fit_cad.md#getting-opencascade)).

### What stays out

Some libraries are in neither pair, each for a reason, and each still answers by name:

- **KaHIP**: no Emscripten port, and a graph partitioner would bloat every bundle; `partition(..., 'kahip')` throws naming `MESHIOPLUSPLUS_WITH_KAHIP`, and `'auto'` resolves to the SFC method.
- **ADIOS2** (`vtx`): no one has made a wasm build of it; `vtx` is absent from `availableFormats()`, and a `.bp` path names ADIOS2.
- **TecIO** (`szplt`): its redistribution terms; absent, and named, the same way.
- **The Kokkos parallel backend**: CMake refuses it under Emscripten; `parallelBackend()` is `seq` or `openmp`.
- **Memory mapping**: there is nothing to map inside MEMFS; the option is accepted and ignored.
- **MMG**: there is no native MMG backend yet ([roadmap §4.3](./roadmap.md#_4-operations)); bringing it to the browser through the full pair is roadmap §6.2.3.

The smoke suite asserts this set exactly on each artifact (`MESHIOPLUSPLUS_WASM_VARIANT=full` and `full-seq` run it on the full pair).

## Publishing a release to npm

Pushing a `v*` tag builds the WASM artifacts and publishes `@meshioplusplus/wasm` at that tag's version through the `publish-npm` job in `.github/workflows/wasm.yml`.

Authentication is npm trusted publishing (OIDC), so there is no long-lived token to expire. It needs a one-time setup by a package maintainer on npmjs.com: on the package's Access page, add a Trusted Publisher for GitHub Actions with owner `loumalouomega`, repository `meshioplusplus`, workflow file name `wasm.yml`, and environment `npm` (matching the job's `environment:`). Until that entry exists, the job falls back to the `NPM_TOKEN` repo secret; a rejected token fails the job's `npm whoami` check with an actionable error instead of npm's cryptic scoped-PUT 404 (the registry masks auth failures on scoped publishes as "not in this registry").

If a publish fails after the artifacts built fine, fix the credential (the trusted-publisher entry, or a fresh `NPM_TOKEN`) and re-run only the failed job with `gh run rerun <run-id> --failed`. The job skips publishing when the version is already on the registry, so re-runs are safe.

## Usage

```js
import { loadMeshioPlusPlus } from "@meshioplusplus/wasm";

const meshio = await loadMeshioPlusPlus();

// Write bytes into the Emscripten virtual filesystem, then read them as a mesh.
const response = await fetch("example.vtu");
meshio.FS.writeFile("/example.vtu", new Uint8Array(await response.arrayBuffer()));
const mesh = meshio.readMesh("/example.vtu");

console.log(mesh.points);          // Float64Array, flat (numPoints * dim)
console.log(mesh.cells[0].type);   // e.g. "triangle"
console.log(mesh.cells[0].data);   // Int32Array connectivity, flat (numCells * nodesPerCell)

// Convert directly (no intermediate JS object), or round-trip through one.
meshio.convert("/example.vtu", "/example.stl");
meshio.writeMesh("/example.msh", mesh, "gmsh");

// What this build can actually read and write, both sorted.
const { readers, writers } = meshio.availableFormats();
```

`regionAdjacency(mesh, regions?)` returns the conforming facets shared by named Cell regions. When the selector list is omitted or empty, the module uses all Cell regions, or multiple cell blocks when fewer than two Cell regions exist. The result carries parent-cell/facet ids and per-facet length or area; see [region adjacency](/region_adjacency).

`findInterface(mesh, regionA, regionB, options?)` returns `{mesh, report}` for conforming or proximity matching; `contactPairs(slaveMesh, slavePoints, masterCells, options?)` returns typed projection arrays; and `splitInterface(mesh, side, options?)` returns `{mesh, report}` after duplicating point fans. Cell, point, facet and subfacet ids are 0-based, matching the JS mesh representation. The `options` shapes are documented in [`index.d.ts`](https://github.com/loumalouomega/meshioplusplus/blob/main/src/wasm/index.d.ts); see also [interfaces and contact](/region_adjacency).

```js
const found = meshio.findInterface(mesh, 'part_a', 'part_b', { mode: 'proximity' });
const pairs = meshio.contactPairs(mesh, 'slave_nodes', 'master_cells', { tolerance: 0.01 });
const split = meshio.splitInterface(mesh, 'crack_sides', { addCohesive: true });
```

### Rendering a mesh

`convertSurface` produces something a surface renderer can draw: a mesh with 3D cells becomes its boundary, anything else passes through, and the result is linearized (a triangle renderer has no mid-side nodes, so `triangle6` connectivity drawn verbatim is visible garbage). Boundary facets inherit their owning cell's data, so a per-cell material or tag still colours correctly.

```js
meshio.convertSurface("/solid.msh", "/solid.vtp");
const vtp = meshio.FS.readFile("/solid.vtp");   // hand this to vtk.js
```

Prefer it over `readMesh` → `extractSkin` → `writeMesh` for anything headed to a renderer: it never materializes a JS mesh, so nothing is copied through the boundary at all, and multi-component arrays keep their exact dtype rather than being widened to `Float64` (they survive the object boundary too since v9.9.0, via `*_components` — see "The mesh object shape" — but at the cost of a copy).

### Applying operations first

`convertSurfaceOps` takes a pipeline and applies it before extracting the surface — all inside C++:

```js
const report = meshio.convertSurfaceOps('/part.msh', '/part.vtp', [
  { op: 'clean', weld: true },
  { op: 'smooth', method: 'taubin', iterations: 20 },
  { op: 'quality' },
]);
console.log(report.steps, report.warnings);
```

Prefer it over chaining the individual operation bindings (`clean`, `smooth`, …): each of those takes and returns a JS `Mesh`, so a pipeline built from them copies the whole mesh across the boundary once per step (and, before v9.9.0, flattened every multi-component array on the first one). An **empty** pipeline is byte-identical to `convertSurface`, which is what lets a viewer use one code path for the plain and the post-operation display — and makes undo a replay of a shortened pipeline rather than a set of inverse operations.

`convertSurfaceOps` takes `{ inFormat, outFormat, keepProvenance, compressVtp }` as its options. A `.vtp` output is zlib-compressed by default; `compressVtp: false` writes it as uncompressed base64 (still binary, never ASCII). The option is ignored for every other output format. Leave it at the default for anything that renders with vtk.js: the uncompressed file carries a 4-byte header (a `header_type="UInt64"` only appears when an array could pass 4 GiB), and vtk.js reads an uncompressed array as `new Float64Array(buffer, headerBytes)`, which throws a `RangeError` at offset 4.

### Handing the surface to a renderer as arrays

`surfaceBuffersOps(inPath, ops, { inFormat, keepProvenance })` runs the same read, pipeline and surface extraction as `convertSurfaceOps` and returns the surface instead of writing it: `{ surface, steps, warnings }`, where `surface` is `{ points, verts, lines, polys, pointData, cellData }`, the arrays a vtk.js `vtkPolyData` is built from. `points` holds x, y and z per point (2-D meshes are padded with z = 0), as a `Float32Array` for float32 points and a `Float64Array` otherwise. `verts`, `lines` and `polys` are `Uint32Array`s in VTK's legacy cell layout (a count, then the ids), filled in the VTP writer's order: vertex blocks, then line blocks, then triangle, quad and polygon blocks, each in block order. `pointData` and `cellData` are lists of `{ name, values, components }`, sorted by name, with cell values in that same cell order. The cell ids are therefore the VTP file's, so `surface:parent_cell` and a picked cell id mean the same in both. Values keep their dtype, except that a 64-bit integer array arrives as an `Int32Array` or `Uint32Array` when every value fits and as a `Float64Array` when one does not, never as a `BigInt64Array`, which vtk.js cannot hold (vtk.js's own VTP reader keeps the low 32 bits instead). `components` is the product of an array's trailing dimensions, so a 3×3 tensor has 9. The [browser viewer](./viewer.md) renders through this call; on a 980,000-triangle surface it takes 0.41 s where writing, reading out and parsing the VTP took 4.08 s ([benchmarks](./benchmarks.md#a-typed-array-channel-for-the-browser-viewer)).

```js
const { surface } = meshio.surfaceBuffersOps('/part.msh', [{ op: 'quality' }], { keepProvenance: true });
const polydata = vtkPolyData.newInstance();
polydata.getPoints().setData(surface.points, 3);
polydata.setPolys(vtkCellArray.newInstance({ values: surface.polys }));
for (const a of surface.cellData)
  polydata.getCellData().addArray(
    vtkDataArray.newInstance({ name: a.name, values: a.values, numberOfComponents: a.components }));
```

This is exactly what the [browser viewer](./viewer.md) does. It is built on this package and is worth reading as a worked example of the whole pipeline — worker, transferable buffers, and vtk.js — as well as being a live client-side format converter you can try at **<https://loumalouomega.github.io/meshioplusplus/viewer/>**.

### Running a settings pipeline

`runPipeline` (v9.11.0) runs a whole [settings document](./pipeline.md) — read → operation chain → write, against MEMFS paths — with **no** surface-extraction tail: what the pipeline produces is what is written.

```js
const report = meshio.runPipeline({
  Input: { Path: '/part.msh' },
  Operations: [
    { Op: 'ConvertCells', Mode: 'simplexify' },
    { Op: 'Gradient', Array: 'temperature' },
    { Op: 'Quality' },
  ],
  Output: { Path: '/part.vtu' },
});
```

It accepts the parsed object, the JSON text, or a MEMFS path ending in `.json` (the wrapper `JSON.parse`s the string forms — the wasm binary carries no JSON parser of its own, so this surface never needed the nlohmann submodule). Note the vocabulary is the settings.json **PascalCase** one (`Op`/`RemoveOrphans`), while `convertSurfaceOps` keeps its pre-existing camelCase op specs — the two dispatch through the same core engine (`apply_pipeline_step`), differing only in spelling, so they cannot drift. The returned report is `{steps: [{op, ...counters}], warnings: []}` with PascalCase counter keys; parsing is strict (unknown ops/keys error by name), and the multi-mesh ops (`Merge`, `Split`, ...) are rejected pointing at the CLI verbs.

## The mesh object shape

Unlike the Python bindings (which hand numpy a zero-copy view straight into the C++ buffer), WASM linear memory and the JS heap are different address spaces, so every value crossing the boundary is copied once. `readMesh` returns, and `writeMesh` accepts, a plain object:

```ts
{
  points: Float64Array,       // flat, row-major: numPoints * dim
  dim: number,                 // 2 or 3
  cells: [
    { type: string, data: Int32Array, nodesPerCell: number }
    // one entry per cell block, data flat row-major: numCells * nodesPerCell
  ],
  point_data?: { [name: string]: DataArray },
  cell_data?: { [name: string]: DataArray[] },   // one array per cell block
  field_data?: { [name: string]: DataArray },

  // Per-array component count, for any array that is not a scalar (v9.9.0).
  point_data_components?: { [name: string]: number },
  cell_data_components?:  { [name: string]: number },   // per ARRAY, not per block
  field_data_components?: { [name: string]: number },
}
```

This deliberately mirrors the Python `Mesh`'s structure (points, a list of cell blocks, `cell_data` as one array per block). Cell connectivity is always `Int32Array` — the C++ core's connectivity dtype is Int64, but node/point counts for any mesh a browser can reasonably hold fit comfortably in 32 bits, and `Int32Array` is far more ergonomic in JS than `BigInt64Array`; a connectivity value that does not fit `Int32Array` throws a catchable `Error` naming it, instead of silently wrapping.

### Data array dtypes

`points` is always `Float64Array` and connectivity is always `Int32Array` (see above), but `point_data`/`cell_data`/`field_data` arrays carry their own dtype as `DataArray = Float32Array | Float64Array | Int8Array | Int16Array | Int32Array | BigInt64Array | Uint8Array | Uint16Array | Uint32Array | BigUint64Array`, instead of every array always widening to `Float64Array` (as it did before v11.2.0). `writeMesh`/`convert`-side callers may hand over any of those classes (or a plain `Array`, treated as `Float64`); every block of one named `cell_data` array must share the same class, or the call throws naming the array.

What actually comes back out of `readMesh` is bounded by the mesh's own in-memory representation, which canonicalizes *within kind* — every floating-point dtype to `Float64Array`, every integer dtype to `BigInt64Array` — rather than preserving the exact input class. Concretely: a `Uint8Array` or `Int32Array` array you write is accepted, but a mesh you then read back has it as `BigInt64Array`, with every value exact (not the double it silently became before v11.2.0). `BigInt64Array` elements are JS `bigint`, not `number`; `Number(x)` or `Array.from(arr, Number)` converts them where you need ordinary numbers (see `XdmfTimeSeriesWriter.writeDataArrays`, which does this for you on the values it accepts).

### Multi-component (vector / tensor) data arrays

Data arrays are flat, and **a flat typed array carries no shape** — the same problem `xdmfSeriesWriteDataArrays`' own `components` argument solves. Each of the three data maps therefore has a sibling `*_components` object giving the per-entity width of any array that is not a scalar:

```js
const mesh = {
  points, dim: 3, cells,
  point_data: { velocity: new Float64Array(n * 3) },   // interleaved
  point_data_components: { velocity: 3 },
};
meshio.writeMesh('/out.med', mesh, 'med');
const back = meshio.readMesh('/out.med', 'med');
back.point_data_components.velocity;   // 3
```

- **An absent name means one component.** A caller that never sets `*_components` gets exactly the pre-v9.9.0 scalar behaviour, and `readMesh` only writes an entry for genuinely multi-component arrays — so a scalar-only mesh comes back with three empty objects and no existing consumer sees a new key.
- **`cell_data_components` is one value per array, not per block**: every block of a named `cell_data` array must agree on its component count, which is the uniform mesh API's own invariant.
- A length that is not a multiple of the declared count, or a non-positive/non-integer count, is a catchable `Error` naming the array.

Before v9.9.0 there was no such metadata: an `(n,3)` array re-entered C++ as `(3n,1)`, so writing a vector field to MED produced a file the reader itself rejected (`"field data size does not match its declared shape"`), and every object-based operation silently passed the array through instead of gathering it, leaving stale values of the wrong length. The path-based `convert`/`convertSurface`/`convertSurfaceOps` calls were never affected, because they never materialize a JS mesh at all.

**Ragged cell blocks** (polygon/polyhedron blocks with a varying node count per cell, e.g. MED Voronoi polygons or OpenFOAM general polyhedra) cross the boundary as flat CSR arrays instead of a single rectangular `data`/`nodesPerCell` pair, since embind has no efficient representation for a nested array of arrays:

- **Polygon** (1-level ragged, jagged rows): `{type, data, rowOffsets}` — `data` is every row's node ids concatenated, `rowOffsets` is each cell's start index into `data` (length `numCells + 1`).
- **Polyhedron** (2-level ragged, cell → faces → node ids): `{type, data, faceOffsets, cellOffsets}` — `data` is every face's node ids concatenated, `faceOffsets` is each face's start index into `data` (length `totalFaces + 1`), `cellOffsets` is each cell's start index into the face list (length `numCells + 1`).

```js
// A triangle and a 4-gon in one polygon block.
const poly = {
  points: new Float64Array([0, 0, 0, 1, 0, 0, 1, 1, 0, 2, 0, 0, 2, 1, 0]),
  dim: 3,
  cells: [{ type: 'polygon', data: new Int32Array([0, 1, 2, 1, 3, 4, 2]), rowOffsets: new Int32Array([0, 3, 7]) }],
};
m.writeMesh('/ragged.med', poly, 'med'); // MED is the ragged-polygon-capable writer (POG/POG2)
```

These are the same three arrays the [C API](/c_api) hands out through its `mio_poly_conn` snapshot and the Fortran/Julia/R bindings expose in each language's own shape — see [Polyhedra and ragged cells](/polyhedra) for the shared vocabulary and the winding rule. Malformed offsets (non-monotonic, or running past the end of `data`) throw naming the offending array rather than reading out of range.

`readMesh` reports whichever shape the source held; `writeMesh` accepts either, and the target format's own writer decides what it can represent. Polyhedron blocks are written by `vtu` (`VTK_POLYHEDRON` face streams), `ensight` (`nfaced`), `cgns` (`NFACE_n`), `med` (`MED_POLYHEDRON`, since v9.19.0) and `openfoam`; a format that cannot hold one (legacy `vtk`, `vtp`, …) throws naming the format rather than silently dropping data. The JS boundary itself carries polyhedron blocks unchanged through operations like `clean`/`merge`/`convertCells`.

### Property sets

`Begin Properties <id>` blocks (Kratos material data; currently only MDPA populates them) ride on the mesh object itself, like regions — `mesh.propertySets` is an array of `{id, values}`, where each `values` entry is `{key, values, text, isTable, components?}`. Exactly one of `values` (a `Float64Array`, a table's `components` giving its row width when `> 1`) and `text` carries the value: a plain number is a one-element `values`; an inline `Begin Table` is `values` shaped `(n, k)` with `isTable` set and `key` holding the table header's arguments verbatim; anything else (a constitutive-law name, a bracketed vector/matrix) is kept verbatim in `text`, which is what makes an unrecognized value lossless.

```js
const withProps = {
  ...mesh,
  propertySets: [
    { id: 1, values: [{ key: 'YOUNG_MODULUS', values: new Float64Array([2.1e11]), text: '', isTable: false }] },
  ],
};
meshio.writeMesh('/out.mdpa', withProps, 'mdpa');
meshio.readMesh('/out.mdpa', 'mdpa').propertySets;   // the same, read back
```

Property sets are keyed by **id, not entity index**, so shape-preserving operations (`clean`, `smooth`, `transform`, `attachQuality`, the `data*` ops) carry them through unchanged while restructuring/multi-input ones (`merge`, `cropBbox`/`cropPlane`/`cropPredicate`, `split`, `partition`, `diff`) do not — the same rule the C++ core documents in `mesh_api.hpp`.

### Side channel (info)

Some formats carry metadata a generic `Mesh` cannot represent at all: OpenFOAM patch names/types, MED field units and family groups, MDPA entity names, Ansys (and UNV, whose groups also arrive as `mesh.regions` since v15.6.0) named point and cell sets, Gmsh `$Entities` bounding-entity tags, Exodus info/QA records. `readMeshSelective(path, {info: true})` attaches it as the result's `.info` for the formats that have one; `writeMesh(path, mesh, format, {info})` writes it back (or, with no explicit `options.info`, reuses `mesh.info` when its own `format` field matches the write target — a `{info: true}` read round-trips through a write with no extra plumbing):

```js
const mesh = meshio.readMeshSelective('/case.foam', { format: 'openfoam', info: true });
mesh.info;   // { format: 'openfoam', patches: [{ familyId, names, type? }, ...] }
meshio.writeMesh('/out.foam', mesh, 'openfoam', { info: mesh.info });   // or just `{ info: mesh.info }` omitted -- mesh.info.format already matches
```

Each format's `info` shape (see `doc/wasm.md`'s TypeScript-adjacent field names below match `index.d.ts`'s `MeshInfo` union exactly):

| Format | Shape | Writable |
| --- | --- | --- |
| `openfoam` | `{patches: [{familyId, names, type?}]}` | yes |
| `med` | `{pointTags, cellTags, meshName, description, unitTime, unitCoords, pointTagGroups, cellTagGroups, skippedConstructs, fieldUnits, stepMeta, fieldTimeValues}` (`skippedConstructs`/`fieldUnits`/`stepMeta` are read-side, lenient-mode diagnostics only; `fieldTimeValues` is always filled) | yes |
| `mdpa` | `{entityNames: [{name, isCondition}], skippedConstructs}` — properties are `mesh.propertySets`, not here | yes Since v16.27.0 it also carries `modelPartData`, `tables`, `geometries`, `meshBlocks`, `subModelParts` and `rawBlocks` (see [MDPA](formats/mdpa.md#the-blocks-the-mesh-cannot-hold-v16-26-0)). |
| `ansysinp` / `unv` | `{pointSets: {name: number[]}, cellSets: {name: number[][]}}` (identical shape on both) | yes |
| `gmsh` | `{boundingEntities: number[][], periodic: {entity: number[], affine: number[], nodePairs: BigInt64Array}[]}` | yes (4.1 or `gmsh22`) |
| `exodus` | `{infoRecords: string[]}` | **no** (read-only; there is no Info-bearing Exodus writer) |
| `mfem` | `{patches: [{degrees, knots, spaceDim, controlPoints, weights, fields}], fieldNames, fieldComponents, topology: {dimension, numVertices, patchVertices, patchAttributes, boundaryVertices, boundaryAttributes, edges} \| null, cellPatch, cellSpan, boundaryCellPatch}` — the patches of an `MFEM NURBS mesh` (v16.40.0), flat arrays with their widths beside them, 0-based ([MFEM](formats/mfem.md#the-patches-themselves-the-nurbs-side-channel)) | yes (as an `MFEM NURBS mesh`, while the mesh still is their sampling) |
| `cad_json` (and, in the full artifact, `step` and `iges`, read only) | `{patches: [{degrees, knots, controlPoints, weights}], breps, faces: [{id, name, patch, swappedNormal, trimmed, rational, reducedKnots, controlPointIds, loops: [{type, trims, extra}], embeddedEdges, surfaceExtra, extra}], edges: [{id, name, curve, topology, extra}], vertices, extra}` — the trimmed NURBS B-rep of a [Kratos CAD JSON](formats/cad_json.md) file (`CadJsonInfo`, with `BrepCurve`, `BrepTrim` and `BrepRef`); ids are -1 where the file names the entity instead, every index is 0-based and the `extra` strings hold the JSON members meshio++ does not model. In the [full artifact](#lean-and-full-builds) a `step` or `iges` file read with `{info: true}` returns its B-rep in this shape too (still `format: 'cad_json'`); the lean artifact has no reader for either | yes, as CAD JSON (`step` and `iges` info is read-only) |

Every shape includes a `format` discriminator matching the table's left column. `pointsOnly`/`arrays`/`timeStep` (`readMeshSelective`'s other selective-read options) reach the read only for the formats whose info-bearing reader takes them natively (`med`/`mdpa`/`gmsh`/`exodus`/`unv`/`mfem`); `openfoam`/`ansysinp` have no selective-read path with or without `info`, exactly as they have none without it. `info: true` for a format with no side channel is silently ignored (no `.info` on the result) rather than thrown, so a caller can always pass it and check `mesh.info` itself; writing `info` for a format with no Info-bearing writer (`exodus`, or any other format) throws naming it. There is no generic `registry_read_info`/`registry_write_info` hook in the C++ core reaching every flat binding — this is a WASM-only dispatch table calling each format's own Info-bearing reader/writer directly, the same per-format special-casing the Python bindings already do; a shared core hook is a roadmap remainder.

MDPA's `info.subModelParts` entries carry `name`, `data`, `tables`, `geometryIds` and `constraintIds`; the last two retain raw file ids, ordering and duplicates and are not remapped by operations. `writeMesh(path, mesh, 'mdpa', {info})` restores nested membership-only parts. Constraints themselves remain opaque `rawBlocks` content. Omitting the two new fields on write is equivalent to empty lists.

## Mesh operations

Region-backed sets↔integer data conversions use `setsToData(mesh, "point", {dataName, joinChar: "-", order: []})` and `dataToSets(mesh, "cell", key)`, returning a new mesh without changing geometry. Both accept only point/cell locations. Omit `order` for native region-name order; otherwise it must list every distinct set name once. Later overlaps win, uncovered rows are -1, and converted regions/source fields are removed; side regions survive. See [sets/data semantics](data_manage.md#sets--integer-data).

Besides read/write/convert, the module exposes the mesh operations on the same mesh-object shape: `extractSurface(mesh [, recordParentIds])`, `extractSkin(mesh [, linearize])`, `attachQuality(mesh)`, `sniffFormat(path)`, `computeBandwidth(mesh)`, and `reorder(mesh, method)` (method `"rcm"`, `"morton"`, or `"hilbert"`). `reorder` returns `{ mesh, nodePermutation, cellPermutations }` — the renumbered mesh plus the applied old→new permutations (`Int32Array`s). See [Reordering / renumbering](./reorder.md). Comparison is exposed as `meshesEqual(meshA, meshB, atol, rtol, unordered)` (boolean) and `diff(meshA, meshB, atol, rtol, unordered)` (a report object with `verdict`, `points`, `blocks`, `pointData`/`cellData`/`fieldData`); named point/cell sets are not compared. See [Mesh comparison (diff)](./diff.md). Merging is exposed as `merge([meshA, meshB, ...], weld, atol, sourceTag, dataPolicy, dropDuplicateCells)` (`dataPolicy` is `"intersection"` or `"fill"`), returning a new mesh object; named point/cell sets are not carried. See [Merge / combine](./merge.md). The editing/stats bundle is exposed too: `transform(mesh, matrix, rotateVectorData)` (a 16-element row-major matrix) → mesh; `clean(mesh, weld, atol, removeOrphans, dropDegenerate, dropDuplicateCells)` → `{ mesh, pointsWelded, ... }`; `cropBbox(mesh, lo, hi, mode, recordIds)` / `cropPlane(mesh, point, normal, mode, recordIds)` (`mode` `"all"`/`"any"`) / `cropPredicate(mesh, array, compare, value, recordIds)` (a scalar `cell_data` comparison, `"<"`/`"<="`/`">"`/`">="`/`"=="`/`"!="`; **no `mode`** — a per-cell value has nothing for an all/any rule to reduce) → mesh; `split(mesh, by, tagName)` → an array of `{ key, mesh }`; and `stats(mesh)` → an object of geometric measures (`bboxMin`/`bboxMax`/`extent`/`centroid`, `cellTypeCounts`, `totalArea`, `signedVolume`, `unsignedVolume`, `numInverted`). Element-representation conversion is exposed as `convertCells(mesh, mode, recordParentIds)` with `mode` `"linearize"`, `"simplexify"`, or `"elevate"`, returning a new mesh; a polyhedron block under `"simplexify"` and the full-Lagrange targets (`quad9`/`hexahedron27`) under `"elevate"` throw a catchable `Error`. **Polyhedral refinement** is exposed as `subdivide(mesh, recordParentIds)`: one polyhedral child per face of every eligible 3D cell, connected to a new interior point, with no per-type template table — tabulated types (reduced to corners for a quadratic variant) and existing polyhedron blocks are handled uniformly, and the result is automatically conforming, unlike `refine`. There is no point map, unlike `convertCells` (subdivide never prunes or renumbers a point); a cell whose faces are not a closed orientable surface throws a catchable `Error`. See [subdivide](/subdivide). **Polyhedral coarsening**, the many-to-one counterpart, is exposed as `agglomerate(mesh, targetGroupSize, returnMaps, options)` (`options`: `{mergeCoplanarFaces, coplanarAngle, minSphericity}`; with `returnMaps` the result adds `numFacesMerged` and `numRejected`): greedy seed-and-grow over the mesh's shared-face dual, absorbing face-adjacent neighbours by accumulated shared-face area until each group reaches `targetGroupSize` (default 8) members, then emitting one polyhedron per group whose faces are exactly its external boundary — conserving volume exactly, since internal faces are simply dropped rather than re-triangulated. Like `subdivide` there is no point map (points are never pruned or renumbered — `clean(mesh, ..., removeOrphans: true)` is the follow-up for a minimal point set), and a non-manifold input (a face shared by three or more cells) throws a catchable `Error` naming the face. It is also a `convertSurfaceOps`/`runPipeline` pipeline step (`{op: 'agglomerate', targetGroupSize}`), reached through the same generic `pipe_op_table()` dispatch every other step goes through with no extra WASM code. See [agglomerate](/agglomerate). Uniform refinement is exposed as `refine(mesh, levels, recordParentIds, options)`, subdividing every cell into same-type children (`triangle`/`quad` into 4, `tetra`/`wedge`/`hexahedron` into 8) with shared mid-entity nodes, so the result has no hanging nodes; a higher-order cell, a `pyramid`, or a ragged block throws a catchable `Error`. The optional fourth argument selects a **subset** to refine — `{cells, region, array, compare, value, closure, recordLevels, recordHierarchy}`, at most one selector — in which case the hanging nodes that leaves are resolved by the closure and, for `'redgreen'` and `'propagate'`, the output is still conforming (`'balanced'` deliberately keeps them and reports each in `refine:hanging`); the `convertSurfaceOps` pipeline op `{op: 'refine', ...}` takes the same fields, where the comparison is spelled `compare` because `op` is the step's own discriminant. `recordHierarchy` attaches the persistent `refine:cell_id`/`refine:parent_id` cell_data — a link between the meshes a multigrid caller keeps across passes, not a tree inside one — and forces `refine:entity` (the prolongation stencil) even when the closure leaves no hanging node; see [refine](/refine#refinecell_id-and-refineparent_id). **Green-element undo**, the two-mesh op restoring `refine`'s transitional cells to their coarse parent, is exposed as `undoGreen(coarse, fine)` → `{ mesh, numGroupsUndone, numCellsRemoved }` — a lookup and substitution read verbatim from `coarse` (no template inversion, no winding repair), needing `fine` to carry `refine:cell_id`/`refine:parent_id`/`refine:level` (i.e. it must come from a `refine(coarse, ..., {recordHierarchy: true, recordLevels: true})` call). The six reserved `refine:*` arrays are dropped from the output. Being a two-mesh op it is deliberately **not** reachable as a `convertSurfaceOps` pipeline step — `{op: 'undoGreen'}` throws the same "unknown operation" catchable `Error` `Merge`/`Interpolate`/`Split`/`Diff` already throw there, since none of the two-mesh ops ever reach `pipeline_op_table()`. See [green-element undo](/undo_green). **Mass-preserving field transfer**, the sibling of `interpolate` that conserves `sum(value * measure)` over the region two meshes share rather than sampling pointwise, is exposed as `conservativeInterpolate(source, target, arrays, defaultValue, onConflict)` → mesh (always Float64 arrays). Both meshes are simplexified first, accepting ragged/polyhedron blocks for free; unlike `interpolate`, an empty `arrays` means every source point_data **and** cell_data array. Like `interpolate` it is a two-mesh op and so is **not** reachable as a `convertSurfaceOps` pipeline step (`{op: 'conservativeInterpolate'}` throws the same "unknown operation" catchable `Error`). See [conservative interpolation](/conservative_interpolate). Partitioning is exposed as `partition(mesh, nparts, method, imbalance, mode, seed, recordIds, ghostLayers, weightsKey)` → an array of `{ partId, mesh }` (exactly `nparts` entries, blocks kept 1:1 with the input, unlike `split`) and `partitionLabels(mesh, nparts, method, imbalance, mode, seed, weightsKey)` → one label array per cell block. **Only the SFC method exists in the WASM build** — KaHIP is never compiled in (no Emscripten port, and it would bloat the bundle), so `method: "kahip"` throws a catchable `Error` naming `MESHIOPLUSPLUS_WITH_KAHIP` and `"auto"` always resolves to SFC — a deliberate constraint, recorded in the roadmap's [Non-goals](/roadmap#non-goals-and-decisions-taken). Smoothing is exposed as `smooth(mesh, method, iterations, lambda, mu, fixBoundary, preserveFeatures, featureAngle, guardInversion, frozen)` with `method` `"taubin"` (the default, shrink-free), `"laplacian"` (stronger per pass, but shrinking) or `"odt"` (tet-only: each free interior vertex moves toward the volume-weighted average of its incident tets' circumcenters), returning `{ mesh, numNodesMoved, maxDisplacement, numSkippedInversion }`; only the point coordinates move, so connectivity and every data value come through unchanged. Boundary nodes, feature nodes and the nodes of blocks with unknown edge topology are pinned by default; `frozen` is an optional array of 0-based point ids to pin outright, unioned with those pins, and an out-of-range id throws a catchable `Error` naming it. Volume decimation is exposed as `decimateVolume(mesh, ratio, targetCells, maxError, placement, preserveBoundary, preserveFeatures, featureAngle, frozen)` (`decimate`'s tetrahedral-mesh sibling: quadric-error-metric tet-edge collapse, exactly one of `ratio`/`targetCells`/`maxError` non-negative) → `{ mesh, tetsRemoved, pointsRemoved, collapsesRejected, maxErrorApplied }`; `preserveBoundary` defaults to `false` here (unlike `decimate`'s `true`), since a volume mesh's outer surface is usually interior geometry a solver still wants simplified rather than a boundary to protect. See [decimate](/decimate) and [decimate_volume](/decimate_volume). The two cutters are exposed as `slice(mesh, origin, normal, recordParentIds)` — the planar cross-section, one dimension below the cut cells — and `isosurface(mesh, array, isovalues, component, recordParentIds)`, its data-driven sibling: the level set of a scalar `point_data` array. `isovalues` accepts a number or an array (several contours land in one mesh, tagged per cell with `iso:value` and `iso:index`), `component` is negative for the row magnitude, and naming a `cell_data` array throws a catchable `Error` — cell data is piecewise constant and has no level set, so convert it with `dataCellToPoint` first.

`computeNormals(mesh [, pointNormals, cellNormals, weight, splitAngle, recordParentIds, region])` returns `{ mesh, numIsolated, numUndefined, numDegenerate, numSplitPoints, numAddedPoints, quality }`: point (and optionally cell) normals of a surface written as the `normals` array, with vertices split at creases when `splitAngle` is a number of degrees in `[0, 180]` (a negative `splitAngle`, the default, means no split — the copies a split makes are appended after the original points and cells keep their numbering). `weight` is `"angle"` or `"area"`. It never reorients: check `quality.inconsistentPairs`. See [normals](/normals).

`checkQuality(mesh [, require, maxInverted, maxDegenerate])` returns `{ passed, numCells, numInverted, numDegenerate, checks, summary }`: the [quality gate](/mesh_quality#quality-gate), every check with its violations and worst cell.

`featureEdges(mesh [, featureAngle, feature, boundary, nonManifold, inconsistent, region])` returns `{ mesh, numFeature, numBoundary, numNonManifold, numInconsistent }`: the sharp, open, non-manifold and inconsistently wound edges as a `line` mesh with `feature:kind`/`feature:angle` (see [feature edges](/feature_edges)). `hausdorffDistance(a, b [, faceSamples, regionA, regionB])` returns `{ distance, aToB, bToA, meanAToB, rmsAToB, meanBToA, rmsBToA, numSamplesA, numSamplesB, worstPointA, worstPointB }` (see [Hausdorff distance](/hausdorff)). `editRegions(mesh, edits)` applies `{op, inputs, output, dim, tag, keepInputs}` edits in order and returns the edited mesh (see [editing regions](/regions#editing-regions)). `matchPeriodicNodes(mesh, slave, master, matrix [, atol, requireComplete])` returns `{ slave, master, unmatched, numFixed, maxResidual }` with 0-based `Int32Array` ids, `matrix` being 16 numbers, row-major (see [periodic node pairs](/periodic)).

`segment(mesh [, featureAngle, normalTolerance, curvatureTolerance, minPatchCells, region, prefix])` returns `{ mesh, numPatches, numMerged, numFeatureEdges, numBoundaryEdges, numNonManifoldEdges }`: the input with `segment:patch` cell data, `segment:adjacency` field data and one cell region per patch ([segment](/segment)). `fitPrimitives(mesh [, patchData, tolerance, relative, minInlierFraction, types, maxIterations])` returns `{ mesh, numPatches, numPlane, numSphere, numCylinder, numCone, numTorus, numUnfitted }`, with `types` an array of `'plane'`, `'sphere'`, `'cylinder'`, `'cone'`, `'torus'` ([fit primitives](/fit_primitives)). <a id="fitcad"></a>`fitCad(mesh [, options])` turns those fits into an analytic CAD solid ([fit CAD](/fit_cad)) in the [full artifact](#lean-and-full-builds), which links OpenCASCADE; the lean one throws naming `{full: true}`. `options` is `{stepPath, tolerance, exactEdges, types, meshDeflection, maxDemotionRounds, freeform, freeformDegree, freeformMaxControl, freeformSplit, fillets, threads}` (an unknown key is an error), and the result carries `mesh` (the solid re-meshed, with the `fit_cad:*` cell data), `valid`, the volumes, `maxDeviation`, every face count of the native report in camelCase (`numPlane`, `numFillet`, `numHelical`, `numChamfer`, `numVertexBlend`, `numConicBlend`, `numTaperedThread`, ...) and `info`, the solid's B-rep in CAD JSON's shape, which `writeMesh(path, mesh, 'cad_json', {info})` writes and the `iga*` operations take. With `stepPath` the solid is written to MEMFS as AP214 STEP.

`igaRefine(info [, options])`, `igaSample(info [, options])` and `igaEvaluate(info, patch, params [, options])` act on the NURBS model in `readMeshSelective(path, {info: true}).info` (`format: 'mfem'` or `'cad_json'`): `igaRefine` takes `{subdivisions, elevate}` (a number, or one per direction) and the sampling options `{sampleSubdivisions, trimSamples, edges, trim}` and returns `{ mesh, info }`, the refined model sampled and the refined model for `writeMesh(path, mesh, format, {info})`; `igaSample` returns the same for the model as it is; `igaEvaluate` takes a flat `params` array (one parametric coordinate per direction of the patch, per point), `{derivatives, trimSamples}` and returns `{ spaceDim, points, fields, fieldComponents, derivatives?, normals?, inside? }`. See [isogeometric layer](/iga).

`tessellate(mesh [, levels, curved, fields, recordStencil])` (defaults `2`, `true`, `true`, `false`) returns `{ mesh, numSourcePoints, numSourceCells, numCurvedSourceCells, numPassThroughSourceCells, numPoints, numCells, numBoundaryFacets, numInteriorFacets, numFacetsWithBadCount }`: every `triangle6`, `quad8`, `quad9`, `tetra10` and `hexahedron27` cell sampled on a reference lattice through its own quadratic shape functions, curved surface cells becoming one `triangle` block and curved volume cells one `tetra` block after the unchanged pass-through blocks ([tessellation](/tessellation)). The mesh carries `tessellate:source_point` (Int64 point data, -1 for a new point), `tessellate:source_cell` and `tessellate:sub_index` (Int64 cell data); `fields` interpolates point data and copies cell data onto the children, and `recordStencil` adds the `(P, 27)` `tessellate:stencil`/`tessellate:weights` point data. The facet counts cover the output `tetra` block, so `numFacetsWithBadCount` is 0 for a conforming tessellation; `levels` must be an integer in 0..20.

`gradient(mesh, array, operator, method, location, output, component, overwrite)` returns `{ mesh, numSkipped, numFallback }` — the gradient, divergence or curl of a `point_data` field. Unlike the cutters it changes no geometry: it attaches one array and hands the mesh back. **Its result is an `(n, 3)` or `(n, 9)` array, so its width travels in the returned mesh's `point_data_components` / `cell_data_components` sibling maps** (see "The mesh object shape" above) — a JS caller that rebuilds a mesh by hand must carry them, or the array re-enters C++ flattened. Note that `component` is negative for **every** component here, deliberately the opposite of `isosurface`'s sentinel. It is also a `convertSurfaceOps` pipeline step (`{op: 'gradient', array, operator, method, location, output}`).

`hessian(mesh, array, method, location, output, overwrite)` returns `{ mesh, numSkipped, numFallback }` — the Hessian (second derivative) of a **scalar** `point_data` field, `gradient`'s companion one order further. A composition of TWO `gradient` calls, not a new numerical kernel: the field is differentiated once (point location, regardless of the caller's own `location`), then that `(n, 3)` gradient is differentiated again with the default gradient operator, producing `(n, 9)` — the flattened row-major 3x3 Hessian. `method` is forwarded to BOTH internal passes. Like `gradient`, the result's `(n, 9)` width travels in the returned mesh's `point_data_components` / `cell_data_components` sibling maps. A field that is at most linear has an exactly zero Hessian everywhere — the one mesh-shape-independent guarantee; a genuinely quadratic field's composition is exact on a structured/symmetric mesh away from its own boundary and approximate on an irregular mesh (see [second derivatives](./hessian.md)). Input must have exactly one component. It is also a `convertSurfaceOps` pipeline step (`{op: 'hessian', array, method, location, output}`).

`estimateError(mesh, array, method, marking, markingValue, output, marked, overwrite)` returns `{ mesh, globalError, numSkipped, numMarked }` — the Zienkiewicz-Zhu recovery-based error indicator of a `point_data` field, a composition of `gradient` with the measure-weighted point↔cell averaging round trip, not a new kernel. Like `gradient` it changes no geometry: `error:zz` (Float64) is always attached, and `error:marked` (Int64 0/1) too when `marking` is not `"none"` — so `refine`'s own `where` selector needs no change at all. Cells that cannot be evaluated read NaN in `error:zz` and `0` (never NaN) in `error:marked`, counted in `numSkipped`. It is also a `convertSurfaceOps` pipeline step (`{op: 'estimateError', array, method, marking, markingValue, output, marked}`). See [error estimation](./error.md), [field derivatives](./gradient.md), [second derivatives](./hessian.md), [transform](./transform.md), [clean](./clean.md), [crop](./crop.md), [split](./split.md), [stats](./stats.md), [cell conversion](./convert_cells.md), [polyhedral refinement](./subdivide.md), [polyhedral coarsening](./agglomerate.md), [refine](./refine.md), [green-element undo](./undo_green.md), [partitioning](./partition.md), [smoothing](./smooth.md), [slicing](./slice.md), and [isosurfaces](./isosurface.md).

::: tip Reachable from `loadMeshioPlusPlus()` since v7.4.0
Before v7.4.0 the geometry operations above were bound in the WASM module but **not forwarded by the package wrapper**, so they were unreachable through `loadMeshioPlusPlus()` (only file I/O and the `data_*` operations were). They are all forwarded now.

### Index maps

`clean`, `decimate`, `decimateVolume`, `convertCells`, `refine`, `cropBbox`, `cropPlane` and `cropPredicate` accept a trailing `returnMaps` argument (default `false`, so the return shape is unchanged from before v11.2.0): when `true`, the result also carries `pointMap` (`Int32Array`, input point index → output point index, `-1` if pruned/collapsed away) and `cellMaps` (one `Int32Array` per **input** cell block, input cell index → output index within the corresponding output block, `-1` if dropped) — `decimate`'s own `-1`/survivor convention is documented on [decimate](/decimate). `subdivide` and `agglomerate` also take `returnMaps`, but their shape differs: `subdivide` never prunes or renumbers a point, so there is no `pointMap` at all, only `cellMaps`; `agglomerate` returns a single **flat** `cellMap` (global input cell index → global output cell index) instead of per-block `cellMaps`, since an output cell's group membership does not depend on which input block it came from. `split` and `partition` attach `pointMap`/`cellMaps` **per piece** when `returnMaps` is set — each piece's maps are relative to the whole input mesh, not to the piece. `merge` attaches `pointMaps`/`cellMaps` **per input mesh** (not per input block) when `returnMaps` is set: input `i`'s local point/cell index → the output index.

Without `returnMaps`, use the `recordIds`/`recordParentIds` flags instead, which attach the same provenance as ordinary data arrays (or `partitionLabels` for the raw per-cell assignment) — `smooth` needs neither, since it never adds, removes or renumbers a node or a cell. `extractSurface`/`extractSkin` have no index map in the C++ core at all (tracked as a roadmap remainder — see [roadmap](/roadmap)); `undoGreen`'s per-block cell maps are likewise not exposed here, a documented flat-binding gap shared with its C API/Fortran counterparts.
:::

The [data operations](./data_operations.md) — which act on `point_data` / `cell_data` / `field_data` and never modify the geometry — are exposed as `dataDrop(mesh, location, names, ignoreMissing)`, `dataKeep(...)`, `dataRename(mesh, location, from, to)`, `dataPointToCell(mesh, names, suffix)`, `dataCellToPoint(mesh, names, weight, suffix)`, `dataCalc(mesh, expression, location, outputName, overwrite)`, `dataCondition(mesh, location, names, mode, lo, hi, scope, nanPolicy, nanReplacement, suffix)`, `tensorInvariants(mesh, location, names, outputs, prefix, suffix, overwrite)` — each returning a new mesh — and `dataInfo(mesh)`, which returns an array of per-array summary objects (`location`, `name`, `dtype`, `shape`, `numBlocks`, `numEntries`, `numComponents`, `numValues`, `min`, `max`, `mean`, `minPerComponent`/`maxPerComponent`/`meanPerComponent`, `numNan`, `numInf`, `numFinite`, `inconsistentBlocks`). Enumerations cross as strings: `location` is `"point"`/`"cell"`/`"field"`, `weight` is `"uniform"`/`"measure"`, `mode` is `"clamp"`/`"normalize"`/`"standardize"`, `scope` is `"component"`/`"magnitude"`, `nanPolicy` is `"ignore"`/`"replace"`/`"fail"`, and `outputs` is a JS array of any of `"mises"`/`"principal"`/`"hydrostatic"`/`"deviatoric"` (empty means all four). A malformed `dataCalc` expression throws a catchable `Error`. See [data operations](./data_operations.md), [array management](./data_manage.md), [averaging](./data_average.md), [expressions](./data_calc.md), [conditioning](./data_condition.md), [tensor invariants](./tensor_invariants.md) and [data summary](./data_info.md).

`dataIntegrate(mesh, arrays)` returns an array of per-array field-integral objects (`name`, `numComponents`, `domain`, `regions`) — a cell-measure-weighted total/mean of one or more `cell_data` arrays, `gradient`'s integration counterpart. `domain` and each entry in `regions` share one shape (`numCells`, `numSkipped`, and parallel `totalPerComponent`/`meanPerComponent`/`domainMeasurePerComponent`/`numNanPerComponent` arrays — the `dataInfo` `*PerComponent` convention), `regions` carrying one entry per named Cell region present on the mesh (not a partition: a cell in two regions contributes fully to both). `arrays` empty/undefined means every `cell_data` array. Unlike `gradient`/`estimateError` it is **read-only** — like `dataInfo`, it changes no geometry and is not a `convertSurfaceOps`/pipeline step. A `point_data`-only name throws a catchable `Error` naming `dataPointToCell` as the fix. See [field integration](./field_integration.md).

## Transient (time-series) XDMF

`createXdmfTimeSeriesWriter(path, { dataFormat, gzipLevel })` is the **one stateful** thing in this API: every other binding is a pure function over a mesh object, but a time series writes the mesh **once** and then appends one cheap step per solve, and its `.xdmf` light data can only be written when the collection is complete. See [XDMF time series](./xdmf_time_series.md).

```javascript
const w = m.createXdmfTimeSeriesWriter('/series.xdmf');   // 'HDF' by default
w.writePointsCells(mesh);                                 // the static grid, once
for (let k = 0; k < nsteps; ++k) {
    w.writeData(k * dt, stepMesh(k));   // only point_data/cell_data are used
}
w.close();                              // the files appear HERE

const xdmf = m.FS.readFile('/series.xdmf');
const h5 = m.FS.readFile('/series.h5'); // 'HDF' writes TWO files — see below
```

The returned object has `writePointsCells(mesh)`, `writeData(time, mesh)`, `finalize()`, `numSteps()`, `finalized()` and `close()`. `close()` finalizes if needed and releases the handle; it is safe to call twice and safe to call from a `finally`. After it, every method throws a catchable `Error` — the underlying handle is an index into a module-local table, not a pointer, so a stale handle can never be a use-after-free.

::: warning `'HDF'` writes TWO files into the virtual filesystem
The `.xdmf` at the path you gave **and** its sibling heavy-data file `<path minus extension>.h5` (a sibling of the `.xdmf`, not a file in the current directory — this differs from the pure-Python `TimeSeriesWriter`). Copy both out of `FS`; an `.xdmf` without its `.h5` is unreadable. `'XML'` writes one self-contained file; `'Binary'` writes the `.xdmf` plus one `<path minus extension><n>.bin` per array. Nothing is on the filesystem until `finalize()`/`close()` runs.
:::

`'HDF'` works here: the shipped artifact links a wasm32 HDF5 (v8.0.0+). Read a series back with the ordinary `readMesh` (which resolves the temporal collection structurally and gives you the first step), `readMeshSelective(path, { timeStep: k })` for a particular one (negative counts from the end), or `readMetadata(path).timeValues` to see the steps without loading any payload.

::: tip Why a handle and not a class
The raw embind surface is an opaque integer handle plus seven free functions (`xdmfSeriesCreate`, `xdmfSeriesWritePointsCells`, `xdmfSeriesWriteData`, `xdmfSeriesFinalize`, `xdmfSeriesNumSteps`, `xdmfSeriesFinalized`, `xdmfSeriesFree`), deliberately **not** an embind `class_`, and the object above is the wrapper's ergonomic face of it. `@meshioplusplus/wasm` never hands JS a live C++ object — `NDArray`, `CellBlock` and `Mesh` are all internal — so an embind class instance, with the Emscripten-specific `.delete()` it comes with, would be the API's one exception. Free functions also all go through the same C++ error wrapper, so a `WriteError` arrives as a readable JS `Error`; a bound *member* function surfaces as a bare, message-less `WebAssembly.Exception`. It matches the [C API](/c_api)'s `mio_xdmf_series*` too, so the two flat bindings describe the same object the same way.
:::

This is not a registry format (there is no single `(path, mesh)` call for the registry to name), so it does **not** change the format counts below and does not appear in `availableFormats()`.

## Format support

**As of v8.0.0 the WASM build ships every format the C++ core has**, including the five that need HDF5 or netCDF — there is no longer a WASM-specific format gap. That is 81 readable and 76 writable format keys in the lean pair, and 83 readable in the [full pair](#lean-and-full-builds), which also reads `step` and `iges` (in the lean pair `gmsh22` and `step` are write-only; `ansys_rst`, `ansys_rst_cyclic` and `lsdyna_d3plot` read-only, each for the reason the [conformance page](./conformance.md#formats-missing-a-direction) gives; `libmesh` gained its writer and `radioss_anim` joined in v16.11.0, `lsdyna_binout` and `radioss_th` in v16.12.0, `radioss` and `marc` gained writers in v16.17.0, `nastran_op2`, `lsdyna_binout` and `abaqus_fil` after them). The two formats of v16.13.0 are not in it: `vtx` needs ADIOS2 and `szplt` TecIO, neither of which the Emscripten build links, so both report the missing library by name. `gid` gained a reader in v10.19.0 and is now both. `.vts`/`.vtr`/`.vtm` (v11.6.0, tier B4) added three of each, none needing a `js_bindings.cpp` change at all -- see the format list above. `.pvd`/`.pvtu`/`.pvtp` (v15.0.0) added three more of each, again with no `js_bindings.cpp` change: they are registered in `registry.cpp` and reached through the existing read options. `.pcd` and `.xyz` (v15.1.0, point clouds) added one of each the same way; `.k`/`.key`/`.dyn` (v15.2.0, LS-DYNA keyword decks) did too, with `*INCLUDE` files resolved next to the deck in MEMFS; `.frd` (v15.3.0, CalculiX results) added one **read-only** key the same way, its increments reachable through the existing `timeStep` read option; `.h5` (v15.7.0, MSC Nastran HDF5 results) added another **read-only** key, on the wasm32 HDF5 already linked, its result domains reachable through `timeStep` too; `.mail` (v16.0.0, Code_Aster meshes) added one of each; `.mphbin` (v16.1.0, COMSOL binary meshes) added one of each; `elmer` and `.feb` (v16.2.0, Elmer mesh directories and FEBio input) added one of each, the Elmer writer creating its directory in MEMFS and the reader finding it by sniffing the directory, and `.xplt` (FEBio plot files) one **read-only** key, its states reachable through `timeStep` and its compressed states inflated by Emscripten's zlib port; `.rst`/`.rth` (v16.3.0, Ansys MAPDL results) added one more **read-only** key, its result sets reachable through `timeStep`; `.xda`/`.xdr`, `.fil` and `.rad` (v16.7.0: libMesh, Abaqus results, OpenRadioss) added three **read-only** keys, the `.fil` increments reachable through `timeStep`, and Z88 one of each, found by its fixed file name; `marc` (`.dat` files that open as MSC Marc decks), `.t19` (Marc post files, increments through `timeStep`) and `ansys_rst_cyclic` (by name) (v16.8.0) added three more **read-only** keys; `lsdyna_d3plot` (a family written into MEMFS as `d3plot`, `d3plot01`...: found by its file name, its states through `timeStep`) and `.op2` (Nastran OP2 results, steps through `timeStep`) (v16.9.0) added two more; `.pcd` is also the one format whose LZF stream (`binary_compressed`) is decoded in WASM with no dependency at all. `openfoam` became writable in v9.20.0; `gid` (the GiD postprocess format, on a vendored gidpost) needed no WASM-specific binding at all — it reaches every registry-driven consumer, including this one, purely by being registered in the shared `registry.cpp` writer/extension tables, and `-sUSE_ZLIB=1` supplies the `<zlib.h>` its vendored C sources need at the compiler (not just the linker).

> **Four Python-only formats are absent here.** `pmsh`, `zarr`, `cae` and `usd` (v10.35.0, the physics-ML data path — see [formats](./formats.md)) are registered in the Python layer only, not in the shared C++ dispatch registry, so this surface does not see them. Two of the four are directory stores and one needs `usd-core`; all are export targets for a training pipeline rather than interchange formats. `neuroglancer` has a native reader and writer (unreleased) and is in the list below, with no extension, so pass its format name.

`abaqus`, `abaqus_fil` (Abaqus results `.fil`, written since the `.fil` writer joined; every increment is a step, see [Abaqus `.fil`](./formats/abaqus_fil.md)), `ansys`, `ansys_rst` (Ansys MAPDL results, written since the writer joined; every result set is a step, see [Ansys results](./formats/ansys_rst.md)), `ansys_rst_cyclic` (the full rotor of a static cyclic `.rst`, **read-only**, by name), `ansysinp` (read/write; `CMBLOCK` components become regions, see [Ansys coded database](./formats/ansysinp.md)), `avsucd`, `cgns`, `code_aster` (Code_Aster `.mail`; `GROUP_MA`/`GROUP_NO` become regions, see [Code_Aster](./formats/code_aster.md)), `dex`, `dolfin`, `elmer` (Elmer mesh directories; written into MEMFS, read back with no format by sniffing the directory — see [Elmer](./formats/elmer.md)), `ensight` (EnSight Gold geometry, `.case`/`.geo`, ASCII + C-binary, and Fortran binary read since v16.17.0), `exodus`, `femap` (Femap neutral `.neu`; every output set is a step, see [Femap](./formats/femap.md)), `febio` (FEBio `.feb` meshes; surfaces become side regions, see [FEBio](./formats/febio.md)), `flac3d` (including `ZGROUP`/`FGROUP` cell groups, which ride on the mesh object as named regions — see [FLAC3D](./formats/flac3d.md#data-mapping)), `flux`, `frd` (CalculiX results, ASCII or binary written too; every increment is a step, see [CalculiX results](./formats/frd.md)), `freefem`, `gid` (GiD postprocess; *writing* rides a vendored gidpost and needs zlib, *reading* is meshio++'s own code and needs nothing for the ascii flavour), `gltf` (`.glb`/`.gltf` — written as the surface in glTF 2.0 with the default options, and read back by a general glTF 2.0 reader, see [glTF](./formats/gltf.md)), `gmsh`, `h5m`, `hmf`, `ip`, `libmesh` (libMesh `.xda`/`.xdr`, read and, since v16.11.0, written; gzip is read in both artifacts, bzip2 read and written in the [full](#lean-and-full-builds) one; see [libMesh](./formats/libmesh.md)), `lsdyna` (LS-DYNA keyword decks; `*PART`/`*SET_*` become regions, see [LS-DYNA](./formats/lsdyna.md)), `lsdyna_d3plot` (LS-DYNA state databases, **read-only**; the file named `d3plot` and its numbered family, every state a step, see [LS-DYNA d3plot](./formats/lsdyna_d3plot.md)), `marc` (MSC Marc input decks; a `.dat` that opens as a Marc deck reads as one, Tecplot's otherwise; written by name since v16.17.0, `writeMesh(path, mesh, 'marc')`, since a `.dat` alone writes Tecplot — see [Marc](./formats/marc.md)), `marc_t19` (Marc formatted post files `.t19`, written too; every increment is a step), `mdpa` (Kratos; mesh-level blocks only — see [MDPA](./formats/mdpa.md#c-core)), `med`, `medit`, `mfem` (MFEM `.mesh`; a `.mesh` whose first line names an MFEM mesh is read as one, see [MFEM](./formats/mfem.md)), `mff`, `mfm`, `mphbin`, `neuroglancer` (precomputed mesh fragments, triangles only, no extension: pass the format name), `mphtxt` (COMSOL meshes, binary and text; Selections become regions, see [COMSOL](./formats/mphtxt.md)), `nastran` (any bulk-data deck; HyperMesh components and OptiStruct `SET`s become regions, see [Nastran](./formats/nastran.md)), `nastran_h5` (MSC Nastran HDF5 results, written too; every result domain is a step, see [MSC Nastran HDF5](./formats/nastran_h5.md)), `nastran_op2` (Nastran OP2 results; every subcase, mode or time is a step, written too, see [Nastran OP2](./formats/nastran_op2.md)), `netgen`, `obj`, `off`, `openfoam` (read/write; writing creates a `constant/polyMesh` directory in MEMFS), `patran` (Patran 2 neutral `.pat`/`.out`; components become regions, see [Patran](./formats/patran.md)), `pcd` (PCL point clouds; reads `ascii`/`binary`/`binary_compressed`, writes `ascii`/`binary`), `permas`, `ply`, `pvd` (a time-indexed collection of `.vtu` steps), `pvtp`, `pvtu` (parallel indices over per-part pieces; index and pieces sit side by side in MEMFS), `radioss` (OpenRadioss starter `_0000.rad`, read and, since v16.17.0, written; `#include` files resolved next to the deck in MEMFS, see [Radioss](./formats/radioss.md)), `radioss_anim` (OpenRadioss animation files `<run>A001`…, v16.11.0, written since; see [OpenRadioss animation files](./formats/radioss_anim.md)), `radioss_th` (OpenRadioss time-history files `<run>T01`, read v16.12.0 and written since from the `radioss_th:*` field data, each output a step through `timeStep`; see [OpenRadioss time-history files](./formats/radioss_th.md)), `lsdyna_binout` (LS-DYNA's binout, read v16.12.0, written since; see [LS-DYNA binout](./formats/lsdyna_binout.md)), `step` (a closed surface written as a faceted CAD solid in both artifacts, and read, as its triangulation and B-rep, in the [full](#lean-and-full-builds) one; see [STEP](./formats/step.md)), `iges` (read in the full artifact only, see [IGES](./formats/iges.md)), `stl`, `su2`, `svg` (2D drawing; reads back only the drawings meshio++ wrote, as planar meshes, and refuses an unmarked file — see [SVG](./formats/svg.md)), `tecplot` (`.dat`/`.tec`, and binary `.plt` read-only), `tetgen`, `tikz` (reads back only the figures meshio++ wrote, as planar meshes, and refuses an unmarked file — see [TikZ](./formats/tikz.md); 2D LaTeX visualization), `triangle` (`.poly` by default; see the ambiguous-extensions table for `.node`/`.ele`), `ugrid`, `unv` (`.unv`/`.uff`; groups become regions, results are steps — see [UNV](./formats/unv.md)), `vti` (VTK XML ImageData: a regular lattice, so writing needs one — see [VTI](./formats/vti.md)), `vtk`, `vtkhdf` (HDF5-backed: UnstructuredGrid, PolyData and the two composite types, including polyhedra, partitions and time; `piece`/`timeStep` select what to read), `vtp`, `vtu` (zlib compression works via Emscripten's built-in port), `wkt`, `xdmf` (XML, Binary **and** HDF), `xplt` (FEBio plot files, written too; every state is a step, see [FEBio plot files](./formats/xplt.md)), `xyz` (headerless ASCII point clouds, see [XYZ](./formats/xyz.md)), `z88` (Z88 `z88i1.txt`, found by file name; its `z88o2.txt`/`z88o3.txt` results are read from beside it, see [Z88](./formats/z88.md)). The three field-only formats (`dex`, `ip`, `mff`) read/write geometry-less meshes (field values in `point_data`).

Ask the loaded module rather than trusting this list — it is generated from the same registry the build actually links:

```javascript
const { readers, writers } = m.availableFormats();
```

`readMesh(path, 'gmsh')` reads non-periodic Gmsh 2.2, 4.0 and 4.1. The 4.0 path supports ASCII and binary with either 4- or 8-byte producer counts, independent of WASM's own integer widths. For `$Periodic`, use `readMeshSelective(path, {format: 'gmsh', info: true})`; info-less reads refuse the section rather than lose it. `mesh.info.periodic` contains raw `entity` tags `[dimension, slave, master]`, `affine` (0 or 16 coefficients) and flattened `(N, 2)` `nodePairs` as **0-based point rows**, keeping ordering and duplicates. `writeMesh` automatically restores Gmsh info to 4.1 (`gmsh`) or 2.2 (`gmsh22`); callers may also supply `{info}` explicitly. Operations do not remap the channel, so rebuild it after point/topology edits. See [Gmsh](formats/gmsh.md).

### What the HDF5/netCDF formats cost, and what changed

- **The `.wasm` is ~6.3 MB sequential / ~6.7 MB threaded, up from ~2.3 MB before HDF5/netCDF** (the published npm tarball is roughly a third of that, since the binary compresses about 3:1 and browsers fetch it compressed). libhdf5, libnetcdf and libcgns are statically linked, and they are real bytes — though cgnslib itself is the cheap one, adding about **290 KB** (measured, v9.22.0), because it is a thin layer over the HDF5 already there. If you do not need these formats, build your own artifact with `./build/configure-wasm.sh --without-hdf5 --build`, or keep HDF5 and drop only the CGNS MLL with `--without-cgnslib` (see below) — the JS API is identical either way, and `availableFormats()` reports the smaller set while `hasCgnslib()` reports whether the MLL is present.
- **Breaking: `.xdmf` now writes an HDF companion file.** The registry's XDMF writer default follows the build, exactly as it does natively: with HDF5 present it emits `Format="HDF"` heavy data into a sibling `<base>.h5` and leaves only the XML skeleton in the `.xdmf`. A caller that used to pull one file out of the virtual filesystem must now pull **two**. Reading is unaffected (all three data formats are read).
- **MED writes plain fields, but not the enhanced ones.** The C++ MED writer handles ordinary `point_data`/`cell_data` — the single-timestep `CHA` common case, including multi-component fields since v9.9.0 — so a field-carrying mesh writes fine here. What it cannot do is the constructs the Python reference writer alone implements (units, multi-timestep metadata, component names supplied via `med:nom`, named profiles, ELNO/ELGA), and this build has no Python to defer to, so those throw by name. MED geometry, `point_tags`/`cell_tags`, families, named regions and `med:num` global numbering are all written normally. See [MED quirks](./formats/med.md#quirks-limitations).
- **`sniffFormat` still cannot identify them.** A plain HDF5 magic number says nothing about which of `med`/`cgns`/`h5m`/`hmf`/`vtkhdf`/XDMF-HDF a file is, so the sniffer deliberately never claims it. Use the extension, or pass `format` explicitly.

### Ambiguous extensions

Some extensions are shared by more than one format. `readMesh`/`writeMesh`/ `convert` all take an optional trailing `format` argument (or an `{inFormat, outFormat}` options object for `convert`) to disambiguate, mirroring Python's `file_format=` kwarg:

| Extension | Default format | Pass `format=` to select instead |
|-----------|-----------------|-----------------------------------|
| `.msh` | `gmsh` | `"ansys"`, `"freefem"` |
| `.inp` | `abaqus` | `"ansysinp"` |
| `.node` / `.ele` | `tetgen` | `"triangle"` (2D Triangle pairs) |

## Threads (OpenMP)

The package ships **two** native artifacts and picks one at load time:

- `meshioplusplus_wasm` — the **sequential** build. Always loadable, anywhere.
- `meshioplusplus_wasm_mt` — the **threaded** build, compiled with the OpenMP parallel backend over Emscripten's Wasm threads (pthreads + `SharedArrayBuffer`). The core's parallel loops — every [mesh operation](#mesh-operations) run through `convertSurfaceOps`, and VTU zlib compression — run multi-threaded here.

`loadMeshioPlusPlus()` **auto-selects**: the threaded build under Node (where `SharedArrayBuffer` is always available) and in a **cross-origin-isolated** browser context, the sequential build otherwise. Force one with the `variant` option:

```js
const meshio = await loadMeshioPlusPlus({}, { variant: "mt" });  // "auto" (default) | "mt" | "seq"
meshio.parallelBackend();  // "openmp" for the threaded build, "seq" for the sequential one
```

**A browser page must be [cross-origin isolated](https://developer.mozilla.org/en-US/docs/Web/API/crossOriginIsolated) to load the threaded build.** `SharedArrayBuffer` requires the document to be served with:

```
Cross-Origin-Opener-Policy: same-origin
Cross-Origin-Embedder-Policy: require-corp
```

A threaded module **cannot instantiate** without this — and there is no in-artifact fallback — which is exactly why both variants ship and the loader chooses. Where you cannot set those headers (e.g. static hosts like GitHub Pages), a COOP/COEP [service worker](https://github.com/gzuidhof/coi-serviceworker) is the usual workaround; it is what the meshio++ browser viewer's Pages demo uses. If a page is not isolated, `loadMeshioPlusPlus()` transparently loads the sequential build instead, so correctness never depends on the headers — only speed does.

When you pass a `locateFile` override to relocate the `.wasm` (a bundler/CDN setup), return the URL matching the requested filename, since the two variants ask for different binaries:

```js
const meshio = await loadMeshioPlusPlus({
    locateFile: (path) => (path.includes("_mt") ? mtWasmUrl : wasmUrl),
});
```

Under Node no headers are needed — Wasm threads use `worker_threads`. The threaded artifact pre-spawns a worker pool of `navigator.hardwareConcurrency` (falling back to 8 where `navigator` is absent, e.g. Node < 21, then growing on demand).

If instantiation fails — a missing `.wasm` file, a `locateFile` that resolved to the wrong URL, or an environment that cannot host the requested variant — `loadMeshioPlusPlus()` rejects with a `MeshioPlusPlusLoadError` (extends `Error`) instead of a bare Emscripten abort. It carries `variant`, `glue`, `requestedFile` and `resolvedUrl` (what Emscripten asked `locateFile` for and what it returned) and `cause` (the underlying error). A `locateFile` that returns the sequential `.wasm` for the threaded glue (or vice versa) can link successfully and only then report the wrong `parallelBackend()`; `loadMeshioPlusPlus()` checks for this and rejects with the same error type rather than returning a module that silently runs the other variant.

## Known v1 limitations

- **No zero-copy.** Every array is copied once crossing the JS/WASM boundary (see above) — for very large meshes this has a real memory/time cost that the Python bindings' numpy views avoid.
- **Per-format write options are limited to encoding/codec/float format.** `writeMesh`/`convert` accept `{encoding, codec, floatFormat}` (see "Write options" below); there is no gzip *level* control or a VTK 4.2 vs 5.1 selector, since `WriteOptions` has neither on the C++ side (gzip level 4 is a fixed registry default; `vtk42`/`vtk51` are separate format keys, not a `vtk` option). Omitting `options` entirely reproduces the exact pre-v11.2.0 default for every format (e.g. `vtu` binary+zlib, `stl` ASCII, `gmsh` the 4.1 binary format; `stl`/`ply` extract and write the boundary **skin** of a 3D volume mesh — the Python default, see [Skin extraction](./extract_skin.md) — and `svg`/`tikz` render 3D meshes with the default isometric camera).

### Write options

`writeGltf(path, mesh, options)` exposes all native glTF options independently of generic `writeMesh`: `{colorBy, cmap, component, vmin, vmax, splitAngle, upAxis, scale, container, normalWeight, normals, fields, recenter, byRegion, unlit, nanColor}`. Components are 0-based, omitted selects magnitude, and the return value lists every written MEMFS path (`.gltf` plus `.bin`, or `.glb`). For example, `meshio.writeGltf('/colored.glb', mesh, {colorBy: 'temperature', cmap: 'turbo', upAxis: 'z', scale: 0.001})`. Native reads also accept XDMF2/3 with absolute DataItem references and Netgen name/periodic tables; Netgen periodic arrays remain in `field_data` with 1-based file node ids, and `.vol.gz` needs a zlib-enabled build.

`writeMesh(path, mesh, format, options)` and `convert(inPath, outPath, {..., encoding, codec, floatFormat})` take an options object honoured through the same core entry point (`registry_write_ex`) the C API and native CLI use:

```js
meshio.writeMesh('/out.vtu', mesh, '', { encoding: 'ascii', codec: 'none' });
meshio.convert('/in.vtu', '/out.stl', { encoding: 'binary' });
```

- `encoding` — `'ascii'` or `'binary'`; throws for a format with only one variant (`frd` has both: `'binary'` writes the layout `ccx` uses for `*NODE OUTPUT`, see [CalculiX results](./formats/frd.md#binary-output)). `'raw_appended'` (since v16.21.0) writes a `vtu`'s arrays as raw binary in one `<AppendedData>` section and throws for any other format.
- `codec` — `'none'`, `'zlib'`, `'lz4'`, or `'zstd'` for VTK-XML (`vti`/`vtu`/`vtp`), or `'lzf'` for PCD `binary_compressed`. LZF needs no optional library; combining it with ASCII or a non-PCD format throws.
- `floatFormat` — a `printf`-style float format for ASCII writers that take one (e.g. `".16e"`, the default).

An option the target format cannot honour throws a catchable `Error` naming the format and the option, rather than silently ignoring it — the same rule the CLI's own `--encoding`/`--codec`/`--float-format` flags follow.

Both calls **return every virtual-FS path the write touched** (new or changed since the call started), sorted: `['/out.vtu']` for a single-file writer, `['/out.h5', '/out.xdmf']` for `xdmf`'s HDF companion, or the five `constant/polyMesh/*` files (plus the `case.foam` marker) for `openfoam`. This is computed by diffing the output's directory before and after the write (there is no core-level "list of paths written" yet — a documented roadmap remainder), so it is exact for any present or future writer with no binding change, at the cost of a directory walk per call.
- **Named regions are carried** (since v8.1.0). They ride on the mesh object itself — `mesh.regions` is an array of `{ name, kind, dim, tag, entries }` — so `readMesh` / `writeMesh` / `convert` carry them with no extra call, and nothing new had to be forwarded by the wrapper. `kind` is `'point'`, `'cell'` (global block-major cell indices) or `'side'` (`(cell, facet)` pairs). The Phase-1 formats (gmsh, abaqus, and MED since v9.6.0 — one region per `FAS`/`GRO` group name) map onto them fully; Exodus reads but does not yet write them. See [Named regions](./regions.md).
- **Side-channel data is exposed since v11.2.0.** `openfoam`'s patch names/types, MED field units/step metadata, MDPA properties/entity names, and the `ansysInp`/`unv` set channels all reach JS via `readMeshSelective(path, {info: true})` and `writeMesh(path, mesh, format, {info})` — see [Side channel (info)](#side-channel-info) above.
- **Data arrays are `Float64` or `BigInt64`, not their exact source width.** Since v11.2.0 an integer material id no longer becomes a double crossing the boundary (see "Data array dtypes" above), but the mesh's own storage still canonicalizes *within kind* rather than preserving e.g. `Int32Array` vs. `Uint8Array` exactly — every integer dtype comes back as `BigInt64Array`, every float dtype as `Float64Array`. Multi-component (vector/tensor) arrays *are* supported since v9.9.0 via the `*_components` objects (see "The mesh object shape" above); before that they were flattened to N unrelated scalars in both directions. The path-based `convert`/`convertSurface` calls still avoid the boundary entirely, and so preserve exact dtypes as well as shapes.

## Building from source

Requires the [Emscripten SDK](https://emscripten.org/docs/getting_started/downloads.html):

```sh
git clone https://github.com/emscripten-core/emsdk.git
cd emsdk && ./emsdk install latest && ./emsdk activate latest
source ./emsdk_env.sh
cd ../meshioplusplus  # this repo
./build/configure-wasm.sh --build
node tests/wasm/smoke.mjs
```

`build/configure-wasm.sh --build` builds **both** artifacts by default (`--seq-only` skips the threaded one). Each is configured with `-DMESHIOPLUSPLUS_BUILD_PYTHON=OFF` (no Python/pybind11 involved) and `-DMESHIOPLUSPLUS_MESH_BACKEND=NATIVE` (the fastest [in-memory mesh backend](cpp_backends.md) — canonical Float64/Int64 storage, so the embind typed-array boundary needs no dtype dispatch; the JS API shape is unchanged, and `meshBackend()` on the loaded module reports `"native"`). The two differ only in the parallel backend: the sequential variant uses `-DMESHIOPLUSPLUS_PARALLEL_BACKEND=SEQ`, the threaded one `-DMESHIOPLUSPLUS_PARALLEL_BACKEND=OPENMP -DMESHIOPLUSPLUS_WASM_THREADS=ON` (which adds `-pthread` and a pre-spawned worker pool — see [Threads (OpenMP)](#threads-openmp)). They build in separate trees (`build/wasm-<type>` and `build/wasm-<type>-mt`) because `-pthread` is a whole-translation-unit property. See `--help` for `--without-zlib`, `--without-hdf5`, `--without-netcdf`, `--deps-prefix` and `--build-type`. CI (`.github/workflows/wasm.yml`) builds and smoke-tests both on PRs touching the wasm surface, and publishes to npm on `v*` tags.

`--full` adds the [full pair](#lean-and-full-builds) (`--full-only` builds only it; `--without-occt` leaves OpenCASCADE out of it, for a codecs-only full build that stays MIT), configured with `-DMESHIOPLUSPLUS_WASM_FULL=ON` in its own trees (`build/wasm-<type>-full[-mt]`). Only under that option does CMake look for OpenCASCADE, zstd and lz4 on a wasm target, in the prefixes `build/build-wasm-deps.sh --full` builds, one per threading flavour (OpenCASCADE's reference counts are `std::atomic`, and an object built without `-pthread` cannot join the threaded link); bzip2 is Emscripten's port, like zlib. A full build never drops a library it was asked for: one not found is a configure error, not a smaller artifact. OpenCASCADE is a long build, most of an hour per flavour on four cores; CI caches its prefix apart from the lean one, and the release job assembles the package from both builds' artifacts.

```sh
./build/build-wasm-deps.sh --full                       # lean deps + OpenCASCADE, zstd, lz4 (seq and mt)
./build/configure-wasm.sh --full --build                # all four artifacts
MESHIOPLUSPLUS_WASM_VARIANT=full node tests/wasm/smoke.mjs
```

The deps script fetches release tarballs by checksum and the OpenCASCADE tag by a shallow clone verified against its commit; netCDF and CGNS fall back to such a clone of their tags when GitHub's archive host is unreachable, as some proxies make it.

### The HDF5 and netCDF dependencies

Neither is a system library on this target, and meshio++'s CMake never downloads anything — so `build/build-wasm-deps.sh` produces them, and CMake only *finds* the result:

```sh
./build/build-wasm-deps.sh                 # HDF5 1.14.6 + netcdf-c 4.9.3 -> a prefix
./build/build-wasm-deps.sh --print-prefix   # where that prefix is
```

`configure-wasm.sh` runs it automatically the first time (several minutes) and passes the prefix to CMake via `CMAKE_FIND_ROOT_PATH` — the required knob, since the Emscripten toolchain re-roots every `find_package` at its sysroot, which makes `CMAKE_PREFIX_PATH` and `HINTS` alone ineffective. Pass `--deps-prefix` to reuse a prefix you already have.

Both libraries are built static, `-Oz`, against Emscripten's own zlib port, with everything that assumes an OS this target does not have switched off: HDF5 without threadsafe/MPI/plugins and without the ROS3, direct, mirror and subfiling VFDs; netCDF without DAP, DAP4, byterange, NCZarr, S3, libxml2 (it falls back to its bundled `ezxml`), the `dlopen` plugin loader and mmap. Two upstream-neutral fix-ups are applied and documented in the script: a one-hunk patch guarding HDF5's `feclearexcept(FE_INVALID)` (wasm defines no floating-point exception flags at all), and a rewrite of the exported CMake link interfaces, which name imported targets (`ZLIB::ZLIB`, `HDF5::HDF5`, `hdf5::hdf5_hl`) that no single `find_package` in a consumer defines.

**cgnslib is linked in too** (v9.22.0), cross-compiled by the same `build/build-wasm-deps.sh` as a third pinned + SHA256-checked source build. It is strictly an *addition*: meshio++ reads and writes CGNS itself over raw HDF5, including polyhedral `NGON_n`/`NFACE_n` sections since v9.21.0, so what the MLL buys here is **ADF-backed containers** — which are not HDF5 at all, and so are unreachable from the hand-rolled path by construction — and the **CGNS 3.x** section layout. `--without-cgnslib` drops it and keeps everything else.

**The wasm stack is sized for them.** `CMakeLists.txt` links the wasm target with `-sSTACK_SIZE=4MB`, well above Emscripten's 64 KiB default, because HDF5's and netCDF-4's frames overrun it. That overrun is silent — the stack grows down into the static data segment — and cost a real bug during development: one Exodus write clobbered libc++'s locale facets, after which every ASCII reader in the module trapped. `-sSTACK_OVERFLOW_CHECK=1` is on so a recurrence aborts loudly instead.

## Selective reads, metadata, and codecs

`readMeshSelective(path, { format, pointsOnly, arrays, timeStep })` and `readMetadata(path, format)` mirror the Python API; `readerSupportsOptions(format)` reports whether a format has a native selective path. `arrays: null` reads every data array, `arrays: []` reads none.

```javascript
const mesh = m.readMeshSelective('big.vtu', { arrays: ['u'] });
const meta = m.readMetadata('big.vtu');   // meta.fellBackToFullRead

// A multi-step file: 0 (default) is the first step, negative counts from the end.
const last = m.readMeshSelective('run.exo', { format: 'exodus', timeStep: -1 });
m.readMetadata('run.exo', 'exodus').timeValues;  // [0, 0.5, 1] -- always present
```

**Memory mapping is unavailable** here — the Emscripten virtual filesystem has nothing to map, so `FileSource` always uses buffered reads. The option is accepted and ignored.

**zstd and lz4 are in the full artifact only.** The lean pair keeps zlib (`-sUSE_ZLIB=1`) as its one block codec, unchanged and still the default, so every file a WASM build wrote before still round-trips; `codec: 'zstd'` or `'lz4'` with `encoding: 'binary'` there throws naming the codec (the codec applies to binary blocks only, and a VTK XML write is ASCII unless `encoding` says otherwise). The [full pair](#lean-and-full-builds) links both (zstd 1.5.7 and lz4 1.10.0, cross-compiled by `build/build-wasm-deps.sh --full`), and `hasZstd()`/`hasLz4()` say which artifact you have. A zstd-compressed VTK XML file is a meshio++ extension (VTK has no zstd compressor), so vtk.js and ParaView do not read it; meshio++'s native, Python and full WASM readers do. lz4 is VTK's own `vtkLZ4DataCompressor`.

## Software rendering (v16.34.0)

`render` draws a mesh with a deterministic software rasterizer: **no WebGL**, so it works in a Node script or a server-side preview where the [browser viewer](/viewer) cannot. See [terminal rendering](/tui) for what is drawn and every option; the options are one object with camelCase names (`colorBy`, `scaleBar`, `isoLevels`, ...), and an unknown key is an error listing the ones that exist.

```js
const frame = m.render(mesh, { width: 800, height: 600, supersample: 2, colorBy: 'temperature', colorbar: true, clip: [2, 98] });
const image = new ImageData(frame.rgba, frame.width, frame.height);   // rgba is a Uint8ClampedArray
frame.cellIds;   // Float64Array: the input cell at each pixel, -1 for none
frame.range;     // [vmin, vmax], or null when no field is mapped; frame.notes holds the ticks and keys

console.log(m.renderText(mesh, { colorBy: 'temperature' }, { cols: 80, rows: 24, encoding: 'braille' }));
const png = m.renderPng(mesh, { width: 800, height: 600 });           // Uint8Array
```

**Cut-aways (v16.36.0).** `cutaway: [px, py, pz, nx, ny, nz]` (a point and the normal of the side kept), or an array of up to two such arrays, with `cutawayTint`, clips the geometry away and tints the back faces now in view: `m.render(mesh, { view: '-x', cutaway: [0.5, 0, 0, 1, 0, 0] })`.

**Streamlines (v16.37.0).** `streamlines: 'velocity'` follows a vector point array over the mesh's cells, with `streamSeeds` (default 40), `streamLength` (model diagonals each way, default 0.5) and `streamColor`: `m.render(mesh, { streamlines: 'velocity', streamSeeds: 60 })`. The point array is a flat typed array of three numbers per point. An unknown array or a seed count outside 1 to 10000 throws, naming it.

Colours are `0xRRGGBBAA` numbers or `[r, g, b, a?]`. `renderText` takes the cell encodings only (`halfblock`, `quadrant`, `sextant`, `braille`, `ascii`): the Kitty, iTerm2 and Sixel protocols and the interactive terminal loop are CLI features and are not exported, so nothing is added to the bundle but the rasterizer itself.

## v9.1.0 additions

- `readMeshSelective(path, { lenient: true })` — see [`doc/selective_read.md`](selective_read.md).
- `readMeshSelective(path, { piece: k })` — keep one piece of a partitioned file (VTKHDF partitions or composite blocks): `0` is the first, negative counts from the end, `null` (the default) merges every piece into one mesh with one `cell` region per piece; out of range throws naming the piece count. See [`doc/selective_read.md`](selective_read.md#picking-a-piece).
- `readMeshSelective(path, { dropGhosts: true })` — remove the ghost cells (halo) of a partitioned `.pvtu`, `.pvtp` or `.pvd`: every cell with a `vtkGhostType` bit set, and the points only they used. `false` (the default) keeps them; every other reader ignores it. See [`doc/selective_read.md`](selective_read.md#dropping-ghost-cells).
- XDMF series: `w.flush()`, `w.writeDataArrays(time, pointData, cellData, components)`, and `createXdmfTimeSeriesWriter(path, { mode: 'append', autoFlush: true })`. With `flush()` the `.xdmf` appears in MEMFS before `finalize()`, so a partially-written series can be copied out.

## Sequences (transient / multi-file datasets)

A set of MEMFS files — or the steps inside one multi-step file — treated as one ordered dataset, the same surface `convert` and `runPipeline` already work on. See [sequences](sequences.md) for the ordering rule, the time-value precedence and the streaming guarantee.

```js
// Stage the steps in MEMFS, then treat them as one dataset.
for (let i = 0; i < 12; ++i) m.writeMesh(`/seq/out_${i}.vtu`, meshes[i]);

const plan = m.sequenceEntries('/seq/out_*.vtu');
// -> [{path, step, time, timeSource}, ...] in NATURAL-NUMERIC order, so
//    out_9.vtu precedes out_10.vtu (a plain sort gets that backwards).

m.sequenceToTimeseries('/seq/out_*.vtu', '/seq/series.xdmf');   // fan-in  -> 12
const bytes = m.FS.readFile('/seq/series.xdmf');

const paths = m.timeseriesToSequence('/seq/series.xdmf', '/seq/back_{step}.vtu');
// -> ['/seq/back_0000.vtu', ...] fan-out; read them back with Module.FS

// Onto new times (linear blend of the two bracketing steps; at most two live):
m.resampleSequence('/seq/out_*.vtu', '/seq/u_{index}.vtu',
                   { Times: { Start: 0, Stop: 11, Step: 0.5 }, Method: 'linear' });  // -> 23
const half = m.blendSteps(meshes[0], meshes[1], 0.5);   // one blended step

// A whole transient post-processing run: the chain applies to EVERY step.
m.runPipeline({
  Version: 1,
  Input: { Pattern: '/seq/out_*.vtu' },
  Operations: [{ Op: 'Quality' }],
  Output: { Path: '/seq/post_{step}.vtu' },
});
```

`runPipeline` **routes** a transient document (a `Pattern`/`Paths` input, a `{step}`/`{index}` output, or `Mode`/`Parallel`/`Workers`) to the sequence driver automatically; a plain single-mesh document takes the unchanged path, so nobody has to know which kind of document they hold.

The pattern language is deliberately just `*` and `?` — **no** `**`, no `[set]`, and the directory component is literal — identical to the core's, so the browser and the CLIs accept exactly the same words.

Three things that fail **by name** rather than doing something surprising: a fan-in to a format that cannot hold a series (only XDMF can), a fan-out without a `{step}`/`{index}` token, and a multi-step input aimed at a single-step output. None of them silently keeps step 0.

**`Parallel` is accepted and ignored with a warning** in the report: it is a Python-driver feature (a thread pool), which this build does not run. The steps run in order, which is what the streaming guarantee needs anyway.

### Reading a sequence step by step

`openSequence(source, options)` is `sequenceEntries`' **stateful** counterpart: it plans the sequence once (no heavy data read, same `source`/`options`), then lets you read one step's mesh at a time — at most one mesh alive, whatever the step count, the same streaming guarantee `sequenceToTimeseries`/`timeseriesToSequence` give you without a JS-side loop managing state of your own:

```js
const seq = m.openSequence('/seq/out_*.vtu');
seq.count;                    // 12
seq.entry(3);                 // {path, step, time, timeSource} -- entry(i) === entries()[i]
seq.read(3);                  // the full Mesh at entry 3
seq.read(3, { pointsOnly: true, arrays: ['temperature'], lenient: true });   // as readMeshSelective
seq.close();                  // release the handle; safe to call twice / from a `finally`
```

`read`'s options are `readMeshSelective`'s (`pointsOnly`/`arrays`/`lenient`); the entry's own step index is always used, so there is no separate `timeStep` option. Every method after `close()` throws a catchable `Error` naming the handle, not a use-after-free. `sequenceEntries` (a single call, the whole plan at once) is still the right choice when you want every entry's metadata up front without reading any mesh at all — `openSequence` is for when you also intend to read some or all of the meshes, one at a time.

## Regular grids and signed distance (v9.24.0, `computeSdf` v9.25.0)

`grid(dims, origin, spacing)`, `voxelize(mesh, resolution, ...)`, `surfaceWatertightCheck(mesh)`, `sampleDistance(surface, points)` and `distanceToSurface(query, surface)`. `grid` is the only binding in the package that takes no input mesh — it creates one. `sampleDistance` takes a **flat** `[x0,y0,z0, x1,y1,z1, …]` array and returns a `Float64Array`.

`computeSdf(surface, structure, resolution, ...)` does both halves in one call, returning `{ mesh, dims, origin, spacing, maxDepth, numBanded, quality }`. `structure: 'octree'` refines only near the surface, and sizes itself from `rootResolution`/`maxDepth` — passing `resolution` or `cellSize` with it is an error, since its finest cell is already determined. Its output is 1-irregular (it has hanging nodes).

`{ Op: 'Voxelize' }` and `{ Op: 'ComputeSdf' }` also work as pipeline / `convertSurfaceOps` steps, and are the only steps that *replace* their input's geometry rather than transforming it.

The generated grid carries an `sdf:*` numeric `field_data` header describing itself. No file format persists arbitrary `field_data`, so write it as [`.vti`](formats/vti.md), whose `Origin`/`Spacing`/`WholeExtent` attributes are the same information.

See [`doc/voxelize.md`](voxelize.md) and [`doc/sdf.md`](sdf.md).
