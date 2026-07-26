# WebAssembly / JavaScript

The C++ core also compiles to WebAssembly and ships as an npm package, [`@meshioplusplus/wasm`](https://www.npmjs.com/package/@meshioplusplus/wasm), for reading and writing meshes in the browser or Node.js. (It is one of two "flat" bindings over the same core and shared format-dispatch registry — the other is the [C API](/c_api).)

## Install

```sh
npm install @meshioplusplus/wasm
```

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

### Rendering a mesh

`convertSurface` produces something a surface renderer can draw: a mesh with 3D
cells becomes its boundary, anything else passes through, and the result is
linearized (a triangle renderer has no mid-side nodes, so `triangle6`
connectivity drawn verbatim is visible garbage). Boundary facets inherit their
owning cell's data, so a per-cell material or tag still colours correctly.

```js
meshio.convertSurface("/solid.msh", "/solid.vtp");
const vtp = meshio.FS.readFile("/solid.vtp");   // hand this to vtk.js
```

Prefer it over `readMesh` → `extractSkin` → `writeMesh` for anything headed to
a renderer: it never materializes a JS mesh, so **multi-component (vector and
tensor) arrays survive** — the flat JS representation cannot carry them (see
"Known v1 limitations" below).

### Applying operations first

`convertSurfaceOps` takes a pipeline and applies it before extracting the
surface — all inside C++:

```js
const report = meshio.convertSurfaceOps('/part.msh', '/part.vtp', [
  { op: 'clean', weld: true },
  { op: 'smooth', method: 'taubin', iterations: 20 },
  { op: 'quality' },
]);
console.log(report.steps, report.warnings);
```

Prefer it over chaining the individual operation bindings (`clean`, `smooth`,
…): each of those takes and returns a JS `Mesh`, so a pipeline built from them
destroys every multi-component array on the first step. An **empty** pipeline
is byte-identical to `convertSurface`, which is what lets a viewer use one code
path for the plain and the post-operation display — and makes undo a replay of
a shortened pipeline rather than a set of inverse operations.

This is exactly what the [browser viewer](./viewer.md) does. It is built on
this package and is worth reading as a worked example of the whole pipeline —
worker, transferable buffers, and vtk.js — as well as being a live client-side
format converter you can try at
**<https://loumalouomega.github.io/meshioplusplus/viewer/>**.

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
  point_data?: { [name: string]: Float64Array },
  cell_data?: { [name: string]: Float64Array[] },   // one array per cell block
  field_data?: { [name: string]: Float64Array },
}
```

This deliberately mirrors the Python `Mesh`'s structure (points, a list of cell blocks, `cell_data` as one array per block). Cell connectivity is always `Int32Array` — the C++ core's connectivity dtype is Int64, but node/point counts for any mesh a browser can reasonably hold fit comfortably in 32 bits, and `Int32Array` is far more ergonomic in JS than `BigInt64Array`.

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

`readMesh` reports whichever shape the source held; `writeMesh` accepts either, and the target format's own writer decides what it can represent — MED writes ragged polygons directly, but no C++ format writer accepts a polyhedron block yet (a documented, pre-existing gap distinct from the JS boundary itself, which does carry polyhedron blocks correctly through operations like `clean`/`merge`/`convert_cells`), so writing one throws naming the format rather than silently dropping data.

## Mesh operations

Besides read/write/convert, the module exposes the mesh operations on the same
mesh-object shape: `extractSurface(mesh [, recordParentIds])`,
`extractSkin(mesh [, linearize])`, `attachQuality(mesh)`, `sniffFormat(path)`,
`computeBandwidth(mesh)`, and `reorder(mesh, method)` (method `"rcm"`,
`"morton"`, or `"hilbert"`). `reorder` returns
`{ mesh, nodePermutation, cellPermutations }` — the renumbered mesh plus the
applied old→new permutations (`Int32Array`s). See
[Reordering / renumbering](./reorder.md). Comparison is exposed as
`meshesEqual(meshA, meshB, atol, rtol, unordered)` (boolean) and
`diff(meshA, meshB, atol, rtol, unordered)` (a report object with `verdict`,
`points`, `blocks`, `pointData`/`cellData`/`fieldData`); named point/cell sets
are not compared. See [Mesh comparison (diff)](./diff.md). Merging is exposed as
`merge([meshA, meshB, ...], weld, atol, sourceTag, dataPolicy, dropDuplicateCells)`
(`dataPolicy` is `"intersection"` or `"fill"`), returning a new mesh object;
named point/cell sets are not carried. See [Merge / combine](./merge.md). The
editing/stats bundle is exposed too: `transform(mesh, matrix, rotateVectorData)`
(a 16-element row-major matrix) → mesh; `clean(mesh, weld, atol, removeOrphans,
dropDegenerate, dropDuplicateCells)` → `{ mesh, pointsWelded, ... }`;
`cropBbox(mesh, lo, hi, mode, recordIds)` / `cropPlane(mesh, point, normal, mode,
recordIds)` (`mode` `"all"`/`"any"`) → mesh; `split(mesh, by, tagName)` → an
array of `{ key, mesh }`; and `stats(mesh)` → an object of geometric measures
(`bboxMin`/`bboxMax`/`extent`/`centroid`, `cellTypeCounts`, `totalArea`,
`signedVolume`, `unsignedVolume`, `numInverted`). Element-representation
conversion is exposed as `convertCells(mesh, mode, recordParentIds)` with `mode`
`"linearize"`, `"simplexify"`, or `"elevate"`, returning a new mesh; a
polyhedron block under `"simplexify"` and the full-Lagrange targets
(`quad9`/`hexahedron27`) under `"elevate"` throw a catchable `Error`. Uniform
refinement is exposed as `refine(mesh, levels, recordParentIds)`, subdividing
every cell into same-type children (`triangle`/`quad` into 4,
`tetra`/`wedge`/`hexahedron` into 8) with shared mid-entity nodes, so the result
has no hanging nodes; a higher-order cell, a `pyramid`, or a ragged block throws
a catchable `Error`. Partitioning is exposed as `partition(mesh, nparts, method,
imbalance, mode, seed, recordIds, ghostLayers, weightsKey)` → an array of
`{ partId, mesh }` (exactly `nparts` entries, blocks kept 1:1 with the input,
unlike `split`) and `partitionLabels(mesh, nparts, method, imbalance, mode,
seed, weightsKey)` → one label array per cell block. **Only the SFC method
exists in the WASM build** — KaHIP is never compiled in (no Emscripten port,
and it would bloat the bundle), so `method: "kahip"` throws a catchable `Error`
naming `MESHIOPLUSPLUS_WITH_KAHIP` and `"auto"` always resolves to SFC.
Smoothing is exposed as `smooth(mesh, method, iterations, lambda, mu,
fixBoundary, preserveFeatures, featureAngle, guardInversion)` with `method`
`"taubin"` (the default, shrink-free) or `"laplacian"` (stronger per pass, but
shrinking), returning `{ mesh, numNodesMoved, maxDisplacement,
numSkippedInversion }`; only the point coordinates move, so connectivity and
every data value come through unchanged. Boundary nodes, feature nodes and the
nodes of blocks with unknown edge topology are pinned by default; the caller
`frozen` mask of the C++ API is not exposed here. The two cutters are exposed as
`slice(mesh, origin, normal, recordParentIds)` — the planar cross-section, one
dimension below the cut cells — and `isosurface(mesh, array, isovalues,
component, recordParentIds)`, its data-driven sibling: the level set of a scalar
`point_data` array. `isovalues` accepts a number or an array (several contours
land in one mesh, tagged per cell with `iso:value` and `iso:index`), `component`
is negative for the row magnitude, and naming a `cell_data` array throws a
catchable `Error` — cell data is piecewise constant and has no level set, so
convert it with `dataCellToPoint` first. See
[transform](./transform.md), [clean](./clean.md),
[crop](./crop.md), [split](./split.md), [stats](./stats.md),
[cell conversion](./convert_cells.md), [refine](./refine.md),
[partitioning](./partition.md), [smoothing](./smooth.md),
[slicing](./slice.md), and [isosurfaces](./isosurface.md).

::: tip Reachable from `loadMeshioPlusPlus()` since v7.4.0
Before v7.4.0 the geometry operations above were bound in the WASM module but
**not forwarded by the package wrapper**, so they were unreachable through
`loadMeshioPlusPlus()` (only file I/O and the `data_*` operations were). They are
all forwarded now. The index maps the C++ core returns for
`cropBbox`/`cropPlane`/`split`/`convertCells`/`refine`/`partition` are still not
carried across the JS boundary — use the `recordIds`/`recordParentIds` flags,
which attach the same provenance as ordinary data arrays (or `partitionLabels`
for the raw assignment). `smooth` needs none of that — it never adds, removes or
renumbers a node or a cell.
:::

The [data operations](./data_operations.md) — which act on `point_data` /
`cell_data` / `field_data` and never modify the geometry — are exposed as
`dataDrop(mesh, location, names, ignoreMissing)`, `dataKeep(...)`,
`dataRename(mesh, location, from, to)`,
`dataPointToCell(mesh, names, suffix)`,
`dataCellToPoint(mesh, names, weight, suffix)`,
`dataCalc(mesh, expression, location, outputName, overwrite)`,
`dataCondition(mesh, location, names, mode, lo, hi, scope, nanPolicy, nanReplacement, suffix)`
— each returning a new mesh — and `dataInfo(mesh)`, which returns an array of
per-array summary objects (`location`, `name`, `dtype`, `shape`, `numBlocks`,
`numEntries`, `numComponents`, `numValues`, `min`, `max`, `mean`,
`minPerComponent`/`maxPerComponent`/`meanPerComponent`, `numNan`, `numInf`,
`numFinite`, `inconsistentBlocks`). Enumerations cross as strings: `location` is
`"point"`/`"cell"`/`"field"`, `weight` is `"uniform"`/`"measure"`, `mode` is
`"clamp"`/`"normalize"`/`"standardize"`, `scope` is
`"component"`/`"magnitude"`, and `nanPolicy` is
`"ignore"`/`"replace"`/`"fail"`. A malformed `dataCalc` expression throws a
catchable `Error`. See [data operations](./data_operations.md),
[array management](./data_manage.md), [averaging](./data_average.md),
[expressions](./data_calc.md), [conditioning](./data_condition.md) and
[data summary](./data_info.md).

## Transient (time-series) XDMF

`createXdmfTimeSeriesWriter(path, { dataFormat, gzipLevel })` is the **one
stateful** thing in this API: every other binding is a pure function over a
mesh object, but a time series writes the mesh **once** and then appends one
cheap step per solve, and its `.xdmf` light data can only be written when the
collection is complete. See [XDMF time series](./xdmf_time_series.md).

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

The returned object has `writePointsCells(mesh)`, `writeData(time, mesh)`,
`finalize()`, `numSteps()`, `finalized()` and `close()`. `close()` finalizes if
needed and releases the handle; it is safe to call twice and safe to call from a
`finally`. After it, every method throws a catchable `Error` — the underlying
handle is an index into a module-local table, not a pointer, so a stale handle
can never be a use-after-free.

::: warning `'HDF'` writes TWO files into the virtual filesystem
The `.xdmf` at the path you gave **and** its sibling heavy-data file
`<path minus extension>.h5` (a sibling of the `.xdmf`, not a file in the
current directory — this differs from the pure-Python `TimeSeriesWriter`).
Copy both out of `FS`; an `.xdmf` without its `.h5` is unreadable. `'XML'`
writes one self-contained file; `'Binary'` writes the `.xdmf` plus one
`<path minus extension><n>.bin` per array. Nothing is on the filesystem until
`finalize()`/`close()` runs.
:::

`'HDF'` works here: the shipped artifact links a wasm32 HDF5 (v8.0.0+). Read a
series back with the ordinary `readMesh` (which resolves the temporal
collection structurally and gives you the first step), `readMeshSelective(path,
{ timeStep: k })` for a particular one (negative counts from the end), or
`readMetadata(path).timeValues` to see the steps without loading any payload.

::: tip Why a handle and not a class
The raw embind surface is an opaque integer handle plus seven free functions
(`xdmfSeriesCreate`, `xdmfSeriesWritePointsCells`, `xdmfSeriesWriteData`,
`xdmfSeriesFinalize`, `xdmfSeriesNumSteps`, `xdmfSeriesFinalized`,
`xdmfSeriesFree`), deliberately **not** an embind `class_`, and the object
above is the wrapper's ergonomic face of it. `@meshioplusplus/wasm` never hands
JS a live C++ object — `NDArray`, `CellBlock` and `Mesh` are all internal — so
an embind class instance, with the Emscripten-specific `.delete()` it comes
with, would be the API's one exception. Free functions also all go through the
same C++ error wrapper, so a `WriteError` arrives as a readable JS `Error`;
a bound *member* function surfaces as a bare, message-less
`WebAssembly.Exception`. It matches the [C API](/c_api)'s `mio_xdmf_series*`
too, so the two flat bindings describe the same object the same way.
:::

This is not a registry format (there is no single `(path, mesh)` call for the
registry to name), so it does **not** change the format counts below and does
not appear in `availableFormats()`.

## Format support

**As of v8.0.0 the WASM build ships every format the C++ core has**, including the five that need HDF5 or netCDF — there is no longer a WASM-specific format gap. That is 41 readable and 43 writable format keys (`openfoam` is read-only; `svg`, `tikz` and `gmsh22` are write-only).

`abaqus`, `ansys`, `ansysInp` (read/write), `avsucd`, `cgns`, `dex`, `dolfin-xml`, `ensight` (EnSight Gold geometry, `.case`/`.geo`, ASCII + C-binary), `exodus`, `flac3d`, `flux`, `freefem`, `gmsh`, `h5m`, `hmf`, `ip`, `mdpa` (Kratos; mesh-level blocks only — see [MDPA](./formats/mdpa.md#c-core)), `med`, `medit`, `mff`, `mfm`, `mphtxt`, `nastran`, `netgen`, `obj`, `off`, `openfoam` (**read-only**, matching the C++/Python core), `permas`, `ply`, `stl`, `su2`, `svg` (**write-only**, 2D visualization), `tecplot`, `tetgen`, `tikz` (**write-only**, 2D LaTeX visualization), `triangle` (`.poly` by default; see the ambiguous-extensions table for `.node`/`.ele`), `ugrid`, `unv`, `vtk`, `vtp`, `vtu` (zlib compression works via Emscripten's built-in port), `wkt`, `xdmf` (XML, Binary **and** HDF). The three field-only formats (`dex`, `ip`, `mff`) read/write geometry-less meshes (field values in `point_data`).

Ask the loaded module rather than trusting this list — it is generated from the same registry the build actually links:

```javascript
const { readers, writers } = m.availableFormats();
```

### What the HDF5/netCDF formats cost, and what changed

- **The `.wasm` is ~5.5 MB, up from ~2.3 MB** (the published npm tarball is ~1.8 MB, since the binary compresses about 3:1 and browsers fetch it compressed). libhdf5 and libnetcdf are statically linked, and they are real bytes. If you do not need these five formats, build your own artifact with `./build/configure-wasm.sh --without-hdf5 --build` (see below) — the JS API is identical and `availableFormats()` reports the smaller set.
- **Breaking: `.xdmf` now writes an HDF companion file.** The registry's XDMF writer default follows the build, exactly as it does natively: with HDF5 present it emits `Format="HDF"` heavy data into a sibling `<base>.h5` and leaves only the XML skeleton in the `.xdmf`. A caller that used to pull one file out of the virtual filesystem must now pull **two**. Reading is unaffected (all three data formats are read).
- **MED cannot write named fields here.** The C++ MED writer defers a mesh carrying `point_data`/`cell_data` to the Python reference writer (the MED-4.1 `CHA` field metadata is inspected byte-for-byte by tests), and this build has no Python to defer to — so it throws by name instead. MED geometry, `point_tags`/`cell_tags` and families are written normally. See [MED quirks](./formats/med.md#quirks-limitations).
- **`sniffFormat` still cannot identify them.** A plain HDF5 magic number says nothing about which of `med`/`cgns`/`h5m`/`hmf`/XDMF-HDF a file is, so the sniffer deliberately never claims it. Use the extension, or pass `format` explicitly.

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

## Known v1 limitations

- **No zero-copy.** Every array is copied once crossing the JS/WASM boundary (see above) — for very large meshes this has a real memory/time cost that the Python bindings' numpy views avoid.
- **No per-format write options.** Parameterized writers (binary vs ASCII, float format strings, gzip levels, VTK 4.2 vs 5.1) use a fixed default matching that format's own Python reference default (e.g. `vtu` writes binary+zlib, `stl` writes ASCII, `gmsh` writes the 4.1 binary format; `stl`/`ply` extract and write the boundary **skin** of a 3D volume mesh — the Python default, see [Skin extraction](./extract_skin.md) — and `svg`/`tikz` render 3D meshes with the default isometric camera). Per-call overrides may be added in a future release. The standalone `extractSkin` utility is not exposed to JS yet (documented follow-up).
- **Named regions are carried** (since v8.1.0). They ride on the mesh object itself — `mesh.regions` is an array of `{ name, kind, dim, tag, entries }` — so `readMesh` / `writeMesh` / `convert` carry them with no extra call, and nothing new had to be forwarded by the wrapper. `kind` is `'point'`, `'cell'` (global block-major cell indices) or `'side'` (`(cell, facet)` pairs). Only the Phase-1 formats (gmsh, abaqus) map onto them so far. See [Named regions](./regions.md).
- **Remaining side-channel data isn't exposed.** `openfoam`'s cell-tag family names, and the `ansysInp`/`unv` set channels pending their Phase-2 region mapping (all carried through a C++ side-channel struct alongside the `Mesh`, mirroring the Python bindings' `AnsysInfo`/`OpenFoamInfo`), are not yet surfaced to JS.
- **No multi-component data arrays.** `point_data`/`cell_data`/`field_data` cross the boundary as flat, shapeless `Float64Array`s (see "The mesh object shape" above) — there is no field carrying a per-array component count, so a 3-component vector or a 3×3 tensor is indistinguishable from N unrelated scalars on the way in. Every array, including through the [data operations](./data_operations.md), is therefore effectively scalar-only in the JS API; `norm(v)`-style expressions and vector/tensor-aware conditioning need the Python, C API, or Fortran bindings, all of which carry the shape. **This is why the path-based calls exist**: `convert` and `convertSurface` never build a JS mesh, so a file passing through them keeps its vector and tensor arrays intact.

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

### The HDF5 and netCDF dependencies

Neither is a system library on this target, and meshio++'s CMake never downloads anything — so `build/build-wasm-deps.sh` produces them, and CMake only *finds* the result:

```sh
./build/build-wasm-deps.sh                 # HDF5 1.14.6 + netcdf-c 4.9.3 -> a prefix
./build/build-wasm-deps.sh --print-prefix   # where that prefix is
```

`configure-wasm.sh` runs it automatically the first time (several minutes) and passes the prefix to CMake via `CMAKE_FIND_ROOT_PATH` — the required knob, since the Emscripten toolchain re-roots every `find_package` at its sysroot, which makes `CMAKE_PREFIX_PATH` and `HINTS` alone ineffective. Pass `--deps-prefix` to reuse a prefix you already have.

Both libraries are built static, `-Oz`, against Emscripten's own zlib port, with everything that assumes an OS this target does not have switched off: HDF5 without threadsafe/MPI/plugins and without the ROS3, direct, mirror and subfiling VFDs; netCDF without DAP, DAP4, byterange, NCZarr, S3, libxml2 (it falls back to its bundled `ezxml`), the `dlopen` plugin loader and mmap. Two upstream-neutral fix-ups are applied and documented in the script: a one-hunk patch guarding HDF5's `feclearexcept(FE_INVALID)` (wasm defines no floating-point exception flags at all), and a rewrite of the exported CMake link interfaces, which name imported targets (`ZLIB::ZLIB`, `HDF5::HDF5`, `hdf5::hdf5_hl`) that no single `find_package` in a consumer defines.

**The wasm stack is sized for them.** `CMakeLists.txt` links the wasm target with `-sSTACK_SIZE=4MB`, well above Emscripten's 64 KiB default, because HDF5's and netCDF-4's frames overrun it. That overrun is silent — the stack grows down into the static data segment — and cost a real bug during development: one Exodus write clobbered libc++'s locale facets, after which every ASCII reader in the module trapped. `-sSTACK_OVERFLOW_CHECK=1` is on so a recurrence aborts loudly instead.

## Selective reads, metadata, and codecs

`readMeshSelective(path, { format, pointsOnly, arrays, timeStep })` and `readMetadata(path, format)`
mirror the Python API; `readerSupportsOptions(format)` reports whether a format has a native
selective path. `arrays: null` reads every data array, `arrays: []` reads none.

```javascript
const mesh = m.readMeshSelective('big.vtu', { arrays: ['u'] });
const meta = m.readMetadata('big.vtu');   // meta.fellBackToFullRead

// A multi-step file: 0 (default) is the first step, negative counts from the end.
const last = m.readMeshSelective('run.exo', { format: 'exodus', timeStep: -1 });
m.readMetadata('run.exo', 'exodus').timeValues;  // [0, 0.5, 1] -- always present
```

**Memory mapping is unavailable** here — the Emscripten virtual filesystem has nothing to map,
so `FileSource` always uses buffered reads. The option is accepted and ignored.

**zstd and lz4 are compiled out.** Unlike HDF5 and netCDF, which this build now carries, neither
has an Emscripten port and neither is worth a from-source dependency for an optional VTK block
codec. zlib (`-sUSE_ZLIB=1`) is unchanged and remains the default codec, so every file the WASM
build wrote before it still round-trips.

## v9.1.0 additions

- `readMeshSelective(path, { lenient: true })` — see
  [`doc/selective_read.md`](selective_read.md).
- XDMF series: `w.flush()`, `w.writeDataArrays(time, pointData, cellData,
  components)`, and `createXdmfTimeSeriesWriter(path, { mode: 'append',
  autoFlush: true })`. With `flush()` the `.xdmf` appears in MEMFS before
  `finalize()`, so a partially-written series can be copied out.

`MdpaInfo` (MDPA properties bodies and entity names) is not exposed, as for every
flat binding.
