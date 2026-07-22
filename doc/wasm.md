# WebAssembly / JavaScript

The C++ core also compiles to WebAssembly and ships as an npm package, [`@meshioplusplus/wasm`](https://www.npmjs.com/package/@meshioplusplus/wasm), for reading and writing meshes in the browser or Node.js. (It is one of two "flat" bindings over the same core and shared format-dispatch registry — the other is the [C API](/c_api), which native HDF5/netCDF-capable builds can use.)

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

**Ragged cell blocks** (polygon/polyhedron blocks with a varying node count per cell, e.g. MED Voronoi polygons or OpenFOAM general polyhedra) are not representable in this flat shape and are rejected: `readMesh` throws if the file contains one, and there is no way to construct one for `writeMesh`.

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

## Format support

The WASM build ships the 36 formats with no HDF5/netCDF dependency, plus XDMF's XML/Binary data path (not its HDF variant) — 35 readable formats in total, 36 writable (`openfoam` is read-only; `svg` and `tikz` are write-only):

`abaqus`, `ansys`, `ansysInp` (read/write), `avsucd`, `dex`, `dolfin-xml`, `ensight` (EnSight Gold geometry, `.case`/`.geo`, ASCII + C-binary), `flac3d`, `flux`, `freefem`, `gmsh`, `ip`, `medit`, `mff`, `mfm`, `mphtxt`, `nastran`, `netgen`, `obj`, `off`, `openfoam` (**read-only**, matching the C++/Python core), `permas`, `ply`, `stl`, `su2`, `svg` (**write-only**, 2D visualization), `tecplot`, `tetgen`, `tikz` (**write-only**, 2D LaTeX visualization), `triangle` (`.poly` by default; see the ambiguous-extensions table for `.node`/`.ele`), `ugrid`, `unv`, `vtk`, `vtp`, `vtu` (zlib compression works via Emscripten's built-in port), `wkt`, `xdmf` (XML/Binary only). The three field-only formats (`dex`, `ip`, `mff`) read/write geometry-less meshes (field values in `point_data`).

**Not yet supported: `cgns`, `h5m`, `hmf`, `med`, `exodus`.** All five need HDF5 and/or netCDF, which are not built for this target — porting those C libraries to WebAssembly is a separate, materially larger undertaking than the rest of the C++ core (both have autotools/CMake builds assuming a POSIX filesystem and, in HDF5's case, sometimes MPI). They may follow in a future release; there is no runtime fallback the way there is for the Python bindings, since there's no Python present at all in this build.

### Ambiguous extensions

Some extensions are shared by more than one format. `readMesh`/`writeMesh`/ `convert` all take an optional trailing `format` argument (or an `{inFormat, outFormat}` options object for `convert`) to disambiguate, mirroring Python's `file_format=` kwarg:

| Extension | Default format | Pass `format=` to select instead |
|-----------|-----------------|-----------------------------------|
| `.msh` | `gmsh` | `"ansys"`, `"freefem"` |
| `.inp` | `abaqus` | `"ansysinp"` |
| `.node` / `.ele` | `tetgen` | `"triangle"` (2D Triangle pairs) |

## Known v1 limitations

- **No zero-copy.** Every array is copied once crossing the JS/WASM boundary (see above) — for very large meshes this has a real memory/time cost that the Python bindings' numpy views avoid.
- **No per-format write options.** Parameterized writers (binary vs ASCII, float format strings, gzip levels, VTK 4.2 vs 5.1) use a fixed default matching that format's own Python reference default (e.g. `vtu` writes binary+zlib, `stl` writes ASCII, `gmsh` writes the 4.1 binary format; `stl`/`ply` extract and write the boundary **skin** of a 3D volume mesh — the Python default, see [Skin extraction](./extract_skin.md) — and `svg`/`tikz` render 3D meshes with the default isometric camera). Per-call overrides may be added in a future release. The standalone `extractSkin` utility is not exposed to JS yet (documented follow-up).
- **Side-channel data isn't exposed.** `ansysInp`'s `point_sets`/`cell_sets` and `openfoam`'s cell-tag family names (both carried through a C++ side-channel struct alongside the `Mesh`, mirroring the Python bindings' `AnsysInfo`/`OpenFoamInfo`) are not yet surfaced to JS — reading/writing the mesh geometry and data itself works, but these format-specific extras are dropped for now.
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

`build/configure-wasm.sh` always configures with `-DMESHIOPLUSPLUS_BUILD_PYTHON=OFF` (no Python/pybind11 involved), `-DMESHIOPLUSPLUS_PARALLEL_BACKEND=SEQ` (OpenMP/TBB/the parallel STL have no meaningful story on this target yet), `-DMESHIOPLUSPLUS_MESH_BACKEND=NATIVE` (the fastest [in-memory mesh backend](cpp_backends.md) — canonical Float64/Int64 storage, so the embind typed-array boundary needs no dtype dispatch; the JS API shape is unchanged, and `meshBackend()` on the loaded module reports `"native"`), and HDF5/netCDF off. See `--help` for the `--without-zlib`/`--build-type` options. CI (`.github/workflows/wasm.yml`) builds and smoke-tests on every push/PR and publishes to npm on `v*` tags.

## Selective reads, metadata, and codecs

`readMeshSelective(path, { format, pointsOnly, arrays })` and `readMetadata(path, format)`
mirror the Python API; `readerSupportsOptions(format)` reports whether a format has a native
selective path. `arrays: null` reads every data array, `arrays: []` reads none.

```javascript
const mesh = m.readMeshSelective('big.vtu', { arrays: ['u'] });
const meta = m.readMetadata('big.vtu');   // meta.fellBackToFullRead
```

**Memory mapping is unavailable** here — the Emscripten virtual filesystem has nothing to map,
so `FileSource` always uses buffered reads. The option is accepted and ignored.

**zstd and lz4 are compiled out**, consistent with the HDF5/netCDF-backed formats: there is no
Emscripten port for either. zlib (`-sUSE_ZLIB=1`) is unchanged and remains the default codec,
so every file the WASM build wrote before it still round-trips.
