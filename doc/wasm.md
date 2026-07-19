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
```

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
`signedVolume`, `unsignedVolume`, `numInverted`). See
[transform](./transform.md), [clean](./clean.md), [crop](./crop.md),
[split](./split.md), and [stats](./stats.md).

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

## Building from source

Requires the [Emscripten SDK](https://emscripten.org/docs/getting_started/downloads.html):

```sh
git clone https://github.com/emscripten-core/emsdk.git
cd emsdk && ./emsdk install latest && ./emsdk activate latest
source ./emsdk_env.sh
cd ../meshioplusplus  # this repo
./build/configure-wasm.sh --build
node wasm/test/smoke.mjs
```

`build/configure-wasm.sh` always configures with `-DMESHIOPLUSPLUS_BUILD_PYTHON=OFF` (no Python/pybind11 involved), `-DMESHIOPLUSPLUS_PARALLEL_BACKEND=SEQ` (OpenMP/TBB/the parallel STL have no meaningful story on this target yet), `-DMESHIOPLUSPLUS_MESH_BACKEND=NATIVE` (the fastest [in-memory mesh backend](cpp_backends.md) — canonical Float64/Int64 storage, so the embind typed-array boundary needs no dtype dispatch; the JS API shape is unchanged, and `meshBackend()` on the loaded module reports `"native"`), and HDF5/netCDF off. See `--help` for the `--without-zlib`/`--build-type` options. CI (`.github/workflows/wasm.yml`) builds and smoke-tests on every push/PR and publishes to npm on `v*` tags.
