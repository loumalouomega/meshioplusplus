# @meshioplusplus/wasm

[meshio++](https://github.com/loumalouomega/meshioplusplus) mesh I/O, compiled to
WebAssembly. Read 80 and write 75 mesh file formats (VTK, VTU, Gmsh, STL, OBJ,
Nastran, MED, CGNS, Exodus, and more) in the browser or Node.js — the same set
the native library supports, with libhdf5 and libnetcdf compiled in. The full
build adds STEP and IGES reading through OpenCASCADE (82 readable).

Full docs, the format-support table, and known v1 limitations:
[doc/wasm.md](https://loumalouomega.github.io/meshioplusplus/wasm).

## Install

```sh
npm install @meshioplusplus/wasm
```

## Usage

```js
import { loadMeshioPlusPlus } from "@meshioplusplus/wasm";

const meshio = await loadMeshioPlusPlus();

// Write bytes into the virtual filesystem, then read them as a mesh.
const response = await fetch("example.vtu");
meshio.FS.writeFile("/example.vtu", new Uint8Array(await response.arrayBuffer()));
const mesh = meshio.readMesh("/example.vtu");

console.log(mesh.points);          // Float64Array, flat (numPoints * dim)
console.log(mesh.cells[0].type);   // e.g. "triangle"
console.log(mesh.cells[0].data);   // Int32Array connectivity

// Convert directly, or round-trip through a JS mesh object.
meshio.convert("/example.vtu", "/example.stl");
meshio.writeMesh("/example.msh", mesh, "gmsh"); // .msh needs an explicit
                                                 // format when writing
                                                 // ansys/freefem instead
```

See [doc/wasm.md](https://loumalouomega.github.io/meshioplusplus/wasm) for the mesh
object shape, the full list of supported formats, and format-selection rules
for ambiguous extensions (`.msh`, `.inp`).

## The full build: STEP, IGES, fitCad, zstd, lz4, bzip2, KaHIP

The default artifact is lean and MIT. A second, larger artifact links
OpenCASCADE (reading STEP and IGES, and `fitCad`, which turns a fitted skin
into an analytic CAD solid), the zstd and lz4 VTK XML codecs, bzip2, and the
serial KaHIP kaffpa partitioner (`partition` with `method: 'kahip'`). Load
it with `{ full: true }`; the same `variant` choice applies, and only the
artifact you load is fetched:

```js
const meshio = await loadMeshioPlusPlus({}, { full: true });
meshio.hasOcct();                          // true
const part = meshio.readMeshSelective("/part.step", { info: true });
part.info.faces.length;                    // the B-rep, as CAD JSON
```

**Licence:** OpenCASCADE is LGPL-2.1 with the OCCT exception, so the full
artifact is not plain MIT; see
[THIRD_PARTY_NOTICES.md](./THIRD_PARTY_NOTICES.md) for every bundled component
and how the artifact is rebuilt from source. An application that never passes
`full: true` loads only MIT and permissively licensed code.

A bundler never sees the full artifact's import (it would otherwise ship it in
every lean bundle): copy `dist/meshioplusplus_wasm_full*.{mjs,wasm}` yourself
and pass `{ full: true, fullGlue: url }` with a matching `locateFile`.
