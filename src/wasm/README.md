# @meshioplusplus/wasm

[meshio++](https://github.com/loumalouomega/meshioplusplus) mesh I/O, compiled to
WebAssembly. Read and write 41 mesh file formats (VTK, VTU, Gmsh, STL, OBJ,
Nastran, MED, CGNS, Exodus, and more) in the browser or Node.js — the same set
the native library supports, with libhdf5 and libnetcdf compiled in.

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
