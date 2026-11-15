// The WebAssembly heap grows past 2 GB (roadmap 5.13.1).
//
// Emscripten's default ceiling is 2 GB, and the core ran out of it at about 16
// million tetrahedra: `readMesh` threw `std::bad_alloc` (a Kuhn cube of 140
// cells a side, 16.5 million, failed; 130, 13.2 million, passed). The artifacts
// are linked with -sMAXIMUM_MEMORY=4GB, so a mesh of 20 million tetrahedra, whose
// core copy, in-memory file and read-back need well over 2 GB of heap, must
// round-trip. This is that mesh. It takes about 15 seconds and 4 GB of memory,
// so it is not part of smoke.mjs.
//
// Usage: node tests/wasm/memory.mjs [--cells-a-side 150]
//        MESHIOPLUSPLUS_WASM_VARIANT=seq|mt (default seq)

import assert from 'node:assert/strict';
import { loadMeshioPlusPlus } from '../../src/wasm/src/index.mjs';
import { kuhnCube } from './cube.mjs';

const i = process.argv.indexOf('--cells-a-side');
const n = i >= 0 ? Number(process.argv[i + 1]) : 150;
const variant = process.env.MESHIOPLUSPLUS_WASM_VARIANT ?? 'seq';

const m = await loadMeshioPlusPlus({}, { variant });
const mesh = kuhnCube(n);
const tets = mesh.cells[0].data.length / 4;
console.log(`${variant}: ${tets} tetrahedra`);
m.writeMesh('/memory.vtk', mesh, 'vtk', { encoding: 'binary' });
const back = m.readMesh('/memory.vtk', 'vtk');
assert.equal(back.cells[0].data.length, mesh.cells[0].data.length);
assert.equal(back.points.length, mesh.points.length);
console.log('ok - a mesh needing more than 2 GB of heap round-trips');
