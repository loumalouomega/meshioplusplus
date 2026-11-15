// Timings of the JS <-> WebAssembly mesh boundary (roadmap 5.9's probe).
//
// `writeMesh` takes a JS mesh object: every typed array crosses embind into a
// std::vector and then into the core's own buffer, so the cost of ingesting a
// large mesh is the cost of those copies plus the writer. `readMesh` is the
// other direction. This builds a Kuhn tetrahedral cube of about a million
// tetrahedra with one point field and one cell field and times both, best of
// `--repeats`, in milliseconds.
//
// Usage: node tests/wasm/bench.mjs [--cells-a-side 55] [--repeats 3]
//        (after `build/configure-wasm.sh --build` has populated src/wasm/dist)
// MESHIOPLUSPLUS_WASM_VARIANT=seq|mt picks the artifact (default seq: the
// copies are the same in both, and the threaded one only adds start-up).

import { loadMeshioPlusPlus } from '../../src/wasm/src/index.mjs';
import { kuhnCube } from './cube.mjs';

function arg(name, fallback) {
    const i = process.argv.indexOf(name);
    return i >= 0 ? Number(process.argv[i + 1]) : fallback;
}
const n = arg('--cells-a-side', 55);
const repeats = arg('--repeats', 3);
const variant = process.env.MESHIOPLUSPLUS_WASM_VARIANT ?? 'seq';

const m = await loadMeshioPlusPlus({}, { variant });
const mesh = kuhnCube(n);
const nTets = mesh.cells[0].data.length / 4;
const path = '/bench.vtk';

function best(fn) {
    fn(); // warm-up
    let b = Infinity;
    for (let r = 0; r < repeats; r++) {
        const t = performance.now();
        fn();
        b = Math.min(b, performance.now() - t);
    }
    return b;
}

const writeMs = best(() => m.writeMesh(path, mesh, 'vtk', { encoding: 'binary' }));
const readMs = best(() => m.readMesh(path, 'vtk'));
const back = m.readMesh(path, 'vtk');
if (back.cells[0].data.length !== mesh.cells[0].data.length) throw new Error('round trip lost cells');

console.log(`# ${variant}: ${nTets} tetrahedra, ${mesh.points.length / 3} points`);
console.log('step,ms');
console.log(`writeMesh (JS object -> core -> .vtk),${writeMs.toFixed(1)}`);
console.log(`readMesh (.vtk -> core -> JS object),${readMs.toFixed(1)}`);
