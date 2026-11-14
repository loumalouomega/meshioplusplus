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

function arg(name, fallback) {
    const i = process.argv.indexOf(name);
    return i >= 0 ? Number(process.argv[i + 1]) : fallback;
}
const n = arg('--cells-a-side', 55);
const repeats = arg('--repeats', 3);
const variant = process.env.MESHIOPLUSPLUS_WASM_VARIANT ?? 'seq';

// The Kuhn tetrahedralisation of the unit cube: 6 * n^3 tetrahedra.
function kuhnCube(n) {
    const np = n + 1;
    const points = new Float64Array(np * np * np * 3);
    for (let k = 0, p = 0; k < np; k++)
        for (let j = 0; j < np; j++)
            for (let i = 0; i < np; i++, p += 3) {
                points[p] = i / n;
                points[p + 1] = j / n;
                points[p + 2] = k / n;
            }
    const id = (i, j, k) => (k * np + j) * np + i;
    const paths = [
        [[1, 0, 0], [1, 1, 0]], [[1, 0, 0], [1, 0, 1]], [[0, 1, 0], [1, 1, 0]],
        [[0, 1, 0], [0, 1, 1]], [[0, 0, 1], [1, 0, 1]], [[0, 0, 1], [0, 1, 1]],
    ];
    const conn = new Int32Array(6 * n * n * n * 4);
    let c = 0;
    for (let k = 0; k < n; k++)
        for (let j = 0; j < n; j++)
            for (let i = 0; i < n; i++)
                for (const [a, b] of paths) {
                    const v = [[0, 0, 0], a, b, [1, 1, 1]].map(([x, y, z]) => id(i + x, j + y, k + z));
                    // Orient every tetrahedron positively, so the file is a valid mesh.
                    const p = v.map((q) => points.subarray(q * 3, q * 3 + 3));
                    const d = [1, 2, 3].map((m) => [0, 1, 2].map((x) => p[m][x] - p[0][x]));
                    const det =
                        d[0][0] * (d[1][1] * d[2][2] - d[1][2] * d[2][1]) -
                        d[0][1] * (d[1][0] * d[2][2] - d[1][2] * d[2][0]) +
                        d[0][2] * (d[1][0] * d[2][1] - d[1][1] * d[2][0]);
                    if (det < 0) [v[2], v[3]] = [v[3], v[2]];
                    conn.set(v, c);
                    c += 4;
                }
    const nTet = c / 4;
    const u = new Float64Array(np * np * np);
    for (let q = 0; q < u.length; q++) u[q] = points[q * 3] ** 2 + points[q * 3 + 1] * points[q * 3 + 2];
    const material = new Float64Array(nTet).map((_, q) => q % 7);
    return {
        points,
        dim: 3,
        cells: [{ type: 'tetra', data: conn, nodesPerCell: 4 }],
        point_data: { u },
        cell_data: { material: [material] },
        field_data: {},
    };
}

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
