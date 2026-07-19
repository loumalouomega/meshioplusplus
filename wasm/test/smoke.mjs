// Smoke test for @meshioplusplus/wasm, run in CI (.github/workflows/wasm.yml)
// after every build and usable as a live usage example. Round-trips a
// synthetic mesh through 3 representative formats (VTU binary+zlib, STL
// binary, OBJ ascii) plus a plain-text read, and exercises writeMesh/
// readMesh/convert/numNodesPerCell through the package's own public API
// (src/index.mjs) -- not the raw embind glue -- so this is exactly what a
// real consumer would call.
//
// Usage: node wasm/test/smoke.mjs   (after `build/configure-wasm.sh --build`
// has populated wasm/dist/meshioplusplus_wasm.{mjs,wasm})

import assert from 'node:assert/strict';
import { loadMeshioPlusPlus } from '../src/index.mjs';

let failed = false;
function step(name, fn) {
    try {
        fn();
        console.log(`ok - ${name}`);
    } catch (err) {
        failed = true;
        console.error(`NOT OK - ${name}`);
        console.error(err);
    }
}

const m = await loadMeshioPlusPlus();

// A small synthetic tetrahedron + a point/cell data field, built directly as
// a JS mesh object (bypassing file I/O) to test the writeMesh(object) path.
const tet = {
    points: new Float64Array([0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1]),
    dim: 3,
    cells: [{ type: 'tetra', data: new Int32Array([0, 1, 2, 3]), nodesPerCell: 4 }],
    point_data: { temperature: new Float64Array([1, 2, 3, 4]) },
    cell_data: { material: [new Float64Array([7])] },
    field_data: {},
};

step('numNodesPerCell metadata table', () => {
    const table = m.numNodesPerCell();
    assert.equal(table.tetra, 4);
    assert.equal(table.triangle, 3);
});

step('topologicalDimension metadata table', () => {
    const table = m.topologicalDimension();
    assert.equal(table.tetra, 3);
    assert.equal(table.triangle, 2);
});

step('VTU binary+zlib round-trip (object -> file -> object)', () => {
    m.writeMesh('/tet.vtu', tet);
    const back = m.readMesh('/tet.vtu');
    assert.equal(back.points.length, 12);
    assert.equal(back.cells.length, 1);
    assert.equal(back.cells[0].type, 'tetra');
    assert.deepEqual(Array.from(back.cells[0].data), [0, 1, 2, 3]);
    assert.deepEqual(Array.from(back.point_data.temperature), [1, 2, 3, 4]);
    assert.deepEqual(Array.from(back.cell_data.material[0]), [7]);
});

step('STL binary round-trip', () => {
    // STL is a surface-only format (one triangle soup, no volume cells) --
    // writing a "tetra" block to it is a legitimate no-op (matches native
    // meshio++: `mesh.write("x.stl")` on a tetra-only mesh silently produces
    // an empty "solid\nendsolid\n", with a warning), so this needs its own
    // triangle-only mesh rather than reusing `tet`.
    const tri = {
        points: new Float64Array([0, 0, 0, 1, 0, 0, 0, 1, 0]),
        dim: 3,
        cells: [{ type: 'triangle', data: new Int32Array([0, 1, 2]), nodesPerCell: 3 }],
    };
    m.writeMesh('/tri.stl', tri, 'stl');
    const back = m.readMesh('/tri.stl', 'stl');
    assert.equal(back.cells.length, 1);
    assert.equal(back.cells[0].type, 'triangle');
    assert.deepEqual(Array.from(back.cells[0].data), [0, 1, 2]);
});

step('OBJ ascii read from a hand-written file', () => {
    const obj = 'v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n';
    m.FS.writeFile('/tri.obj', obj);
    const mesh = m.readMesh('/tri.obj');
    assert.equal(mesh.points.length, 9);
    assert.equal(mesh.cells.length, 1);
    assert.equal(mesh.cells[0].type, 'triangle');
    assert.deepEqual(Array.from(mesh.cells[0].data), [0, 1, 2]);
});

step('convert() reads one format and writes another directly', () => {
    m.convert('/tri.obj', '/tri.vtk');
    const back = m.readMesh('/tri.vtk');
    assert.equal(back.cells[0].type, 'triangle');
});

step('format collision: .msh defaults to gmsh, explicit override selects ansys', () => {
    m.writeMesh('/tet.msh', tet, 'ansys');
    let threw = false;
    try {
        m.readMesh('/tet.msh', 'gmsh'); // an ansys file is not valid gmsh
    } catch (err) {
        threw = true;
        assert.ok(err instanceof Error);
        assert.ok(err.message.length > 0, 'error should carry a real message');
    }
    assert.ok(threw, 'expected reading an ansys-written .msh as gmsh to throw');
});

step('unknown format raises a catchable Error, not an abort', () => {
    assert.throws(() => m.readMesh('/does/not/exist.obj'), /Could not open file/);
});

// A triangle mesh reused across the extra format round-trips below.
const tri2 = {
    points: new Float64Array([0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0]),
    dim: 3,
    cells: [{ type: 'triangle', data: new Int32Array([0, 1, 2, 0, 2, 3]), nodesPerCell: 3 }],
};

step('PLY binary round-trip', () => {
    m.writeMesh('/tri.ply', tri2, 'ply');
    const back = m.readMesh('/tri.ply', 'ply');
    assert.equal(back.cells.length, 1);
    assert.equal(back.cells[0].type, 'triangle');
    assert.deepEqual(Array.from(back.cells[0].data), [0, 1, 2, 0, 2, 3]);
});

step('OFF ascii round-trip', () => {
    m.writeMesh('/tri.off', tri2, 'off');
    const back = m.readMesh('/tri.off', 'off');
    assert.equal(back.cells[0].type, 'triangle');
    assert.deepEqual(Array.from(back.cells[0].data), [0, 1, 2, 0, 2, 3]);
});

step('GMSH ascii round-trip (tetra, volume format)', () => {
    m.writeMesh('/tet.gmsh.msh', tet, 'gmsh');
    const back = m.readMesh('/tet.gmsh.msh', 'gmsh');
    assert.equal(back.cells[0].type, 'tetra');
    assert.deepEqual(Array.from(back.cells[0].data), [0, 1, 2, 3]);
    assert.equal(back.points.length, 12);
});

step('malformed file raises a catchable Error, not a WASM abort', () => {
    m.FS.writeFile('/bad.vtu', '<?xml version="1.0"?><NotVTK></NotVTK>');
    assert.throws(
        () => m.readMesh('/bad.vtu', 'vtu'),
        (err) => err instanceof Error && err.message.length > 0,
    );
});

if (failed) {
    console.error('\nSMOKE TEST FAILED');
    process.exit(1);
}
console.log('\nSMOKE TEST PASSED');
