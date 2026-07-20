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

// --- data operations -------------------------------------------------------
// These act on the mesh's data arrays; the geometry must come through
// untouched. `tetv` carries point_data.temperature and cell_data.material.
//
// The WASM mesh-object shape crosses point_data/cell_data/field_data as flat,
// shapeless Float64Arrays (see js_bindings.cpp's mesh_to_val/val_to_mesh) --
// there is no way to convey a per-array component count through this API, so
// every array here is a plain scalar. A multi-component (vector/tensor) field
// is a pre-existing, out-of-scope limitation of the whole WASM binding, not
// something the data operations can work around.

const tetv = {
    points: new Float64Array([0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1]),
    dim: 3,
    cells: [{ type: 'tetra', data: new Int32Array([0, 1, 2, 3]), nodesPerCell: 4 }],
    point_data: {
        temperature: new Float64Array([1, 2, 3, 4]),
        pressure: new Float64Array([10, 20, 30, 40]),
    },
    cell_data: { material: [new Float64Array([7])] },
    field_data: {},
};

step('dataInfo summarizes every array', () => {
    const arrays = m.dataInfo(tetv);
    assert.equal(arrays.length, 3);
    const t = arrays.find((a) => a.name === 'temperature');
    assert.equal(t.location, 'point_data');
    assert.equal(t.numEntries, 4);
    assert.equal(t.numComponents, 1);
    assert.equal(t.min, 1);
    assert.equal(t.max, 4);
    assert.equal(t.mean, 2.5);
    assert.equal(t.numNan, 0);
    const p = arrays.find((a) => a.name === 'pressure');
    assert.equal(p.min, 10);
    assert.equal(p.max, 40);
});

step('dataCalc evaluates an expression', () => {
    const out = m.dataCalc(tetv, '2 * temperature + 1', 'point', 'derived', false);
    const derived = out.point_data.derived;
    assert.ok(Math.abs(derived[0] - 3) < 1e-12);
    assert.ok(Math.abs(derived[3] - 9) < 1e-12);
    // Geometry is never modified.
    assert.equal(out.points.length, 12);
    assert.equal(out.cells.length, 1);
});

step('dataCalc rejects an unknown function with a catchable Error', () => {
    assert.throws(
        () => m.dataCalc(tetv, 'log(temperature)', 'point', 'bad', false),
        (err) => err instanceof Error && err.message.length > 0,
    );
});

step('dataDrop / dataKeep / dataRename', () => {
    const dropped = m.dataDrop(tetv, 'point', ['temperature'], false);
    assert.ok(!('temperature' in dropped.point_data));
    assert.ok('pressure' in dropped.point_data);

    const kept = m.dataKeep(tetv, 'point', ['temperature'], false);
    assert.deepEqual(Object.keys(kept.point_data), ['temperature']);

    const renamed = m.dataRename(tetv, 'point', 'temperature', 'T');
    assert.ok('T' in renamed.point_data);
    assert.ok(!('temperature' in renamed.point_data));
});

step('dataPointToCell / dataCellToPoint', () => {
    const toCell = m.dataPointToCell(tetv, ['temperature'], '_c');
    // Single tetra: the mean of {1,2,3,4}.
    assert.ok(Math.abs(toCell.cell_data.temperature_c[0][0] - 2.5) < 1e-12);

    const toPoint = m.dataCellToPoint(tetv, ['material'], 'uniform', '');
    assert.equal(toPoint.point_data.material.length, 4);
    assert.ok(Math.abs(toPoint.point_data.material[0] - 7) < 1e-12);
});

step('dataCondition normalizes to [0, 1]', () => {
    const out = m.dataCondition(
        tetv, 'point', ['temperature'], 'normalize', 0, 1,
        'component', 'ignore', 0, '',
    );
    const t = out.point_data.temperature;
    assert.ok(Math.abs(t[0] - 0) < 1e-12);
    assert.ok(Math.abs(t[3] - 1) < 1e-12);
});


// ---------------------------------------------------------------------------
// Selective reads, file summaries, and the codec build profile
// ---------------------------------------------------------------------------

const seltri = {
    points: new Float64Array([0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0]),
    dim: 3,
    cells: [{ type: 'triangle', data: new Int32Array([0, 1, 2, 0, 2, 3]), nodesPerCell: 3 }],
    point_data: {
        u: new Float64Array([1, 2, 3, 4]),
        v: new Float64Array([10, 20, 30, 40]),
    },
};
m.writeMesh('/selective.vtu', seltri);

step('readMeshSelective: no options reads every array', () => {
    const mesh = m.readMeshSelective('/selective.vtu');
    assert.deepEqual(Object.keys(mesh.point_data).sort(), ['u', 'v']);
});

step('readMeshSelective: pointsOnly keeps geometry, drops data', () => {
    const mesh = m.readMeshSelective('/selective.vtu', { pointsOnly: true });
    assert.equal(Object.keys(mesh.point_data).length, 0);
    // pointsOnly narrows data, not topology.
    assert.equal(mesh.points.length, 12);
    assert.equal(mesh.cells[0].data.length, 6);
});

step('readMeshSelective: arrays subset, and [] vs null', () => {
    assert.deepEqual(
        Object.keys(m.readMeshSelective('/selective.vtu', { arrays: ['u'] }).point_data), ['u']);
    // The distinction that motivates std::optional all the way down: an empty
    // list means *no* arrays, null means *every* array.
    assert.equal(
        Object.keys(m.readMeshSelective('/selective.vtu', { arrays: [] }).point_data).length, 0);
    assert.deepEqual(
        Object.keys(m.readMeshSelective('/selective.vtu', { arrays: null }).point_data).sort(),
        ['u', 'v']);
});

step('readMetadata summarizes without loading the arrays', () => {
    const meta = m.readMetadata('/selective.vtu');
    assert.equal(meta.numPoints, 4);
    assert.equal(meta.numCells, 2);
    assert.equal(meta.cellBlocks.length, 1);
    assert.equal(meta.cellBlocks[0].type, 'triangle');
    assert.equal(meta.cellBlocks[0].nodesPerCell, 3);
    assert.deepEqual(meta.pointDataNames.sort(), ['u', 'v']);
    assert.equal(meta.format, 'vtu');
    // vtu has a native metadata path, so this really was cheap...
    assert.equal(meta.fellBackToFullRead, false);
    // ...and a native summary never decodes the coordinates, so it reports no
    // bbox rather than a fabricated one at the origin.
    assert.ok(!('bboxMin' in meta));
});

step('readMetadata flags a full-read fallback and can then afford a bbox', () => {
    m.writeMesh('/selective.stl', seltri, 'stl');
    const meta = m.readMetadata('/selective.stl');
    assert.equal(meta.fellBackToFullRead, true);
    assert.ok('bboxMin' in meta && 'bboxMax' in meta);
});

step('readerSupportsOptions reports the native paths', () => {
    assert.equal(m.readerSupportsOptions('vtu'), true);
    assert.equal(m.readerSupportsOptions('stl'), false);
});

step('zstd/lz4 are compiled out of the WASM build', () => {
    // Consistent with the HDF5/netCDF-backed formats: no Emscripten port
    // exists, so zlib remains the only block codec here.
    m.writeMesh('/codec.vtu', seltri);
    const text = new TextDecoder().decode(m.FS.readFile('/codec.vtu'));
    assert.ok(!text.includes('vtkZSTDDataCompressor'));
    assert.ok(!text.includes('vtkLZ4DataCompressor'));
});

// The geometry operations below are reached through the package wrapper
// (src/index.mjs), not `Module.*` directly. That is the whole point: an embind
// binding the wrapper does not forward is unreachable from
// loadMeshioPlusPlus(), which is exactly how these ops shipped broken before.

const cube = {
    points: new Float64Array([
        0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1,
    ]),
    dim: 3,
    cells: [
        { type: 'hexahedron', data: new Int32Array([0, 1, 2, 3, 4, 5, 6, 7]), nodesPerCell: 8 },
    ],
    point_data: {},
    cell_data: {},
    field_data: {},
};

step('convertCells simplexify: one hexahedron -> 6 tetra', () => {
    const out = m.convertCells(cube, 'simplexify');
    assert.equal(out.cells.length, 1);
    assert.equal(out.cells[0].type, 'tetra');
    assert.equal(out.cells[0].data.length, 6 * 4);
    // The decomposition reuses the parent's own corner nodes.
    assert.equal(out.points.length, 24);
});

step('convertCells elevate/linearize round-trips', () => {
    const up = m.convertCells(tet, 'elevate');
    assert.equal(up.cells[0].type, 'tetra10');
    assert.equal(up.points.length, (4 + 6) * 3);
    const down = m.convertCells(up, 'linearize');
    assert.equal(down.cells[0].type, 'tetra');
    assert.equal(down.points.length, 12);
});

step('convertCells rejects a full-Lagrange elevate target', () => {
    const quad9 = {
        points: new Float64Array(9 * 3),
        dim: 3,
        cells: [
            {
                type: 'quad9',
                data: new Int32Array([0, 1, 2, 3, 4, 5, 6, 7, 8]),
                nodesPerCell: 9,
            },
        ],
        point_data: {},
        cell_data: {},
        field_data: {},
    };
    assert.throws(() => m.convertCells(quad9, 'elevate'));
});

step('geometry operations are reachable through the wrapper', () => {
    // Regression guard for the v7.2.1 class of bug: every geometry op must be
    // forwarded by src/index.mjs, not merely bound in js_bindings.cpp.
    for (const name of [
        'extractSurface',
        'extractSkin',
        'attachQuality',
        'sniffFormat',
        'reorder',
        'computeBandwidth',
        'diff',
        'meshesEqual',
        'merge',
        'transform',
        'clean',
        'cropBbox',
        'cropPlane',
        'split',
        'convertCells',
        'stats',
        'meshBackend',
    ]) {
        assert.equal(typeof m[name], 'function', `${name} is not forwarded by the wrapper`);
    }
});

step('a forwarded geometry operation actually runs', () => {
    const skin = m.extractSurface(cube);
    assert.equal(skin.cells[0].type, 'quad');
    assert.equal(skin.cells[0].data.length, 6 * 4);

    const s = m.stats(cube);
    assert.equal(s.numPoints, 8);
    assert.equal(s.numInverted, 0);

    const pieces = m.split(cube, 'type');
    assert.equal(pieces.length, 1);
    assert.equal(pieces[0].key, 'hexahedron');

    assert.equal(m.meshesEqual(cube, cube), true);
    assert.equal(typeof m.meshBackend(), 'string');
});

if (failed) {
    console.error('\nSMOKE TEST FAILED');
    process.exit(1);
}
console.log('\nSMOKE TEST PASSED');
