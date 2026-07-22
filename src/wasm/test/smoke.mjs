// Smoke test for @meshioplusplus/wasm, run in CI (.github/workflows/wasm.yml)
// after every build and usable as a live usage example. Round-trips a
// synthetic mesh through 3 representative formats (VTU binary+zlib, STL
// binary, OBJ ascii) plus a plain-text read, and exercises writeMesh/
// readMesh/convert/numNodesPerCell through the package's own public API
// (src/index.mjs) -- not the raw embind glue -- so this is exactly what a
// real consumer would call.
//
// Usage: node src/wasm/test/smoke.mjs   (after `build/configure-wasm.sh --build`
// has populated src/wasm/dist/meshioplusplus_wasm.{mjs,wasm})

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

step('refine: one hexahedron -> 8 hexahedra with 27 nodes', () => {
    const out = m.refine(cube);
    assert.equal(out.cells.length, 1);
    assert.equal(out.cells[0].type, 'hexahedron');
    assert.equal(out.cells[0].data.length, 8 * 8);
    // 8 corners + 12 edge mids + 6 face centres + 1 body centre.
    assert.equal(out.points.length, 27 * 3);
});

step('refine levels=2 matches refining twice', () => {
    const direct = m.refine(tet, 2);
    const twice = m.refine(m.refine(tet));
    assert.equal(direct.points.length, twice.points.length);
    assert.deepEqual(Array.from(direct.cells[0].data), Array.from(twice.cells[0].data));
});

step('refine rejects a cell type with no same-type subdivision', () => {
    const up = m.convertCells(tet, 'elevate');
    assert.throws(() => m.refine(up));
});

step('decimate: collapses a refined cube skin, pinning its creases', () => {
    // The skin of a refined cube: 24 quads -> 48 triangles, with every cube
    // edge/corner vertex a pinned feature; only face-interior vertices go.
    const skin = m.extractSkin(m.refine(cube), true);
    const out = m.decimate(skin, 0.5);
    assert.equal(out.mesh.cells.length, 1);
    assert.equal(out.mesh.cells[0].type, 'triangle');
    assert.ok(out.mesh.cells[0].data.length / 3 < 48);
    assert.ok(out.facesRemoved > 0);
    assert.ok(out.pointsRemoved > 0);
    assert.ok(out.collapsesRejected >= 0);
    assert.ok(out.maxErrorApplied >= 0);
});

step('decimate rejects a volume mesh and a missing criterion', () => {
    assert.throws(() => m.decimate(cube, 0.5), /extract_surface/);
    const skin = m.extractSkin(cube, true);
    assert.throws(() => m.decimate(skin));
});

step('decimate is reachable as a convertSurfaceOps pipeline step', () => {
    // The pipeline runs against the loaded mesh itself, so hand it a surface
    // mesh (decimate refuses volume input by design).
    m.writeMesh('/dec.vtu', m.extractSkin(m.refine(cube), true));
    const out = m.convertSurfaceOps('/dec.vtu', '/dec.vtp', [{ op: 'decimate' }]);
    assert.equal(out.steps[0].op, 'decimate');
    assert.equal(typeof out.steps[0].facesRemoved, 'number');
    assert.ok(m.readMesh('/dec.vtp').cells[0].data.length / 3 < 48);
});

step('smooth: relaxes an interior node while pinning the boundary', () => {
    // A refined cube is a 3x3x3 node lattice: exactly one node (the body
    // centre) is interior, so it is the only one smoothing may move.
    const grid = m.refine(cube);
    const points = Float64Array.from(grid.points);
    const inside = (v) => v > 1e-9 && v < 1 - 1e-9;
    let interior = -1;
    for (let i = 0; i < points.length / 3; ++i) {
        if (inside(points[3 * i]) && inside(points[3 * i + 1]) && inside(points[3 * i + 2])) {
            interior = i;
            break;
        }
    }
    assert.notEqual(interior, -1, 'the refined cube should have an interior node');

    // Pull it off its own centroid so there is something to relax.
    points[3 * interior] += 0.2;
    const perturbed = { ...grid, points };

    const out = m.smooth(perturbed, 'laplacian', 5);
    // Geometry only: same point count, same connectivity, moved coordinates.
    assert.equal(out.mesh.points.length, perturbed.points.length);
    assert.equal(out.mesh.cells.length, perturbed.cells.length);
    assert.deepEqual(
        Array.from(out.mesh.cells[0].data),
        Array.from(perturbed.cells[0].data),
    );
    assert.notEqual(out.mesh.points[3 * interior], perturbed.points[3 * interior]);
    // The interior node is pulled back toward the centre it was moved from.
    assert.ok(
        Math.abs(out.mesh.points[3 * interior] - 0.5) <
            Math.abs(perturbed.points[3 * interior] - 0.5),
    );
    // Every boundary node is pinned, so only that one node moved.
    for (let i = 0; i < points.length / 3; ++i) {
        if (i === interior) continue;
        for (let c = 0; c < 3; ++c)
            assert.equal(out.mesh.points[3 * i + c], perturbed.points[3 * i + c]);
    }

    assert.equal(out.numNodesMoved, 1);
    assert.ok(out.maxDisplacement > 0);
    assert.equal(typeof out.numSkippedInversion, 'number');
});

step('smooth: taubin defaults leave a structured hex block alone', () => {
    // The negative default lambda means "the method's own default" and must
    // reach the core unchanged; an already-relaxed lattice is a fixed point.
    const grid = m.refine(cube);
    const out = m.smooth(grid);
    assert.equal(out.numNodesMoved, 0);
    assert.deepEqual(Array.from(out.mesh.points), Array.from(grid.points));
});

step('smooth rejects an unknown method', () => {
    assert.throws(() => m.smooth(cube, 'not-a-method'));
});

step('interpolate: transfers fields onto a target mesh', () => {
    // Nearest: the single target point sits next to the second source point.
    const src = {
        points: new Float64Array([0, 0, 0, 1, 0, 0]),
        dim: 3,
        cells: [],
        point_data: { f: new Float64Array([10, 20]) },
        cell_data: {},
        field_data: {},
    };
    const tgt = {
        points: new Float64Array([0.9, 0, 0]),
        dim: 3,
        cells: [],
        point_data: {},
        cell_data: {},
        field_data: {},
    };
    const out = m.interpolate(src, tgt);
    assert.deepEqual(Array.from(out.point_data.f), [20]);

    // Barycentric is exact on a linear field (g = x at the cube's corners).
    const cubeSrc = { ...cube, point_data: { g: new Float64Array([0, 1, 1, 0, 0, 1, 1, 0]) } };
    const probe = {
        points: new Float64Array([0.25, 0.5, 0.5]),
        dim: 3,
        cells: [],
        point_data: {},
        cell_data: {},
        field_data: {},
    };
    const bar = m.interpolate(cubeSrc, probe, 'barycentric');
    assert.ok(Math.abs(bar.point_data.g[0] - 0.25) < 1e-12);

    assert.throws(() => m.interpolate(src, tgt, 'not-a-method'));
    assert.throws(() => m.interpolate(src, tgt, 'nearest', ['nope']));
});

step('slice: the mid-plane cross-section of the unit cube has area 1', () => {
    const sectioned = m.slice(cube, [0, 0, 0.5], [0, 0, 1], true);
    // A volume mesh sections into surface faces (triangles here).
    let area = 0;
    let nfaces = 0;
    for (const block of sectioned.cells) {
        const npc = block.nodesPerCell;
        const data = block.data;
        for (let c = 0; c < data.length / npc; ++c) {
            nfaces += 1;
            const p = [];
            for (let k = 0; k < npc; ++k) {
                const n = data[c * npc + k];
                p.push([sectioned.points[n * 3], sectioned.points[n * 3 + 1]]);
            }
            // Shoelace on the projected (x, y) ring — the section lies at z = 0.5.
            let s = 0;
            for (let k = 0; k < p.length; ++k) {
                const a = p[k];
                const b = p[(k + 1) % p.length];
                s += a[0] * b[1] - b[0] * a[1];
            }
            area += Math.abs(s) / 2;
        }
    }
    assert.ok(nfaces > 0, 'section is non-empty');
    assert.ok(Math.abs(area - 1) < 1e-9, `section area ${area} != 1`);
    assert.ok('slice:parent_cell' in sectioned.cell_data, 'parent ids attached');
    // A plane missing the mesh yields an empty section.
    const empty = m.slice(cube, [0, 0, 5], [0, 0, 1], false);
    assert.equal(empty.cells.length, 0);
    assert.throws(() => m.slice(cube, [0, 0, 0], [0, 0, 0]));
});

step('isosurface: the f = z level set of the unit cube has area 1', () => {
    // f = z on the cube's 8 corners, so the isosurface at z = 0.5 is the same
    // unit square slice() cuts there -- the data-driven route to it.
    const field = { ...cube, point_data: { f: new Float64Array([0, 0, 0, 0, 1, 1, 1, 1]) } };
    const iso = m.isosurface(field, 'f', 0.5, -1, true);
    let area = 0;
    let nfaces = 0;
    for (const block of iso.cells) {
        const npc = block.nodesPerCell;
        const data = block.data;
        for (let c = 0; c < data.length / npc; ++c) {
            nfaces += 1;
            let s = 0;
            for (let k = 0; k < npc; ++k) {
                const a = data[c * npc + k];
                const b = data[c * npc + ((k + 1) % npc)];
                s += iso.points[a * 3] * iso.points[b * 3 + 1] -
                     iso.points[b * 3] * iso.points[a * 3 + 1];
            }
            area += Math.abs(s) / 2;
        }
    }
    assert.ok(nfaces > 0, 'contour is non-empty');
    assert.ok(Math.abs(area - 1) < 1e-9, `contour area ${area} != 1`);
    // The contoured field reads back as exactly the isovalue.
    for (const v of iso.point_data.f) assert.equal(v, 0.5);
    assert.ok('iso:parent_cell' in iso.cell_data, 'parent ids attached');
    assert.ok('iso:value' in iso.cell_data && 'iso:index' in iso.cell_data, 'contours tagged');

    // Two isovalues land in one mesh, tagged and in ascending order.
    const two = m.isosurface(field, 'f', [0.75, 0.25]);
    const vals = [].concat(...two.cell_data['iso:value'].map((b) => Array.from(b)));
    assert.deepEqual([...new Set(vals)], [0.25, 0.75], 'ascending, tagged');

    // Out of range is an empty contour, not an error; a cell field has no level set.
    assert.equal(m.isosurface(field, 'f', 5).cells.length, 0);
    assert.throws(() => m.isosurface(cube, 'nope', 0.5));
});

step('partition: refined hexahedra decompose into 2 balanced pieces', () => {
    const grid = m.refine(cube);  // 8 hexahedra
    const pieces = m.partition(grid, 2);
    assert.equal(pieces.length, 2);
    assert.equal(pieces[0].partId, 0);
    assert.equal(pieces[1].partId, 1);
    let total = 0;
    for (const piece of pieces) {
        // Blocks are kept 1:1 with the input (unlike split).
        assert.equal(piece.mesh.cells.length, grid.cells.length);
        total += piece.mesh.cells[0].data.length / 8;
    }
    assert.equal(total, 8);

    const labels = m.partitionLabels(grid, 2);
    assert.equal(labels.length, 1);
    assert.equal(labels[0].length, 8);
    for (const p of labels[0]) {
        assert.ok(p >= 0 && p < 2);
    }
});

step('partition: kahip is compiled out of the WASM build and says so', () => {
    assert.throws(() => m.partition(cube, 2, 'kahip'), /MESHIOPLUSPLUS_WITH_KAHIP/);
});

step('every binding is reachable through the wrapper', () => {
    // Regression guard for the v7.2.1 class of bug: every binding must be
    // forwarded by src/index.mjs, not merely bound in js_bindings.cpp.
    // This list is exhaustive on purpose -- it used to cover only the geometry
    // ops, which left the I/O and data surfaces unguarded.
    for (const name of [
        'readMesh',
        'readMeshSelective',
        'readMetadata',
        'readerSupportsOptions',
        'writeMesh',
        'convert',
        'convertSurface',
        'convertSurfaceOps',
        'numNodesPerCell',
        'topologicalDimension',
        'availableFormats',
        'dataDrop',
        'dataKeep',
        'dataRename',
        'dataPointToCell',
        'dataCellToPoint',
        'dataCalc',
        'dataCondition',
        'dataInfo',
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
        'smooth',
        'interpolate',
        'slice',
        'isosurface',
        'cropBbox',
        'cropPlane',
        'split',
        'convertCells',
        'refine',
        'decimate',
        'partition',
        'partitionLabels',
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

step('availableFormats reports what this build can read and write', () => {
    const { readers, writers } = m.availableFormats();
    assert.ok(Array.isArray(readers) && Array.isArray(writers));
    assert.ok(readers.includes('vtu') && writers.includes('vtu'));
    // The two lists genuinely differ: openfoam is read-only, svg/tikz are
    // write-only. A viewer that assumes one list would offer broken menu items.
    assert.ok(readers.includes('openfoam') && !writers.includes('openfoam'));
    assert.ok(writers.includes('svg') && !readers.includes('svg'));
    // No HDF5/netCDF-backed format is in a wasm build.
    assert.ok(!readers.includes('med') && !writers.includes('med'));
    // Sorted, so a UI can render them without sorting again.
    assert.deepEqual(readers, [...readers].sort());
    assert.deepEqual(writers, [...writers].sort());
});

step('convertSurface turns a volume mesh into its renderable boundary', () => {
    m.writeMesh('/cube.vtu', cube);
    m.convertSurface('/cube.vtu', '/cube-surf.vtp');
    const surf = m.readMesh('/cube-surf.vtp');
    // One hexahedron -> 6 boundary quads, and no 3D cells left.
    assert.equal(surf.cells.length, 1);
    assert.equal(surf.cells[0].data.length, 6 * 4);
    assert.equal(surf.points.length, 8 * 3);
});

step('convertSurface carries cell data onto the boundary facets', () => {
    // extract_surface drops cell data (a facet is not a cell), but for a
    // viewer the useful answer is the owning cell's value -- colouring a solid
    // by its material tag is the common case, and without the parent-id gather
    // the array would simply vanish on the way to the renderer.
    const tagged = { ...cube, cell_data: { material: [new Float64Array([7])] } };
    m.writeMesh('/tagged.vtu', tagged);
    m.convertSurface('/tagged.vtu', '/tagged-surf.vtp');
    const surf = m.readMesh('/tagged-surf.vtp');
    assert.deepEqual(Array.from(surf.cell_data.material[0]), [7, 7, 7, 7, 7, 7]);
    // The provenance array is plumbing; it must not clutter a colour-by menu.
    assert.equal(surf.cell_data['surface:parent_cell'], undefined);
});

step('convertSurface passes a surface mesh through, linearized', () => {
    const tri6 = {
        points: new Float64Array([0, 0, 0, 1, 0, 0, 0, 1, 0, 0.5, 0, 0, 0.5, 0.5, 0, 0, 0.5, 0]),
        dim: 3,
        cells: [{ type: 'triangle6', data: new Int32Array([0, 1, 2, 3, 4, 5]), nodesPerCell: 6 }],
        point_data: {},
        cell_data: {},
        field_data: {},
    };
    m.writeMesh('/tri6.vtu', tri6);
    m.convertSurface('/tri6.vtu', '/tri6-surf.vtp');
    const out = m.readMesh('/tri6-surf.vtp');
    // Linearized to a 3-node triangle: a renderer has no mid-side nodes, and
    // drawing the 6-node connectivity verbatim would be visible garbage.
    assert.equal(out.cells[0].data.length, 3);
});

step('convertSurface keeps data the flat JS mesh would have dropped', () => {
    // The reason this binding exists rather than readMesh -> extractSkin ->
    // writeMesh: a 3-component array cannot survive the JS mesh boundary.
    const vtu = `<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
 <UnstructuredGrid><Piece NumberOfPoints="4" NumberOfCells="1">
  <Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">
   0 0 0  1 0 0  0 1 0  0 0 1</DataArray></Points>
  <PointData>
   <DataArray type="Float64" Name="disp" NumberOfComponents="3" format="ascii">
    1 2 3  4 5 6  7 8 9  10 11 12</DataArray>
  </PointData>
  <Cells>
   <DataArray type="Int64" Name="connectivity" format="ascii">0 1 2 3</DataArray>
   <DataArray type="Int64" Name="offsets" format="ascii">4</DataArray>
   <DataArray type="UInt8" Name="types" format="ascii">10</DataArray>
  </Cells>
 </Piece></UnstructuredGrid>
</VTKFile>`;
    m.FS.writeFile('/vec.vtu', vtu);
    m.convertSurface('/vec.vtu', '/vec-surf.vtp');
    const text = new TextDecoder().decode(m.FS.readFile('/vec-surf.vtp'));
    assert.match(text, /Name="disp"/);
    assert.match(text, /NumberOfComponents="3"/);
});

step('convertSurfaceOps with an empty pipeline is exactly convertSurface', () => {
    // The viewer calls this for both the plain and the post-operation display,
    // so if the two ever diverged, applying and then undoing an operation
    // would silently change the picture.
    m.writeMesh('/eq.vtu', cube);
    m.convertSurface('/eq.vtu', '/eq-a.vtp');
    const report = m.convertSurfaceOps('/eq.vtu', '/eq-b.vtp', []);
    assert.deepEqual(
        Array.from(m.FS.readFile('/eq-a.vtp')),
        Array.from(m.FS.readFile('/eq-b.vtp'))
    );
    assert.deepEqual(report.steps, []);
    assert.deepEqual(report.warnings, []);
});

step('convertSurfaceOps runs each operation and reports its counters', () => {
    m.writeMesh('/ops.vtu', cube);

    const q = m.convertSurfaceOps('/ops.vtu', '/ops-q.vtp', [{ op: 'quality' }]);
    assert.equal(q.steps[0].op, 'quality');
    assert.ok('quality:scaled_jacobian' in m.readMesh('/ops-q.vtp').cell_data);

    const c = m.convertSurfaceOps('/ops.vtu', '/ops-c.vtp', [{ op: 'clean', weld: true }]);
    assert.equal(typeof c.steps[0].pointsWelded, 'number');

    const s = m.convertSurfaceOps('/ops.vtu', '/ops-s.vtp', [
        { op: 'smooth', method: 'laplacian', iterations: 3, fixBoundary: false },
    ]);
    assert.equal(typeof s.steps[0].numNodesMoved, 'number');

    const r = m.convertSurfaceOps('/ops.vtu', '/ops-r.vtp', [{ op: 'refine', levels: 1 }]);
    assert.equal(r.steps[0].op, 'refine');
    // One hexahedron refines to 8, whose boundary is 6*4 = 24 quads.
    assert.equal(m.readMesh('/ops-r.vtp').cells[0].data.length, 24 * 4);

    const p = m.convertSurfaceOps('/ops.vtu', '/ops-p.vtp', [{ op: 'partition', nparts: 2 }]);
    assert.ok('partition:part' in m.readMesh('/ops-p.vtp').cell_data);

    // Operations compose, in order.
    const both = m.convertSurfaceOps('/ops.vtu', '/ops-both.vtp', [
        { op: 'refine', levels: 1 },
        { op: 'quality' },
    ]);
    assert.deepEqual(both.steps.map((x) => x.op), ['refine', 'quality']);
});

step('convertSurfaceOps takes a planar cross-section (slice)', () => {
    m.writeMesh('/sec.vtu', cube);
    // The mid-plane cross-section of the unit cube is a non-empty surface (a
    // unit square of section faces at z = 0.5), rendered directly.
    const r = m.convertSurfaceOps('/sec.vtu', '/sec.vtp', [
        { op: 'section', point: [0.5, 0.5, 0.5], normal: [0, 0, 1] },
    ]);
    assert.equal(r.steps[0].op, 'section');
    assert.equal(r.warnings.length, 0, 'a plane through the mesh must not warn');
    const sec = m.readMesh('/sec.vtp');
    assert.ok(sec.cells.length > 0, 'the section is non-empty');
    for (const block of sec.cells)
        for (let i = 2; i < sec.points.length; i += 3)
            assert.ok(Math.abs(sec.points[i] - 0.5) < 1e-9, 'section lies at z = 0.5');

    // A plane that misses the mesh yields an empty section, and that is warned.
    const miss = m.convertSurfaceOps('/sec.vtu', '/sec-miss.vtp', [
        { op: 'section', point: [0.5, 0.5, 5], normal: [0, 0, 1] },
    ]);
    assert.ok(miss.warnings.length > 0, 'an empty section must warn');
});

step('convertSurfaceOps keeps multi-component data through an operation', () => {
    // The whole reason this binding exists: the same pipeline expressed as
    // readMesh -> smooth -> writeMesh would flatten `disp` on the way through.
    const vtu = `<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
 <UnstructuredGrid><Piece NumberOfPoints="4" NumberOfCells="1">
  <Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">
   0 0 0  1 0 0  0 1 0  0 0 1</DataArray></Points>
  <PointData>
   <DataArray type="Float64" Name="disp" NumberOfComponents="3" format="ascii">
    1 2 3  4 5 6  7 8 9  10 11 12</DataArray>
  </PointData>
  <Cells>
   <DataArray type="Int64" Name="connectivity" format="ascii">0 1 2 3</DataArray>
   <DataArray type="Int64" Name="offsets" format="ascii">4</DataArray>
   <DataArray type="UInt8" Name="types" format="ascii">10</DataArray>
  </Cells>
 </Piece></UnstructuredGrid>
</VTKFile>`;
    m.FS.writeFile('/vec-ops.vtu', vtu);
    m.convertSurfaceOps('/vec-ops.vtu', '/vec-ops.vtp', [{ op: 'quality' }]);
    const text = new TextDecoder().decode(m.FS.readFile('/vec-ops.vtp'));
    assert.match(text, /Name="disp"/);
    assert.match(text, /NumberOfComponents="3"/);
});

step('convertSurfaceOps can keep the provenance array for a picker', () => {
    m.writeMesh('/prov.vtu', cube);
    m.convertSurfaceOps('/prov.vtu', '/prov.vtp', [], { keepProvenance: true });
    assert.ok('surface:parent_cell' in m.readMesh('/prov.vtp').cell_data);
    m.convertSurfaceOps('/prov.vtu', '/prov2.vtp', []);
    assert.equal(m.readMesh('/prov2.vtp').cell_data['surface:parent_cell'], undefined);
});

step('convertSurfaceOps rejects an unknown operation by name', () => {
    m.writeMesh('/bad.vtu', cube);
    assert.throws(
        () => m.convertSurfaceOps('/bad.vtu', '/bad.vtp', [{ op: 'teleport' }]),
        /unknown operation 'teleport'/
    );
});

step('sniffFormat identifies a file by its leading bytes', () => {
    // Deliberately misleading extension: sniffing must go by content.
    m.writeMesh('/sniffme.vtu', tet);
    m.FS.writeFile('/sniffme.dat', m.FS.readFile('/sniffme.vtu'));
    assert.equal(m.sniffFormat('/sniffme.dat'), 'vtu');
    // Only a confident signature match is claimed; anything else returns "".
    m.FS.writeFile('/ambiguous.dat', 'nothing recognizable here\n');
    assert.equal(m.sniffFormat('/ambiguous.dat'), '');
});

if (failed) {
    console.error('\nSMOKE TEST FAILED');
    process.exit(1);
}
console.log('\nSMOKE TEST PASSED');
