import assert from 'node:assert/strict';
import { test } from 'node:test';

import { buildCellTopology, buildPointCells, cellCentroids, cellNormals, cellSection } from '../../../src/viewer/src/selection/topology.ts';
import {
    combine, connectedMask, countMask, growMask, idsOf, invertMask, maskFromIds, maskStats, shrinkMask,
    thresholdCellValues, thresholdPointValues,
} from '../../../src/viewer/src/selection/mask.ts';
import { cellsInShape, facingMask, pointInPolygon, projectPoints } from '../../../src/viewer/src/selection/project.ts';

// A 3x1 strip of unit squares in the z = 0 plane: points 0..7, quads 0..2,
// plus one extra vertex cell and one line cell to exercise section numbering.
const coords = new Float64Array([
    0, 0, 0, 1, 0, 0, 2, 0, 0, 3, 0, 0,
    0, 1, 0, 1, 1, 0, 2, 1, 0, 3, 1, 0,
]);
const polys = [4, 0, 1, 5, 4, 4, 1, 2, 6, 5, 4, 2, 3, 7, 6];
const lines = [2, 0, 3];
const verts = [1, 7];

test('topology numbers verts, then lines, then polys', () => {
    const topo = buildCellTopology(verts, lines, polys);
    assert.equal(topo.numCells, 5);
    assert.deepEqual(topo.sections, { verts: 1, lines: 1, polys: 3 });
    assert.equal(cellSection(topo, 0), 'verts');
    assert.equal(cellSection(topo, 1), 'lines');
    assert.equal(cellSection(topo, 2), 'polys');
    assert.deepEqual(Array.from(topo.points.slice(topo.offsets[2], topo.offsets[3])), [0, 1, 5, 4]);
});

test('centroids and normals of the strip', () => {
    const topo = buildCellTopology([], [], polys);
    const c = cellCentroids(coords, topo);
    assert.deepEqual(Array.from(c.slice(0, 3)), [0.5, 0.5, 0]);
    assert.deepEqual(Array.from(c.slice(6, 9)), [2.5, 0.5, 0]);
    const n = cellNormals(coords, topo);
    for (let i = 0; i < 3; i++) assert.ok(Math.abs(n[i * 3 + 2] - 1) < 1e-6, 'counter-clockwise quad faces +z');
});

test('combine replaces, adds, subtracts and toggles', () => {
    const a = maskFromIds(5, [0, 1]);
    const b = maskFromIds(5, [1, 2]);
    assert.deepEqual(Array.from(idsOf(combine(a, b, 'replace'))), [1, 2]);
    assert.deepEqual(Array.from(idsOf(combine(a, b, 'add'))), [0, 1, 2]);
    assert.deepEqual(Array.from(idsOf(combine(a, b, 'subtract'))), [0]);
    assert.deepEqual(Array.from(idsOf(combine(a, b, 'toggle'))), [0, 2]);
    assert.equal(countMask(invertMask(a)), 3);
    assert.equal(countMask(maskFromIds(3, [-1, 1, 99])), 1);
});

test('grow, shrink and connected walk across shared points', () => {
    const topo = buildCellTopology([], [], polys);
    const pc = buildPointCells(topo, 8);
    const first = maskFromIds(3, [0]);
    assert.deepEqual(Array.from(idsOf(growMask(first, topo, pc))), [0, 1]);
    assert.deepEqual(Array.from(idsOf(growMask(growMask(first, topo, pc), topo, pc))), [0, 1, 2]);
    assert.deepEqual(Array.from(idsOf(connectedMask(first, topo, pc))), [0, 1, 2]);
    // Cell 1 touches the unselected cell 2, so it is the one that goes.
    assert.deepEqual(Array.from(idsOf(shrinkMask(maskFromIds(3, [0, 1]), topo, pc))), [0]);
    assert.equal(countMask(shrinkMask(maskFromIds(3, [0, 1, 2]), topo, pc)), 3);
});

test('statistics skip non-finite values and support magnitude', () => {
    const v = [1, 2, 3, NaN];
    const s = maskStats(v, 1, 0, maskFromIds(4, [0, 1, 2, 3]));
    assert.deepEqual(s, { count: 3, min: 1, max: 3, mean: 2 });
    const vec = [3, 4, 0, 0, 0, 0];
    assert.equal(maskStats(vec, 3, -1, maskFromIds(2, [0])).max, 5);
    assert.equal(maskStats(v, 1, 0, maskFromIds(4, [])).count, 0);
});

test('thresholds on cell values and on point values', () => {
    assert.deepEqual(Array.from(idsOf(thresholdCellValues([0, 5, 10, NaN], 1, 0, 4, 10, 4))), [1, 2]);
    assert.deepEqual(Array.from(idsOf(thresholdCellValues([0, 5, 10], 1, 0, 10, 4, 3))), [1, 2], 'bounds may be given swapped');
    const topo = buildCellTopology([], [], polys);
    const pv = [0, 0, 0, 0, 1, 1, 1, 1]; // top row high
    assert.deepEqual(Array.from(idsOf(thresholdPointValues(pv, 1, 0, 0.5, 2, topo, 'all'))), []);
    assert.deepEqual(Array.from(idsOf(thresholdPointValues(pv, 1, 0, 0.5, 2, topo, 'any'))), [0, 1, 2]);
});

test('projection matches the camera and rect/polygon tests select cells', () => {
    // The row-major matrix from a vtk.js camera at (2,3,10) looking at the
    // origin, aspect 1.5 (see the spike in the plan): the point (1, .5, .2)
    // lands at ndc (0.2316, 0.1362).
    const m = [2.4397, 0, -0.4879, 0, -0.2066, 3.5803, -1.0328, 0, -0.1958, -0.2937, -0.9791, 9.0232, -0.1881, -0.2822, -0.9407, 10.6301];
    const p = projectPoints([1, 0.5, 0.2], m, 600, 400);
    assert.ok(Math.abs(p.x[0] - ((0.2316 + 1) / 2) * 600) < 0.5);
    assert.ok(Math.abs(p.y[0] - ((1 - 0.1362) / 2) * 400) < 0.5);
    assert.equal(p.inFront[0], 1);

    // An orthographic-like identity view: x, y map straight to the canvas.
    const ident = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1];
    const pts = [-0.5, 0.5, 0, 0.5, 0.5, 0, 0.5, -0.5, 0, -0.9, -0.9, 0];
    const proj = projectPoints(pts, ident, 100, 100);
    const left = cellsInShape(proj, { kind: 'rect', rect: { x0: 0, y0: 0, x1: 50, y1: 50 } });
    assert.deepEqual(Array.from(idsOf(left)), [0]);
    const tri = { kind: 'polygon', points: [0, 0, 100, 0, 100, 100] };
    // Upper-right triangle of the canvas holds the cells at (75,25) and (75,75)'s diagonal side.
    assert.ok(countMask(cellsInShape(proj, tri)) >= 1);
    const eligible = maskFromIds(4, [1]);
    assert.deepEqual(Array.from(idsOf(cellsInShape(proj, { kind: 'rect', rect: { x0: 0, y0: 0, x1: 100, y1: 100 } }, eligible))), [1]);
});

test('a point behind the camera is never selected', () => {
    const m = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, -1, 0]; // w = -z
    const proj = projectPoints([0, 0, 1], m, 10, 10); // w = -1 -> behind
    assert.equal(proj.inFront[0], 0);
});

test('point-in-polygon', () => {
    const sq = [0, 0, 10, 0, 10, 10, 0, 10];
    assert.ok(pointInPolygon(5, 5, sq));
    assert.ok(!pointInPolygon(11, 5, sq));
    const concave = [0, 0, 10, 0, 10, 10, 5, 4, 0, 10];
    assert.ok(!pointInPolygon(5, 8, concave));
    assert.ok(pointInPolygon(5, 2, concave));
});

test('facing keeps cells whose normal points at the eye', () => {
    const cent = [0, 0, 0, 0, 0, 0, 5, 5, 5];
    const norm = [0, 0, 1, 0, 0, -1, 0, 0, 0];
    const f = facingMask(cent, norm, [0, 0, 10], null);
    assert.deepEqual(Array.from(f), [1, 0, 1]);
    const g = facingMask(cent, norm, null, [0, 0, -1]);
    assert.deepEqual(Array.from(g), [1, 0, 1]);
});
