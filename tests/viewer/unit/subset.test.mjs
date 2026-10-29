import assert from 'node:assert/strict';
import { test } from 'node:test';

import { subsetCells, takeTuples } from '../../../src/viewer/src/render/subset.ts';

const verts = [1, 7];
const lines = [2, 0, 3];
const polys = [3, 0, 1, 2, 4, 1, 2, 3, 4];

test('subsetCells keeps the chosen cells in each section and maps them back', () => {
    // cells: 0 vert, 1 line, 2 tri, 3 quad
    const s = subsetCells(verts, lines, polys, [1, 0, 0, 1]);
    assert.deepEqual(Array.from(s.verts), [1, 7]);
    assert.deepEqual(Array.from(s.lines), []);
    assert.deepEqual(Array.from(s.polys), [4, 1, 2, 3, 4]);
    assert.deepEqual(Array.from(s.cellMap), [0, 3]);
});

test('an empty mask gives an empty subset, a full one the original', () => {
    const none = subsetCells(verts, lines, polys, [0, 0, 0, 0]);
    assert.equal(none.cellMap.length, 0);
    const all = subsetCells(verts, lines, polys, [1, 1, 1, 1]);
    assert.deepEqual(Array.from(all.polys), polys);
    assert.deepEqual(Array.from(all.cellMap), [0, 1, 2, 3]);
});

test('takeTuples restricts a per-cell array, keeping its type', () => {
    const v = new Float32Array([1, 2, 3, 4, 5, 6, 7, 8]);
    const out = takeTuples(v, 2, [3, 0]);
    assert.ok(out instanceof Float32Array);
    assert.deepEqual(Array.from(out), [7, 8, 1, 2]);
});
