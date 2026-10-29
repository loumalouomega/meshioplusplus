import assert from 'node:assert/strict';
import { test } from 'node:test';

import { blockDim, cornerCount, directIndex, groupLabels, groupMask, hasVolume, mapRegions, smallestFirst } from '../../../src/viewer/src/groups/regionMap.ts';

test('cell types are classified by dimension', () => {
    assert.equal(blockDim('vertex'), 0);
    assert.equal(blockDim('line3'), 1);
    assert.equal(blockDim('triangle6'), 2);
    assert.equal(blockDim('quad'), 2);
    assert.equal(blockDim('polygon'), 2);
    assert.equal(blockDim('hexahedron20'), 3);
    assert.equal(blockDim('wedge'), 3);
    assert.ok(hasVolume([{ type: 'triangle', numCells: 4 }, { type: 'tetra', numCells: 1 }]));
    assert.ok(!hasVolume([{ type: 'tetra', numCells: 0 }, { type: 'triangle', numCells: 2 }]));
    assert.equal(cornerCount('triangle6'), 3);
    assert.equal(cornerCount('quad9'), 4);
    assert.equal(cornerCount('tetra'), 0);
});

// Two tets sharing a face, the first four surface cells belong to tet 0 and
// the last two to tet 1 (the layout the real extract_surface gave in the spike),
// plus one boundary triangle (file cell 2) matching surface cell 0.
const blocks = [{ type: 'tetra', numCells: 2 }, { type: 'triangle', numCells: 1 }];
const regions = [
    { name: 'solidA', kind: 'cell', dim: 3, tag: 1, entries: [0] },
    { name: 'solidB', kind: 'cell', dim: 3, tag: 2, entries: [1] },
    { name: 'inlet', kind: 'cell', dim: 2, tag: 3, entries: [2] },
    { name: 'nodes', kind: 'point', dim: 0, tag: 4, entries: [0, 3, 9] },
    { name: 'faces', kind: 'side', dim: 2, tag: 5, entries: [0, 1, 1, 2] },
];

test('skin: regions map through the owner cell, boundary patches through node keys', () => {
    const parent = [0, 0, 0, 1, 1, 1];
    const result = mapRegions({
        mode: 'skin',
        blocks,
        regions,
        numSurfaceCells: 6,
        parent,
        lowerDimFacets: (g) => (g === 2 ? [0] : undefined),
        surfacePointOf: [0, -1, -1, 1, -1],
    });
    assert.ok(result.ok);
    const by = Object.fromEntries(result.groups.filter((g) => g.source === 'region').map((g) => [g.name, g]));
    assert.deepEqual(Array.from(by.solidA.surfaceCells), [0, 1, 2]);
    assert.deepEqual(Array.from(by.solidB.surfaceCells), [3, 4, 5]);
    assert.deepEqual(Array.from(by.inlet.surfaceCells), [0]);
    assert.equal(by.inlet.mapped, 1);
    // Point entries 0 and 3 are on the skin; entry 9 is out of range.
    assert.deepEqual(Array.from(by.nodes.surfacePoints), [0, 1]);
    assert.equal(by.nodes.mapped, 2);
    // Facet sets are listed, not mapped.
    assert.equal(by.faces.surfaceCells, null);
    assert.match(by.faces.note, /facet/);
    // Blocks: the tets cover all six facets, the triangle block matches facet 0.
    const blockGroups = result.groups.filter((g) => g.source === 'block');
    assert.deepEqual(blockGroups.map((g) => g.name), ['tetra', 'triangle']);
    assert.equal(blockGroups[0].surfaceCells.length, 6);
});

test('skin: a lower-dimension region with no matching facet says so', () => {
    const result = mapRegions({ mode: 'skin', blocks, regions: [regions[2]], numSurfaceCells: 6, parent: [0, 0, 0, 1, 1, 1] });
    const inlet = result.groups[0];
    assert.equal(inlet.surfaceCells, null);
    assert.equal(inlet.mapped, 0);
    assert.ok(inlet.note);
});

test('skin: refuses provenance that does not match the surface', () => {
    const result = mapRegions({ mode: 'skin', blocks, regions, numSurfaceCells: 6, parent: [0, 1] });
    assert.ok(!result.ok);
    assert.ok(result.reason);
    assert.ok(!mapRegions({ mode: 'skin', blocks, regions, numSurfaceCells: 6 }).ok);
});

test('direct: surface cells follow VTK section order, vertices then lines then polygons', () => {
    // File order: a triangle block, a line block, a vertex block.
    const direct = [{ type: 'triangle', numCells: 2 }, { type: 'line', numCells: 2 }, { type: 'vertex', numCells: 1 }];
    const idx = directIndex(direct, 5);
    // vertex file cell 4 is surface cell 0; lines 2,3 are 1,2; triangles 0,1 are 3,4.
    assert.deepEqual(Array.from(idx), [3, 4, 1, 2, 0]);
    assert.equal(directIndex(direct, 6), null);

    const result = mapRegions({
        mode: 'direct',
        blocks: direct,
        regions: [{ name: 'tris', kind: 'cell', dim: 2, tag: 1, entries: [0, 1] }],
        numSurfaceCells: 5,
    });
    assert.ok(result.ok);
    assert.deepEqual(Array.from(result.groups[0].surfaceCells), [3, 4]);
});

test('direct: a count mismatch is reported, not guessed', () => {
    const r = mapRegions({ mode: 'direct', blocks: [{ type: 'quad', numCells: 3 }], regions: [], numSurfaceCells: 6 });
    assert.ok(!r.ok);
});

test('direct: point regions map by identity when the point counts agree', () => {
    const r = mapRegions({
        mode: 'direct',
        blocks: [{ type: 'triangle', numCells: 1 }],
        regions: [{ name: 'p', kind: 'point', dim: 0, tag: 1, entries: [0, 2] }],
        numSurfaceCells: 1,
        numFilePoints: 3,
        numSurfacePoints: 3,
    });
    assert.deepEqual(Array.from(r.groups[0].surfacePoints), [0, 2]);
});

test('group labels and masks', () => {
    const result = mapRegions({ mode: 'skin', blocks, regions: regions.slice(0, 2), numSurfaceCells: 6, parent: [0, 0, 0, 1, 1, 1] });
    const labels = groupLabels(result.groups, 6, [0, 1]);
    assert.deepEqual(Array.from(labels), [0, 0, 0, 1, 1, 1]);
    // Labels are group indices; groups not in the order leave -1.
    assert.deepEqual(Array.from(groupLabels(result.groups, 6, [1])), [-1, -1, -1, 1, 1, 1]);
    // Where groups overlap the first in the order wins, so smallest-first lets
    // a patch (group 1) sit on the zone (group 0) that contains it.
    const zone = { surfaceCells: Uint32Array.from([0, 1, 2, 3]) };
    const patch = { surfaceCells: Uint32Array.from([1]) };
    const none = { surfaceCells: null };
    assert.deepEqual(smallestFirst([zone, patch, none]), [2, 1, 0]);
    assert.deepEqual(Array.from(groupLabels([zone, patch, none], 5, smallestFirst([zone, patch, none]))), [0, 1, 0, 0, -1]);
    assert.deepEqual(Array.from(groupMask(result.groups[1], 6)), [0, 0, 0, 1, 1, 1]);
});
