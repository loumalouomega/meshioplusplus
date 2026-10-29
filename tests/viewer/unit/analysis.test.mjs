import assert from 'node:assert/strict';
import { test } from 'node:test';

import { binRange, histogram } from '../../../src/viewer/src/analysis/histogram.ts';
import {
    arrayColumns,
    csvChunks,
    csvEscape,
    filterRows,
    makeTable,
    sortedOrder,
} from '../../../src/viewer/src/analysis/table.ts';

const table = () =>
    makeTable(
        4,
        [{ key: 'id', name: 'id', get: (i) => i }],
        [
            { name: 's', values: [3, NaN, 1, 2], components: 1 },
            { name: 'v', values: [3, 4, 0, 0, 0, 0, 1, 0], components: 2 },
        ]
    );

test('vector arrays get a column per component and a magnitude', () => {
    const cols = arrayColumns({ name: 'v', values: [3, 4, 0, 0], components: 2 });
    assert.deepEqual(cols.map((c) => c.name), ['v_x', 'v_y', '|v|']);
    assert.equal(cols[2].get(0), 5);
    assert.deepEqual(arrayColumns({ name: 'w', values: new Array(8).fill(0), components: 4 }).map((c) => c.name).slice(0, 4), ['w_0', 'w_1', 'w_2', 'w_3']);
});

test('sorting is stable, descending reverses and non-finite values go last either way', () => {
    const t = table();
    assert.deepEqual([...sortedOrder(t, 1, 'asc')], [2, 3, 0, 1]);
    assert.deepEqual([...sortedOrder(t, 1, 'desc')], [0, 3, 2, 1]);
});

test('filtering keeps rows in range and drops non-finite ones', () => {
    assert.deepEqual([...filterRows(table(), 1, 1, 2)], [2, 3]);
});

test('csv has a header, escapes names and blanks non-finite values', () => {
    const t = table();
    const text = [...csvChunks(t, [0, 1], [0, 1])].join('');
    assert.equal(text, 'id,s\n0,3\n1,\n');
    assert.equal(csvEscape('a,"b"'), '"a,""b"""');
    assert.equal([...csvChunks(t, [0, 1, 2, 3], [0], 2)].length, 3, 'header plus two chunks');
});

test('histogram bins, ignores non-finite values and reports statistics', () => {
    const h = histogram([0, 1, 2, 3, 4, NaN], 1, { bins: 4 });
    assert.equal(h.count, 5);
    assert.equal(h.nonFinite, 1);
    assert.deepEqual([...h.counts], [1, 1, 1, 2]);
    assert.equal(h.mean, 2);
    assert.ok(Math.abs(h.std - Math.sqrt(2)) < 1e-12);
    assert.deepEqual(binRange(h, 0), [0, 1]);
    assert.deepEqual(binRange(h, 3), [3, 4]);
});

test('histogram of a vector uses the magnitude, and a constant array fills one bin', () => {
    const h = histogram([3, 4, 0, 0], 2, { bins: 2 });
    assert.equal(h.max, 5);
    const c = histogram([7, 7, 7], 1, { bins: 5 });
    assert.deepEqual([...c.counts], [3, 0, 0, 0, 0]);
    assert.equal(histogram([NaN], 1).count, 0);
});

import { PREF_FIELDS, parsePrefs } from '../../../src/viewer/src/state/prefs.ts';

test('preferences keep known controls with the right type and drop the rest', () => {
    const p = parsePrefs({ edges: true, colormap: 'Viridis', opacity: '55', 'light-cull': 'yes', bogus: 1, 'color-bands': 8 });
    assert.deepEqual(p, { edges: true, colormap: 'Viridis', opacity: '55' });
    assert.deepEqual(parsePrefs(null), {});
    assert.deepEqual(parsePrefs([1, 2]), {});
    assert.deepEqual(parsePrefs({ opacity: 'abc' }), {});
    assert.ok(PREF_FIELDS.findIndex((f) => f.id === 'colormap-safe') < PREF_FIELDS.findIndex((f) => f.id === 'colormap'));
});
