/**
 * Pattern/id suggestion for the add-case flow: a suggested pattern must
 * match the selection exactly against the whole listing — never a file the
 * user did not pick. Run with `node --test` via build-and-run.mjs.
 */
import assert from 'node:assert/strict';
import { test } from 'node:test';

import { suggestSource } from '../../../src/viewer/src/dataset/suggest.ts';

test('one file is a Path', () => {
    assert.deepEqual(suggestSource(['runs/a.vtu'], ['runs/a.vtu']), {
        Path: 'runs/a.vtu',
    });
});

test('a numeric family becomes a Pattern', () => {
    const files = ['c/out_0.vtu', 'c/out_1.vtu', 'c/out_2.vtu', 'c/readme.txt'];
    const picked = ['c/out_1.vtu', 'c/out_0.vtu', 'c/out_2.vtu'];
    assert.deepEqual(suggestSource(picked, files), { Pattern: 'c/out_*.vtu' });
});

test('a pattern that would catch unpicked files degrades to Paths', () => {
    const files = ['c/out_0.vtu', 'c/out_1.vtu', 'c/out_2.vtu'];
    const picked = ['c/out_0.vtu', 'c/out_2.vtu']; // out_1 not picked
    assert.deepEqual(suggestSource(picked, files), {
        Paths: ['c/out_0.vtu', 'c/out_2.vtu'],
    });
});

test('mixed directories degrade to Paths in natural order', () => {
    const picked = ['b/out_10.vtu', 'a/out_2.vtu'];
    assert.deepEqual(suggestSource(picked, picked), {
        Paths: ['a/out_2.vtu', 'b/out_10.vtu'],
    });
});
