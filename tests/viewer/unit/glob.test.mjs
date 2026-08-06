/**
 * The TS glob twin must accept exactly the sequence glob language
 * (doc/sequences.md): `*`/`?` only, `[`/`]` literal, no `**`, directory
 * part literal. Run with `node --test` via build-and-run.mjs.
 */
import assert from 'node:assert/strict';
import { test } from 'node:test';

import {
    globMatch,
    isPattern,
    nameMatch,
    naturalCompare,
} from '../../../src/viewer/src/dataset/glob.ts';

test('star and question-mark semantics', () => {
    assert.ok(nameMatch('out_*.vtu', 'out_0001.vtu'));
    assert.ok(nameMatch('out_*.vtu', 'out_.vtu')); // * matches the empty run
    assert.ok(nameMatch('out_????.vtu', 'out_0001.vtu'));
    assert.ok(!nameMatch('out_???.vtu', 'out_0001.vtu'));
    assert.ok(nameMatch('*', 'anything'));
    assert.ok(nameMatch('a*b*c', 'aXXbYYc'));
    assert.ok(!nameMatch('a*b*c', 'aXXbYY'));
});

test('brackets are literal characters, not sets', () => {
    assert.ok(nameMatch('a[bc].vtu', 'a[bc].vtu'));
    assert.ok(!nameMatch('a[bc].vtu', 'ab.vtu'));
});

test('the directory part is literal', () => {
    assert.ok(globMatch('runs/c42/out_*.vtu', 'runs/c42/out_7.vtu'));
    assert.ok(!globMatch('runs/c42/out_*.vtu', 'runs/c43/out_7.vtu'));
    assert.ok(!globMatch('runs/*/out_1.vtu', 'runs/c42/out_1.vtu')); // no dir globbing
    assert.ok(globMatch('out_*.vtu', 'out_7.vtu'));
    assert.ok(!globMatch('out_*.vtu', 'sub/out_7.vtu')); // * never crosses '/'
});

test('isPattern', () => {
    assert.ok(isPattern('a*.vtu'));
    assert.ok(isPattern('a?.vtu'));
    assert.ok(!isPattern('a[0].vtu'));
});

test('natural-numeric ordering', () => {
    const names = ['out_10.vtu', 'out_9.vtu', 'out_100.vtu', 'out_2.vtu'];
    assert.deepEqual(
        [...names].sort(naturalCompare),
        ['out_2.vtu', 'out_9.vtu', 'out_10.vtu', 'out_100.vtu'],
    );
    // huge digit runs never overflow (compared on the digits, never parsed)
    const big = '9'.repeat(40);
    const bigger = `1${'0'.repeat(40)}`;
    assert.ok(naturalCompare(`a${big}`, `a${bigger}`) < 0);
    // zero-padding ties break deterministically
    assert.ok(naturalCompare('out_01', 'out_1') !== 0);
});
