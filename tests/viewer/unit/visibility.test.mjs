import assert from 'node:assert/strict';
import { test } from 'node:test';

import { Visibility } from '../../../src/viewer/src/app/visibility.ts';

test('layers compose with AND and clear independently', () => {
    const applied = [];
    const v = new Visibility((m) => applied.push(m ? Array.from(m) : null));
    assert.equal(v.combined(), null);
    v.set('a', Uint8Array.from([1, 1, 0, 1]));
    assert.deepEqual(applied.at(-1), [1, 1, 0, 1]);
    v.set('b', Uint8Array.from([0, 1, 1, 1]));
    assert.deepEqual(applied.at(-1), [0, 1, 0, 1]);
    assert.equal(v.hiddenCount(), 2);
    v.set('a', null);
    assert.deepEqual(applied.at(-1), [0, 1, 1, 1]);
    v.set('b', null);
    assert.equal(applied.at(-1), null);
    assert.equal(v.hiddenCount(), 0);
});

test('listeners hear about every change and can unsubscribe', () => {
    const v = new Visibility(() => {});
    let n = 0;
    const off = v.onChange(() => n++);
    v.set('a', Uint8Array.from([1, 0]));
    v.reset();
    off();
    v.set('a', Uint8Array.from([1, 0]));
    assert.equal(n, 2);
});

test('pruneTo drops layers made for a different surface', () => {
    const v = new Visibility(() => {});
    v.set('a', Uint8Array.from([1, 0, 1]));
    v.set('b', Uint8Array.from([1, 1]));
    v.pruneTo(3);
    assert.ok(v.get('a'));
    assert.equal(v.get('b'), null);
});
