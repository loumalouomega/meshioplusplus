/**
 * Manifest diffing: the entry-keyed structural comparison and the LCS line
 * diff over the stable serialization.
 *
 * Run with `node --test` via tests/viewer/unit/build-and-run.mjs.
 */
import assert from 'node:assert/strict';
import { test } from 'node:test';

import {
    MAX_LINE_DIFF_CELLS,
    diffLines,
    diffManifests,
    stableStringify,
    summarizeDiff,
} from '../../../src/viewer/src/dataset/diff.ts';
import { parseManifest, stringifyManifest } from '../../../src/viewer/src/dataset/manifest.ts';

const A = parseManifest({
    Version: 1,
    Name: 'campaign',
    Metadata: { solver: 'kratos', mesh: 'coarse' },
    Entries: [
        { Id: 'c42', Source: { Path: 'c42.vtu' }, Split: 'train', Tags: ['a'], Metadata: { Re: 100, Ma: 0.3 } },
        { Id: 'pair', Source: { Paths: ['a.vtu', 'b.vtu'] } },
        { Id: 'same', Source: { Path: 's.vtu' }, Notes: 'kept' },
    ],
});

const B = parseManifest({
    Version: 1,
    Name: 'campaign v2',
    Metadata: { mesh: 'coarse', solver: 'kratos' }, // same content, other key order
    Entries: [
        { Id: 'c42', Source: { Path: 'c42.vtu' }, Split: 'valid', Tags: ['a', 'b'], Metadata: { Ma: 0.3, Re: 100 } },
        { Id: 'same', Source: { Path: 's.vtu' }, Notes: 'kept' },
        { Id: 'new', Source: { Path: 'n.vtu' } },
    ],
});

test('the structural diff is keyed by entry id and ignores key order', () => {
    const diff = diffManifests(A, B);
    assert.deepEqual(diff.header, [{ field: 'name', before: 'campaign', after: 'campaign v2' }]);
    assert.deepEqual(diff.added, ['new']);
    assert.deepEqual(diff.removed, ['pair']);
    assert.equal(diff.changed.length, 1);
    assert.equal(diff.changed[0].id, 'c42');
    assert.deepEqual(
        diff.changed[0].fields.map((f) => f.field),
        ['split', 'tags'], // metadata is equal despite the key order
    );
    assert.deepEqual(summarizeDiff('A', 'B', diff), {
        a: 'A',
        b: 'B',
        headerChanged: 1,
        added: 1,
        removed: 1,
        changed: 1,
    });
    assert.deepEqual(diffManifests(A, A), { header: [], added: [], removed: [], changed: [] });
});

test('stableStringify sorts keys at every level', () => {
    assert.equal(stableStringify({ b: [{ z: 1, a: 2 }], a: null }), '{"a":null,"b":[{"a":2,"z":1}]}');
});

test('the line diff marks the edited lines only', () => {
    const lines = diffLines('a\nb\nc\nd', 'a\nx\nc\nd\ne');
    assert.deepEqual(
        lines.map((l) => `${l.kind}:${l.text}`),
        ['same:a', 'del:b', 'add:x', 'same:c', 'same:d', 'add:e'],
    );
    assert.ok(diffLines('same\ntext', 'same\ntext').every((l) => l.kind === 'same'));
    // a real manifest edit: one split changed
    const before = stringifyManifest(A);
    const doc = JSON.parse(before);
    doc.Entries = doc.Entries.map((e) => (e.Id === 'c42' ? { ...e, Split: 'valid' } : e));
    const after = stringifyManifest(parseManifest(doc));
    const changed = diffLines(before, after).filter((l) => l.kind !== 'same');
    assert.deepEqual(
        changed.map((l) => `${l.kind}:${l.text.trim()}`),
        ['del:"Split": "train",', 'add:"Split": "valid",'],
    );
});

test('a diff too large for the LCS table is declined, not attempted', () => {
    const n = Math.ceil(Math.sqrt(MAX_LINE_DIFF_CELLS)) + 1;
    const big = Array.from({ length: n }, (_, i) => `line ${i}`).join('\n');
    assert.equal(diffLines(big, `${big}\nextra`), null);
});
