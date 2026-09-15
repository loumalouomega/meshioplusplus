/**
 * Overview cards: construction from a manifest's text (a broken manifest is
 * a flagged card, never a missing one), sorting, and the fixed split ->
 * colour assignment.
 *
 * Run with `node --test` via tests/viewer/unit/build-and-run.mjs.
 */
import assert from 'node:assert/strict';
import { test } from 'node:test';

import {
    OTHER_COLOR,
    SPLIT_COLORS,
    UNASSIGNED_COLOR,
    balanceOf,
    buildCard,
    sortCards,
    splitPaletteOf,
} from '../../../src/viewer/src/dataset/overview.ts';

const GOOD = JSON.stringify({
    Version: 1,
    Name: 'campaign',
    Description: 'Re sweep',
    Entries: [
        { Id: 'a', Source: { Path: 'a.vtu' }, Split: 'train', Tags: ['re100', 'coarse'], Group: 'cyl/laminar' },
        { Id: 'b', Source: { Path: 'b.vtu' }, Split: 'train', Tags: ['re100'] },
        { Id: 'c', Source: { Path: 'c.vtu' }, Split: 'test' },
        { Id: 'd', Source: { Path: 'd.vtu' }, Group: 'cyl/turbulent' },
    ],
});

test('a card summarizes a manifest', () => {
    const card = buildCard('dataset.json', GOOD, 1700000000000);
    assert.equal(card.name, 'campaign');
    assert.equal(card.description, 'Re sweep');
    assert.equal(card.numEntries, 4);
    assert.deepEqual(card.splits, { train: 2, test: 1, '': 1 });
    assert.deepEqual(card.tags, ['coarse', 're100']);
    assert.deepEqual(card.groups, ['cyl/laminar', 'cyl/turbulent']);
    assert.equal(card.lastModified, 1700000000000);
    assert.equal(card.parseError, null);
    assert.equal(card.health, null);
    assert.equal(card.dirty, false);
});

test('a broken manifest is a flagged card, never a missing one', () => {
    const card = buildCard('broken.json', '{"Version": 1, "Bogus": 1}', null);
    assert.equal(card.numEntries, 0);
    assert.match(card.parseError, /unknown key 'Bogus'/);
    const notJson = buildCard('nope.json', '{', null);
    assert.match(notJson.parseError, /not valid JSON/);
});

test('cards sort by name, entries, modification time and health', () => {
    const a = buildCard('a.json', GOOD, 100);
    a.name = 'zeta';
    const b = buildCard('b.json', GOOD, 300);
    b.name = 'alpha';
    b.numEntries = 1;
    const broken = buildCard('c.json', '{', 200);
    const bad = buildCard('d.json', GOOD, null);
    bad.name = 'mid';
    bad.health = { badEntries: ['x'] };
    const cards = [a, b, broken, bad];
    assert.deepEqual(sortCards(cards, 'name').map((c) => c.path), ['b.json', 'c.json', 'd.json', 'a.json']);
    assert.deepEqual(sortCards(cards, 'entries').map((c) => c.path), ['d.json', 'a.json', 'b.json', 'c.json']);
    assert.deepEqual(sortCards(cards, 'modified').map((c) => c.path), ['b.json', 'c.json', 'a.json', 'd.json']);
    // worst first: broken, bad, unscanned, then healthy (none here)
    assert.deepEqual(sortCards(cards, 'health').map((c) => c.path), ['c.json', 'd.json', 'b.json', 'a.json']);
    // the input is not reordered
    assert.deepEqual(cards.map((c) => c.path), ['a.json', 'b.json', 'c.json', 'd.json']);
});

test('split colours are assigned in fixed order, never cycled', () => {
    const palette = splitPaletteOf(['valid', 'train', '', 'test', 'train']);
    assert.equal(palette.get('test'), SPLIT_COLORS[0]);
    assert.equal(palette.get('train'), SPLIT_COLORS[1]);
    assert.equal(palette.get('valid'), SPLIT_COLORS[2]);
    assert.equal(palette.get(''), UNASSIGNED_COLOR);
    const many = splitPaletteOf(Array.from({ length: 10 }, (_, i) => `s${String(i).padStart(2, '0')}`));
    assert.equal(many.get('s07'), SPLIT_COLORS[7]);
    assert.equal(many.get('s08'), OTHER_COLOR);
    assert.equal(many.get('s09'), OTHER_COLOR);
});

test('a card balance keeps unassigned last', () => {
    assert.deepEqual(balanceOf({ '': 1, valid: 1, train: 2 }), [
        { split: 'train', count: 2, fraction: 0.5 },
        { split: 'valid', count: 1, fraction: 0.25 },
        { split: '', count: 1, fraction: 0.25 },
    ]);
    assert.deepEqual(balanceOf({}), []);
});
