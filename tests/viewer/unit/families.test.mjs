/**
 * The model-family table the launch form is driven from: a pure twin of
 * `physicsnemo/_train.py`'s `_FAMILIES` (the names are pinned equal, in
 * order, from the Python side), and the control-to-argument rules.
 *
 * Run with `node --test` via tests/viewer/unit/build-and-run.mjs.
 */
import assert from 'node:assert/strict';
import { test } from 'node:test';

import {
    FAMILIES,
    FAMILY_NAMES,
    familyArgs,
    familyFor,
    fieldValue,
} from '../../../src/viewer/src/dataset/families.ts';

test('the five families, in the Python order, each in a known block', () => {
    assert.deepEqual(FAMILY_NAMES, ['meshgraphnet', 'srresnet', 'fno', 'afno', 'deeponet']);
    for (const family of FAMILIES) {
        assert.ok(['Graph', 'Grid', 'Operator'].includes(family.block), family.name);
        assert.match(family.panelId, /^t-.*-opts$/);
    }
    assert.equal(familyFor('deeponet').block, 'Operator');
    assert.throws(() => familyFor('gpt'), /unknown model family 'gpt' \(known: meshgraphnet, srresnet, fno, afno, deeponet\)/);
});

test('every control id and argument key is unique within its family, and panels are distinct', () => {
    const panels = new Set(FAMILIES.map((f) => f.panelId));
    assert.equal(panels.size, FAMILIES.length);
    const ids = new Set();
    for (const family of FAMILIES) {
        const keys = family.fields.map((f) => f.key);
        assert.equal(new Set(keys).size, keys.length, family.name);
        for (const field of family.fields) {
            assert.ok(!ids.has(field.id), `${field.id} is reused`);
            ids.add(field.id);
        }
    }
});

test('field values: numbers, triples, pairs, axes, names, and blanks', () => {
    assert.equal(fieldValue({ key: 'n', id: 'x', kind: 'number' }, ' 32 '), 32);
    assert.equal(fieldValue({ key: 'n', id: 'x', kind: 'number' }, ''), undefined);
    assert.throws(() => fieldValue({ key: 'n', id: 'x', kind: 'number', required: true }, ''), /n is required/);
    assert.throws(() => fieldValue({ key: 'n', id: 'x', kind: 'number' }, 'many'), /must be a number/);
    assert.deepEqual(fieldValue({ key: 'resolution', id: 'x', kind: 'triple' }, '16, 8,4'), [16, 8, 4]);
    assert.throws(() => fieldValue({ key: 'resolution', id: 'x', kind: 'triple' }, '16,8'), /3 positive integers/);
    assert.throws(() => fieldValue({ key: 'resolution', id: 'x', kind: 'triple' }, '16,8,0'), /3 positive integers/);
    assert.deepEqual(fieldValue({ key: 'patch_size', id: 'x', kind: 'pair' }, '8,8'), [8, 8]);
    assert.throws(() => fieldValue({ key: 'patch_size', id: 'x', kind: 'pair' }, '8,8,8'), /2 positive integers/);
    assert.equal(fieldValue({ key: 'squeeze', id: 'x', kind: 'axis' }, '2'), 2);
    assert.equal(fieldValue({ key: 'squeeze', id: 'x', kind: 'axis' }, ''), undefined);
    assert.throws(() => fieldValue({ key: 'squeeze', id: 'x', kind: 'axis' }, '3'), /world axis 0, 1 or 2/);
    assert.deepEqual(fieldValue({ key: 'parameters', id: 'x', kind: 'names' }, 'Load, Modulus,'), ['Load', 'Modulus']);
    assert.throws(() => fieldValue({ key: 'parameters', id: 'x', kind: 'names' }, ' , '), /at least one Metadata key/);
    assert.equal(fieldValue({ key: 'trunk', id: 'x', kind: 'choice' }, 'budget'), 'budget');
});

test('familyArgs posts only the chosen family keys and omits blank optionals', () => {
    const values = {
        't-fno-resolution': '15,15,1',
        't-fno-squeeze': '',
        't-fno-latent': '32',
        't-fno-modes': '12',
        't-fno-layers': '4',
    };
    const args = familyArgs(familyFor('fno'), (id) => values[id] ?? '');
    assert.deepEqual(args, {
        resolution: [15, 15, 1],
        latent_channels: 32,
        num_fno_modes: 12,
        num_fno_layers: 4,
    });
    assert.ok(!('squeeze' in args));
    // afno cannot run without a squeeze: it is refused here, by name
    assert.throws(
        () => familyArgs(familyFor('afno'), (id) => (id === 't-afno-resolution' ? '15,15,1' : '')),
        /squeeze is required/,
    );
    const deep = familyArgs(familyFor('deeponet'), (id) =>
        ({ 't-deep-params': 'Load,Modulus', 't-deep-trunk': 'points', 't-deep-width': '64' })[id] ?? '',
    );
    assert.deepEqual(deep, { parameters: ['Load', 'Modulus'], trunk: 'points', width: 64 });
});
