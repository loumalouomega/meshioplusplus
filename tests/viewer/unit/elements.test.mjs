import assert from 'node:assert/strict';
import { test } from 'node:test';

import { elementFromRow, elementKind, relate } from '../../../src/viewer/src/analysis/elements.ts';

const stats = { rms: 0, sigma: 0, formError: 0, numCells: 1, numUsed: 3, numRejected: 0 };
const row = (name, type, params) =>
    elementFromRow(name, type, [...params, ...new Array(10 - params.length).fill(0)], stats);

const close = (a, b, eps = 1e-9) => assert.ok(Math.abs(a - b) < eps, `${a} != ${b}`);

test('patch types map to kinds and an unfitted patch is refused', () => {
    assert.equal(elementKind(1), 'plane');
    assert.equal(elementKind(3), 'cylinder');
    assert.equal(elementKind(0), null);
    assert.throws(() => row('F1', 0, []), /no primitive/);
});

test('a row keeps the radius, the normalised axis and the cone angle in degrees', () => {
    const bore = row('F1', 3, [1, 2, 3, 0, 0, 2, 2.5]);
    assert.deepEqual(bore.direction, [0, 0, 1]);
    assert.equal(bore.radius, 2.5);
    const cone = row('F2', 4, [0, 0, 0, 0, 1, 0, 0, 0, 30]);
    close(cone.angleDeg, 30);
    assert.equal(cone.radius, 0);
});

test('two parallel planes give their gap; crossing planes give no distance', () => {
    const a = row('A', 1, [0, 0, 0, 0, 0, 1]);
    const b = row('B', 1, [5, 5, 3, 0, 0, -1]);
    const r = relate(a, b);
    close(r.angleDeg, 0);
    close(r.distance, 3);
    const c = row('C', 1, [0, 0, 0, 1, 0, 0]);
    const s = relate(a, c);
    close(s.angleDeg, 90);
    assert.equal(s.distance, null);
});

test('skew axes give the common perpendicular; parallel axes their spacing', () => {
    const x = row('X', 3, [0, 0, 0, 1, 0, 0, 1]);
    const y = row('Y', 3, [0, 0, 4, 0, 1, 0, 1]);
    const r = relate(x, y);
    close(r.angleDeg, 90);
    close(r.distance, 4);
    const x2 = row('X2', 3, [7, 3, 4, 2, 0, 0, 1]);
    close(relate(x, x2).distance, 5);
});

test('an axis against a plane: perpendicular has no distance, parallel has one', () => {
    const floor = row('A', 1, [0, 0, 0, 0, 0, 1]);
    const post = row('P', 3, [1, 1, 0, 0, 0, 1, 0.5]);
    const r = relate(floor, post);
    close(r.angleDeg, 90);
    assert.equal(r.distance, null);
    const rail = row('R', 3, [0, 0, 2, 1, 0, 0, 0.5]);
    const s = relate(rail, floor);
    close(s.angleDeg, 0);
    close(s.distance, 2);
    assert.equal(s.what, 'axis to plane');
});

test('a sphere is measured from its centre', () => {
    const ball = row('S', 2, [3, 4, 0, 0, 0, 0, 1]);
    const origin = row('O', 2, [0, 0, 0, 0, 0, 0, 1]);
    const r = relate(ball, origin);
    assert.equal(r.angleDeg, null);
    close(r.distance, 5);
    const z = row('Z', 3, [0, 0, 0, 0, 0, 1, 1]);
    close(relate(z, ball).distance, 5);
});
