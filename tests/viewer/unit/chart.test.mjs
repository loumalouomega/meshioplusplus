/**
 * The training chart's pure helpers: tick selection on both scales and the
 * compact tick formatting. The drawing itself is asserted end-to-end by the
 * Playwright spec (a real `<svg>` in a real layout).
 *
 * Run with `node --test` via tests/viewer/unit/build-and-run.mjs.
 */
import assert from 'node:assert/strict';
import { test } from 'node:test';

import { formatTick, logTicks, niceTicks } from '../../../src/viewer/src/dataset/chart.ts';

test('niceTicks picks round steps covering the range', () => {
    assert.deepEqual(niceTicks(0, 10, 5), [0, 2, 4, 6, 8, 10]);
    const ticks = niceTicks(3, 97, 5);
    assert.ok(ticks[0] >= 3 && ticks[ticks.length - 1] <= 97);
    assert.ok(ticks.every((t, i) => i === 0 || t > ticks[i - 1]));
    assert.ok(niceTicks(0, 1, 5).length >= 4);
    // degenerate inputs never throw
    assert.deepEqual(niceTicks(5, 5), [5]);
    assert.deepEqual(niceTicks(NaN, 1), []);
});

test('logTicks are decades inside the range, thinned to the budget', () => {
    assert.deepEqual(logTicks(1e-4, 1), [1e-4, 1e-3, 0.01, 0.1, 1]);
    const many = logTicks(1e-9, 1e3, 4);
    assert.ok(many.length <= 5 && many.every((v) => v > 0));
    assert.equal(logTicks(0.002, 0.5).at(0), 1e-3);
    // non-positive ranges are simply empty (the caller falls back to linear)
    assert.deepEqual(logTicks(0, 1), []);
    assert.deepEqual(logTicks(-1, -0.1), []);
});

test('formatTick stays compact across magnitudes', () => {
    assert.equal(formatTick(0), '0');
    assert.equal(formatTick(0.25), '0.25');
    assert.equal(formatTick(12), '12');
    assert.equal(formatTick(0.000123), '1.2e-4');
    assert.equal(formatTick(4.2e-9), '4.2e-9');
    assert.equal(formatTick(1234567), '1.2e6');
    assert.equal(formatTick(NaN), '');
});
