/**
 * The legend must match the render.
 *
 * `colorAt` drives the DOM legend's CSS gradient and `makeColorTransferFunction`
 * drives the 3D colouring. If they ever disagree, the legend lies about what
 * you are looking at — which no amount of screenshot review would reliably
 * catch, and which this test catches exactly.
 *
 * Run with `node --test` (built in; no test framework dependency).
 */
import assert from 'node:assert/strict';
import { test } from 'node:test';

import {
    COLORMAPS,
    colorAt,
    cssGradient,
    formatValue,
    makeColorTransferFunction,
} from '../../../src/viewer/src/render/colormaps.ts';

/** vtk.js quantizes its internal table, so allow a small tolerance. */
const TOL = 0.02;

/** `getColor` writes into an out-parameter rather than returning. */
function ctfColor(ctf, x) {
    const rgb = [0, 0, 0];
    ctf.getColor(x, rgb);
    return rgb;
}

test('colorAt agrees with the colour transfer function at the same t', () => {
    for (const name of COLORMAPS) {
        const ctf = makeColorTransferFunction(name, 0, 1);
        for (let i = 0; i <= 20; i++) {
            const t = i / 20;
            const fromCtf = ctfColor(ctf, t);
            const fromStops = colorAt(name, t);
            for (let c = 0; c < 3; c++) {
                assert.ok(
                    Math.abs(fromCtf[c] - fromStops[c]) < TOL,
                    `${name} at t=${t} channel ${c}: ` +
                        `legend ${fromStops[c].toFixed(3)} vs render ${fromCtf[c].toFixed(3)}`
                );
            }
        }
    }
});

test('the transfer function spans a non-unit range', () => {
    const ctf = makeColorTransferFunction('Viridis', -50, 150);
    // Endpoints must land on the map's own endpoints, not be clamped to 0..1.
    assert.deepEqual(
        ctfColor(ctf, -50).map((v) => v.toFixed(2)),
        colorAt('Viridis', 0).map((v) => v.toFixed(2))
    );
    assert.deepEqual(
        ctfColor(ctf, 150).map((v) => v.toFixed(2)),
        colorAt('Viridis', 1).map((v) => v.toFixed(2))
    );
});

test('a degenerate range still produces a usable transfer function', () => {
    // A constant field has min === max; the CTF must not collapse to nothing.
    const ctf = makeColorTransferFunction('Turbo', 7, 7);
    const color = ctfColor(ctf, 7);
    assert.equal(color.length, 3);
    assert.ok(color.every((v) => Number.isFinite(v)));
});

test('colorAt clamps outside 0..1 instead of extrapolating', () => {
    for (const name of COLORMAPS) {
        assert.deepEqual(colorAt(name, -1), colorAt(name, 0));
        assert.deepEqual(colorAt(name, 2), colorAt(name, 1));
    }
});

test('cssGradient emits one stop per colour stop, in order', () => {
    const css = cssGradient('Cool to Warm');
    assert.match(css, /^linear-gradient\(to top, /);
    assert.equal(css.split('rgb(').length - 1, 3);
    assert.match(css, /0\.0%/);
    assert.match(css, /100\.0%/);
});

test('formatValue switches to exponential outside [0.01, 1000)', () => {
    assert.equal(formatValue(0), '0.000');
    assert.equal(formatValue(1.5), '1.500');
    assert.equal(formatValue(999.5), '999.500');
    assert.equal(formatValue(1000), '1.00e+3');
    assert.equal(formatValue(0.001), '1.00e-3');
    assert.equal(formatValue(-1234), '-1.23e+3');
    assert.equal(formatValue(NaN), '—');
});
