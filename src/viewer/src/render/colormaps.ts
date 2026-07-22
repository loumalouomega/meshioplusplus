/**
 * Colour maps as ordered RGB stop lists.
 *
 * One source of truth: the same stops build the vtk.js colour transfer
 * function (the 3D colouring) and the CSS gradient (the DOM legend). That is
 * the only way the legend can be *proved* to match the render — and it is
 * tested exactly that way, by comparing {@link colorAt} against the CTF's own
 * `getColor` at the same parameter, with no screenshot involved.
 *
 * Hand-rolled rather than using vtk.js's `ColorMapsLite`: those presets are
 * `RGBPoints` arrays that would have to be re-derived into a CSS gradient
 * anyway, `applyColorMap` clobbers the mapping range, and the preset JSON is
 * ~100 KB in a build with a tight single-file budget.
 *
 * Nine stops for the perceptual maps, not five. A five-stop Viridis is
 * visibly wrong through the teal midsection, and each extra stop costs ~25
 * bytes.
 */
import vtkColorTransferFunction from '@kitware/vtk.js/Rendering/Core/ColorTransferFunction';

import type { Vector3 } from '../types';

/** `[t, r, g, b]` with every channel in 0..1. */
export type ColorStop = readonly [number, number, number, number];

export interface Colormap {
    readonly name: string;
    readonly stops: readonly ColorStop[];
}

const COOL_TO_WARM: readonly ColorStop[] = [
    [0.0, 0.231, 0.298, 0.753],
    [0.5, 0.865, 0.865, 0.865],
    [1.0, 0.706, 0.016, 0.149],
];

const VIRIDIS: readonly ColorStop[] = [
    [0.0, 0.267, 0.005, 0.329],
    [0.125, 0.283, 0.141, 0.458],
    [0.25, 0.254, 0.265, 0.53],
    [0.375, 0.208, 0.372, 0.553],
    [0.5, 0.164, 0.471, 0.558],
    [0.625, 0.128, 0.567, 0.551],
    [0.75, 0.135, 0.659, 0.518],
    [0.875, 0.267, 0.749, 0.441],
    [1.0, 0.993, 0.906, 0.144],
];

const PLASMA: readonly ColorStop[] = [
    [0.0, 0.051, 0.03, 0.528],
    [0.125, 0.255, 0.014, 0.615],
    [0.25, 0.418, 0.001, 0.658],
    [0.375, 0.562, 0.052, 0.641],
    [0.5, 0.692, 0.166, 0.564],
    [0.625, 0.799, 0.281, 0.47],
    [0.75, 0.881, 0.393, 0.383],
    [0.875, 0.949, 0.518, 0.296],
    [1.0, 0.94, 0.975, 0.131],
];

const TURBO: readonly ColorStop[] = [
    [0.0, 0.19, 0.072, 0.232],
    [0.125, 0.246, 0.446, 0.938],
    [0.25, 0.164, 0.71, 0.891],
    [0.375, 0.106, 0.899, 0.715],
    [0.5, 0.394, 0.982, 0.417],
    [0.625, 0.71, 0.94, 0.207],
    [0.75, 0.933, 0.759, 0.184],
    [0.875, 0.983, 0.454, 0.111],
    [1.0, 0.729, 0.104, 0.019],
];

const GRAYSCALE: readonly ColorStop[] = [
    [0.0, 0.0, 0.0, 0.0],
    [1.0, 1.0, 1.0, 1.0],
];

const RAINBOW: readonly ColorStop[] = [
    [0.0, 0.0, 0.0, 1.0],
    [0.25, 0.0, 1.0, 1.0],
    [0.5, 0.0, 1.0, 0.0],
    [0.75, 1.0, 1.0, 0.0],
    [1.0, 1.0, 0.0, 0.0],
];

const MAPS: readonly Colormap[] = [
    { name: 'Cool to Warm', stops: COOL_TO_WARM },
    { name: 'Viridis', stops: VIRIDIS },
    { name: 'Plasma', stops: PLASMA },
    { name: 'Turbo', stops: TURBO },
    { name: 'Grayscale', stops: GRAYSCALE },
    { name: 'Rainbow', stops: RAINBOW },
];

export const COLORMAPS: readonly string[] = MAPS.map((m) => m.name);
export const DEFAULT_COLORMAP = 'Cool to Warm';

function mapByName(name: string): Colormap {
    return MAPS.find((m) => m.name === name) ?? MAPS[0];
}

/** A colour transfer function over `[min, max]`, built from the map's stops. */
export function makeColorTransferFunction(
    name: string,
    min: number,
    max: number
): vtkColorTransferFunction {
    const ctf = vtkColorTransferFunction.newInstance();
    const span = max - min;
    for (const [t, r, g, b] of mapByName(name).stops) {
        // A degenerate range would put every stop at the same x, leaving the
        // CTF with a single control point and no gradient at all.
        ctf.addRGBPoint(span === 0 ? min + t : min + t * span, r, g, b);
    }
    // Non-finite values are excluded from the range; give them a neutral grey
    // rather than whichever end of the colormap they happen to clamp to.
    ctf.setNanColor(0.55, 0.55, 0.55, 1);
    return ctf;
}

/** The colour at `t` in 0..1, by the same linear interpolation the CTF uses. */
export function colorAt(name: string, t: number): Vector3 {
    const stops = mapByName(name).stops;
    const x = Math.min(1, Math.max(0, t));
    for (let i = 1; i < stops.length; i++) {
        const [t1, r1, g1, b1] = stops[i];
        if (x > t1) continue;
        const [t0, r0, g0, b0] = stops[i - 1];
        const span = t1 - t0;
        const f = span === 0 ? 0 : (x - t0) / span;
        return [r0 + f * (r1 - r0), g0 + f * (g1 - g0), b0 + f * (b1 - b0)];
    }
    const [, r, g, b] = stops[stops.length - 1];
    return [r, g, b];
}

/** The same stops as a CSS gradient, for the DOM legend. */
export function cssGradient(name: string, direction = 'to top'): string {
    const stops = mapByName(name).stops.map(([t, r, g, b]) => {
        const c = (v: number) => Math.round(v * 255);
        return `rgb(${c(r)}, ${c(g)}, ${c(b)}) ${(t * 100).toFixed(1)}%`;
    });
    return `linear-gradient(${direction}, ${stops.join(', ')})`;
}

/**
 * Format a scalar for the legend.
 *
 * Switches to exponential outside `[0.01, 1000)`: a fixed-point rendering of
 * 1.2e-9 is "0.000", which reads as zero.
 */
export function formatValue(v: number): string {
    if (!Number.isFinite(v)) return '—';
    const a = Math.abs(v);
    if (a !== 0 && (a >= 1000 || a < 0.01)) return v.toExponential(2);
    return v.toFixed(3);
}
