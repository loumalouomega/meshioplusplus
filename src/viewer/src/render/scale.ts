/**
 * Colour-scale helpers: log scaling and banding for the legend.
 *
 * The 3D colouring uses vtk.js's own `LOG10` scale and discretisation; these
 * functions keep the DOM legend in step with them, in pure code.
 */

/** A log scale needs a strictly positive range. */
export function canLogScale(min: number, max: number): boolean {
    return Number.isFinite(min) && Number.isFinite(max) && min > 0 && max > min;
}

/** `count` tick values from `min` to `max`, evenly spaced in log10 when `log`. */
export function legendTicks(min: number, max: number, count: number, log: boolean): number[] {
    const out: number[] = [];
    const n = Math.max(2, count);
    for (let i = 0; i < n; i++) {
        const t = i / (n - 1);
        out.push(log && canLogScale(min, max) ? min * Math.pow(max / min, t) : min + t * (max - min));
    }
    return out;
}

/** The band count offered for a select value ("0" = continuous). */
export function bandCount(value: string): number {
    const n = Number(value);
    return Number.isInteger(n) && n >= 2 && n <= 64 ? n : 0;
}

/**
 * Quantise a 0..1 position to the centre of its band, for drawing the banded
 * legend from the same colour stops. `bands <= 0` leaves it unchanged.
 */
export function bandPosition(t: number, bands: number): number {
    if (bands <= 0) return t;
    const x = Math.min(1, Math.max(0, t));
    const i = Math.min(bands - 1, Math.floor(x * bands));
    return (i + 0.5) / bands;
}
