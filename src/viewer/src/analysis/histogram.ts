/**
 * Histogram and summary statistics of one component of a per-element array.
 * Pure; the Statistics panel draws it and turns a clicked bin into a selection.
 */

export interface Histogram {
    counts: Uint32Array;
    min: number;
    max: number;
    mean: number;
    std: number;
    /** Finite values counted. */
    count: number;
    /** Values that were NaN or infinite and left out. */
    nonFinite: number;
}

export interface HistogramOptions {
    /** `-1` (default) is the magnitude of a multi-component array. */
    component?: number;
    bins?: number;
    /** Fix the range; the default is the data range. */
    range?: [number, number];
}

function valueAt(values: ArrayLike<number>, components: number, component: number, i: number): number {
    if (component < 0 && components > 1) {
        let s = 0;
        for (let k = 0; k < components; k++) s += Number(values[i * components + k]) ** 2;
        return Math.sqrt(s);
    }
    return Number(values[i * components + Math.max(0, component)]);
}

export function histogram(values: ArrayLike<number>, components: number, opts: HistogramOptions = {}): Histogram {
    const bins = Math.max(1, Math.floor(opts.bins ?? 20));
    const component = opts.component ?? -1;
    const n = Math.floor(values.length / components);

    let min = Infinity;
    let max = -Infinity;
    let sum = 0;
    let count = 0;
    let nonFinite = 0;
    for (let i = 0; i < n; i++) {
        const v = valueAt(values, components, component, i);
        if (!Number.isFinite(v)) {
            nonFinite++;
            continue;
        }
        count++;
        sum += v;
        if (v < min) min = v;
        if (v > max) max = v;
    }
    const counts = new Uint32Array(bins);
    if (count === 0) return { counts, min: NaN, max: NaN, mean: NaN, std: NaN, count, nonFinite };

    const mean = sum / count;
    let sq = 0;
    const [lo, hi] = opts.range ?? [min, max];
    const span = hi - lo;
    for (let i = 0; i < n; i++) {
        const v = valueAt(values, components, component, i);
        if (!Number.isFinite(v)) continue;
        sq += (v - mean) ** 2;
        if (v < lo || v > hi) continue;
        // A constant array has no span: every value goes in the first bin.
        const b = span > 0 ? Math.min(bins - 1, Math.floor(((v - lo) / span) * bins)) : 0;
        counts[b]++;
    }
    return { counts, min, max, mean, std: Math.sqrt(sq / count), count, nonFinite };
}

/** The `[from, to]` value range covered by bin `i`. */
export function binRange(h: Histogram, i: number, range?: [number, number]): [number, number] {
    const [lo, hi] = range ?? [h.min, h.max];
    const w = (hi - lo) / h.counts.length;
    return [lo + i * w, i === h.counts.length - 1 ? hi : lo + (i + 1) * w];
}
