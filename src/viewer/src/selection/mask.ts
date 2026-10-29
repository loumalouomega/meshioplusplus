/**
 * Selections as bit-per-cell masks (`Uint8Array`, 1 = selected).
 *
 * A mask is the whole selection: counting, combining, inverting, growing and
 * exporting are plain loops, and a million-cell selection costs a megabyte.
 */
import type { CellTopology, PointCells } from './topology';

export type Mask = Uint8Array;

export function emptyMask(n: number): Mask {
    return new Uint8Array(n);
}

export function countMask(mask: Mask): number {
    let n = 0;
    for (let i = 0; i < mask.length; i++) n += mask[i];
    return n;
}

export function maskFromIds(n: number, ids: ArrayLike<number>): Mask {
    const mask = emptyMask(n);
    for (let i = 0; i < ids.length; i++) {
        const id = ids[i];
        if (id >= 0 && id < n) mask[id] = 1;
    }
    return mask;
}

/** Selected cell ids, ascending. */
export function idsOf(mask: Mask): Uint32Array {
    const out = new Uint32Array(countMask(mask));
    let k = 0;
    for (let i = 0; i < mask.length; i++) if (mask[i]) out[k++] = i;
    return out;
}

export type Combine = 'replace' | 'add' | 'subtract' | 'toggle';

/** Fold a gesture's result into the current selection. */
export function combine(current: Mask, incoming: Mask, how: Combine): Mask {
    if (how === 'replace') return incoming.slice();
    const out = current.slice();
    for (let i = 0; i < out.length; i++) {
        const b = incoming[i];
        if (how === 'add') out[i] = out[i] | b;
        else if (how === 'subtract') out[i] = out[i] & (b ^ 1);
        else out[i] = out[i] ^ b;
    }
    return out;
}

export function invertMask(mask: Mask): Mask {
    const out = new Uint8Array(mask.length);
    for (let i = 0; i < mask.length; i++) out[i] = mask[i] ^ 1;
    return out;
}

/** Every cell that shares a point with a selected cell (one ring outward). */
export function growMask(mask: Mask, topo: CellTopology, pc: PointCells): Mask {
    const out = mask.slice();
    for (let c = 0; c < topo.numCells; c++) {
        if (!mask[c]) continue;
        for (let i = topo.offsets[c]; i < topo.offsets[c + 1]; i++) {
            const p = topo.points[i];
            for (let k = pc.start[p]; k < pc.start[p + 1]; k++) out[pc.cells[k]] = 1;
        }
    }
    return out;
}

/** Remove the outermost ring: keep a cell only if all its neighbours are selected. */
export function shrinkMask(mask: Mask, topo: CellTopology, pc: PointCells): Mask {
    const out = mask.slice();
    for (let c = 0; c < topo.numCells; c++) {
        if (!mask[c]) continue;
        search: for (let i = topo.offsets[c]; i < topo.offsets[c + 1]; i++) {
            const p = topo.points[i];
            for (let k = pc.start[p]; k < pc.start[p + 1]; k++) {
                if (!mask[pc.cells[k]]) {
                    out[c] = 0;
                    break search;
                }
            }
        }
    }
    return out;
}

/** Flood-fill from the selected cells across shared points (the connected piece). */
export function connectedMask(seed: Mask, topo: CellTopology, pc: PointCells): Mask {
    const out = seed.slice();
    const queue: number[] = [];
    for (let c = 0; c < out.length; c++) if (out[c]) queue.push(c);
    while (queue.length) {
        const c = queue.pop() as number;
        for (let i = topo.offsets[c]; i < topo.offsets[c + 1]; i++) {
            const p = topo.points[i];
            for (let k = pc.start[p]; k < pc.start[p + 1]; k++) {
                const d = pc.cells[k];
                if (!out[d]) {
                    out[d] = 1;
                    queue.push(d);
                }
            }
        }
    }
    return out;
}

export interface MaskStats {
    count: number;
    min: number;
    max: number;
    mean: number;
}

/**
 * Statistics of one component of a per-cell array over a selection.
 * `component < 0` means the magnitude of a multi-component array. Non-finite
 * values are skipped, and `count` is how many contributed.
 */
export function maskStats(
    values: ArrayLike<number>,
    components: number,
    component: number,
    mask: Mask
): MaskStats {
    let count = 0;
    let min = Infinity;
    let max = -Infinity;
    let sum = 0;
    for (let c = 0; c < mask.length; c++) {
        if (!mask[c]) continue;
        let v: number;
        if (component < 0 && components > 1) {
            let s = 0;
            for (let k = 0; k < components; k++) s += values[c * components + k] ** 2;
            v = Math.sqrt(s);
        } else {
            v = values[c * components + Math.max(0, component)];
        }
        if (!Number.isFinite(v)) continue;
        count++;
        sum += v;
        if (v < min) min = v;
        if (v > max) max = v;
    }
    return count ? { count, min, max, mean: sum / count } : { count: 0, min: NaN, max: NaN, mean: NaN };
}

/**
 * Cells whose value of a per-cell scalar lies in `[lo, hi]`.
 * A non-finite value never matches, like the extension's threshold.
 */
export function thresholdCellValues(
    values: ArrayLike<number>,
    components: number,
    component: number,
    lo: number,
    hi: number,
    numCells: number
): Mask {
    const [a, b] = lo <= hi ? [lo, hi] : [hi, lo];
    const out = emptyMask(numCells);
    for (let c = 0; c < numCells; c++) {
        let v: number;
        if (component < 0 && components > 1) {
            let s = 0;
            for (let k = 0; k < components; k++) s += values[c * components + k] ** 2;
            v = Math.sqrt(s);
        } else {
            v = values[c * components + Math.max(0, component)];
        }
        if (Number.isFinite(v) && v >= a && v <= b) out[c] = 1;
    }
    return out;
}

/**
 * Cells whose value of a per-POINT scalar lies in `[lo, hi]` at all (`any`) or
 * at every corner (`all`). Mirrors the extension's Nodal threshold rule.
 */
export function thresholdPointValues(
    values: ArrayLike<number>,
    components: number,
    component: number,
    lo: number,
    hi: number,
    topo: CellTopology,
    rule: 'all' | 'any' = 'all'
): Mask {
    const [a, b] = lo <= hi ? [lo, hi] : [hi, lo];
    const out = emptyMask(topo.numCells);
    const inside = (p: number): boolean => {
        let v: number;
        if (component < 0 && components > 1) {
            let s = 0;
            for (let k = 0; k < components; k++) s += values[p * components + k] ** 2;
            v = Math.sqrt(s);
        } else {
            v = values[p * components + Math.max(0, component)];
        }
        return Number.isFinite(v) && v >= a && v <= b;
    };
    for (let c = 0; c < topo.numCells; c++) {
        let hits = 0;
        const n = topo.offsets[c + 1] - topo.offsets[c];
        for (let i = topo.offsets[c]; i < topo.offsets[c + 1]; i++) if (inside(topo.points[i])) hits++;
        if (n > 0 && (rule === 'any' ? hits > 0 : hits === n)) out[c] = 1;
    }
    return out;
}
