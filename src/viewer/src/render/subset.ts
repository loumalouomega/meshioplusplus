/**
 * Cell subsets of a PolyData surface.
 *
 * vtk.js has no per-cell visibility, so hiding cells (hide, isolate, a hidden
 * group) and drawing a selection overlay both work the same way: build a new
 * cell list holding only the cells in a mask. The points are shared with the
 * original, so this costs one pass over the cell arrays and no point copy.
 *
 * `cellMap[i]` is the original cell id of the subset's cell `i`; the picker
 * uses it to report the cell the user actually clicked.
 */
export interface CellSubset {
    verts: Uint32Array;
    lines: Uint32Array;
    polys: Uint32Array;
    /** Subset cell index -> original cell index (verts, then lines, then polys). */
    cellMap: Uint32Array;
}

export function subsetCells(
    verts: ArrayLike<number>,
    lines: ArrayLike<number>,
    polys: ArrayLike<number>,
    keep: ArrayLike<number>
): CellSubset {
    let cell = 0;
    const out: Uint32Array[] = [];
    const mapParts: number[] = [];
    for (const array of [verts, lines, polys]) {
        // Two passes so the result is allocated exactly once.
        let size = 0;
        let c = cell;
        for (let i = 0; i < array.length; i += array[i] + 1, c++) if (keep[c]) size += array[i] + 1;
        const part = new Uint32Array(size);
        let at = 0;
        for (let i = 0; i < array.length; i += array[i] + 1, cell++) {
            if (!keep[cell]) continue;
            const n = array[i];
            for (let k = 0; k <= n; k++) part[at++] = array[i + k];
            mapParts.push(cell);
        }
        out.push(part);
    }
    return { verts: out[0], lines: out[1], polys: out[2], cellMap: Uint32Array.from(mapParts) };
}

type TypedArray =
    | Int8Array
    | Uint8Array
    | Int16Array
    | Uint16Array
    | Int32Array
    | Uint32Array
    | Float32Array
    | Float64Array;

/** The tuples `map` names, in order: a per-cell array restricted to a subset. */
export function takeTuples<T extends TypedArray>(
    values: T,
    components: number,
    map: ArrayLike<number>
): T {
    const Ctor = values.constructor as new (length: number) => T;
    const out = new Ctor(map.length * components);
    for (let i = 0; i < map.length; i++) {
        for (let k = 0; k < components; k++) out[i * components + k] = values[map[i] * components + k];
    }
    return out;
}
