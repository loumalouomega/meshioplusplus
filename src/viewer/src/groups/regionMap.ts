/**
 * Mapping a mesh file's named regions and cell blocks onto the surface the
 * viewer draws.
 *
 * Regions are defined on the *file's* cells (global, block-major indices), but
 * the viewer shows a surface: the skin of a solid, or a surface mesh as it is.
 * Three situations (see the plan notes on `extract_surface`):
 *
 * - **skin** (the file has volume cells): a facet belongs to the region that
 *   holds its owner cell (`surface:parent_cell`). Regions of lower-dimension
 *   cells (boundary patches such as a gmsh physical surface) are not in the
 *   skin at all; they match facets by corner-node key, which the caller
 *   provides through `lowerDimFacets`.
 * - **direct** (no volume cells): the surface is the mesh, numbered in VTK
 *   section order (vertices, lines, then polygons, each in block order).
 * - `side` regions (cell, local facet) are reported but not mapped.
 *
 * Pure typed-array code with no WASM and no DOM, so it is unit-tested in Node
 * and runs in the worker.
 */

export interface BlockShape {
    /** meshio++ cell type name. */
    type: string;
    numCells: number;
}

export interface RawRegion {
    name: string;
    kind: 'point' | 'cell' | 'side';
    dim: number;
    tag: number;
    /** `cell`/`point`: one index per entry; `side`: pairs. Ascending. */
    entries: ArrayLike<number>;
}

export type GroupSource = 'region' | 'block';

export interface RegionInfo {
    name: string;
    kind: 'point' | 'cell' | 'side';
    source: GroupSource;
    dim: number;
    tag: number;
    /** Entries the file holds (cells, points or facets). */
    numEntries: number;
    /** Entries that reached at least one surface cell / point. */
    mapped: number;
    /** Surface cells of the group (ascending), or null when none apply. */
    surfaceCells: Uint32Array | null;
    /** Surface points of a point region (ascending), or null. */
    surfacePoints: Uint32Array | null;
    /** Why nothing was mapped, when that is the case. */
    note?: string;
}

export interface MapInput {
    mode: 'skin' | 'direct';
    blocks: BlockShape[];
    regions: RawRegion[];
    numSurfaceCells: number;
    /** skin: the owner cell (global, file numbering) of each surface cell. */
    parent?: ArrayLike<number>;
    /** skin: surface cells whose corner nodes equal those of lower-dimension cell `g`. */
    lowerDimFacets?: (globalCell: number) => ArrayLike<number> | undefined;
    /** skin: file point -> surface point, -1 where the skin does not use it. */
    surfacePointOf?: ArrayLike<number>;
    numSurfacePoints?: number;
    numFilePoints?: number;
}

export interface MapResult {
    ok: boolean;
    reason?: string;
    /** Named regions of the file, then one entry per cell block. */
    groups: RegionInfo[];
}

// --- cell types ------------------------------------------------------------ //

export function blockDim(type: string): number {
    if (type === 'vertex') return 0;
    if (type.startsWith('line')) return 1;
    if (type.startsWith('triangle') || type.startsWith('quad') || type.startsWith('polygon')) return 2;
    return 3;
}

export function hasVolume(blocks: readonly BlockShape[]): boolean {
    return blocks.some((b) => b.numCells > 0 && blockDim(b.type) === 3);
}

/** Corner count used to compare facets across element orders (triangle6 -> 3). */
export function cornerCount(type: string): number {
    if (type.startsWith('triangle')) return 3;
    if (type.startsWith('quad')) return 4;
    if (type.startsWith('line')) return 2;
    return 0; // all nodes
}

// --- the mapping ------------------------------------------------------------- //

interface BlockIndex {
    starts: number[];
    dims: number[];
    total: number;
}

function indexBlocks(blocks: readonly BlockShape[]): BlockIndex {
    const starts: number[] = [];
    const dims: number[] = [];
    let at = 0;
    for (const b of blocks) {
        starts.push(at);
        dims.push(blockDim(b.type));
        at += b.numCells;
    }
    return { starts, dims, total: at };
}

/** Which block holds global cell `g` (binary search). */
function blockOf(index: BlockIndex, g: number): number {
    let lo = 0;
    let hi = index.starts.length - 1;
    while (lo < hi) {
        const mid = (lo + hi + 1) >> 1;
        if (index.starts[mid] <= g) lo = mid;
        else hi = mid - 1;
    }
    return lo;
}

/**
 * For a surface mesh (no volume cells): file cell -> surface cell, following
 * VTK section order. Returns null if the counts do not add up to the surface
 * the viewer holds (for example a polygon block the worker fanned out).
 */
export function directIndex(blocks: readonly BlockShape[], numSurfaceCells: number): Int32Array | null {
    const index = indexBlocks(blocks);
    const out = new Int32Array(index.total).fill(-1);
    const sections: number[][] = [[], [], []]; // vertices, lines, polygons
    blocks.forEach((b, bi) => {
        const dim = index.dims[bi];
        const section = dim === 0 ? 0 : dim === 1 ? 1 : 2;
        for (let c = 0; c < b.numCells; c++) sections[section].push(index.starts[bi] + c);
    });
    let surface = 0;
    for (const section of sections) for (const g of section) out[g] = surface++;
    return surface === numSurfaceCells ? out : null;
}

/** parent -> surface cells, as CSR. */
function parentIndex(parent: ArrayLike<number>, numCells: number): { start: Uint32Array; cells: Uint32Array } {
    const start = new Uint32Array(numCells + 1);
    for (let i = 0; i < parent.length; i++) {
        const p = parent[i];
        if (p >= 0 && p < numCells) start[p + 1]++;
    }
    for (let p = 0; p < numCells; p++) start[p + 1] += start[p];
    const fill = start.slice(0, numCells);
    const cells = new Uint32Array(start[numCells]);
    for (let i = 0; i < parent.length; i++) {
        const p = parent[i];
        if (p >= 0 && p < numCells) cells[fill[p]++] = i;
    }
    return { start, cells };
}

function sortedUnique(values: number[]): Uint32Array {
    const out = Uint32Array.from(values);
    out.sort();
    let n = 0;
    for (let i = 0; i < out.length; i++) if (i === 0 || out[i] !== out[i - 1]) out[n++] = out[i];
    return out.slice(0, n);
}

export function mapRegions(input: MapInput): MapResult {
    const index = indexBlocks(input.blocks);
    const topDim = Math.max(0, ...input.blocks.filter((b) => b.numCells > 0).map((b) => blockDim(b.type)));

    let direct: Int32Array | null = null;
    let owned: { start: Uint32Array; cells: Uint32Array } | null = null;
    if (input.mode === 'direct') {
        direct = directIndex(input.blocks, input.numSurfaceCells);
        if (!direct) {
            return {
                ok: false,
                reason: 'the file\'s cells do not line up with the surface that is drawn',
                groups: [],
            };
        }
    } else {
        if (!input.parent || input.parent.length !== input.numSurfaceCells) {
            return { ok: false, reason: 'the surface has no owner-cell provenance', groups: [] };
        }
        owned = parentIndex(input.parent, index.total);
    }

    /** Surface cells for one file cell, pushed into `into`. Returns whether any. */
    const surfaceCellsOf = (g: number, into: number[]): boolean => {
        if (g < 0 || g >= index.total) return false;
        if (direct) {
            const s = direct[g];
            if (s < 0) return false;
            into.push(s);
            return true;
        }
        const dim = index.dims[blockOf(index, g)];
        if (dim === topDim) {
            const a = owned!.start[g];
            const b = owned!.start[g + 1];
            for (let k = a; k < b; k++) into.push(owned!.cells[k]);
            return b > a;
        }
        const matched = input.lowerDimFacets?.(g);
        if (!matched || matched.length === 0) return false;
        for (let k = 0; k < matched.length; k++) into.push(matched[k]);
        return true;
    };

    const fromCells = (
        name: string,
        source: GroupSource,
        dim: number,
        tag: number,
        cells: Iterable<number>,
        count: number
    ): RegionInfo => {
        const surface: number[] = [];
        let mapped = 0;
        for (const g of cells) if (surfaceCellsOf(g, surface)) mapped++;
        const unique = sortedUnique(surface);
        return {
            name,
            kind: 'cell',
            source,
            dim,
            tag,
            numEntries: count,
            mapped,
            surfaceCells: unique.length ? unique : null,
            surfacePoints: null,
            note: unique.length ? undefined : 'none of its cells reach the drawn surface',
        };
    };

    const groups: RegionInfo[] = [];

    for (const r of input.regions) {
        if (r.kind === 'cell') {
            groups.push(fromCells(r.name, 'region', r.dim, r.tag, Array.from(r.entries), r.entries.length));
        } else if (r.kind === 'point') {
            const surface: number[] = [];
            let mapped = 0;
            const map = input.surfacePointOf;
            const identity =
                !map && input.numFilePoints !== undefined && input.numFilePoints === input.numSurfacePoints;
            for (let i = 0; i < r.entries.length; i++) {
                const p = identity ? r.entries[i] : map ? map[r.entries[i]] : -1;
                if (p !== undefined && p >= 0) {
                    surface.push(p);
                    mapped++;
                }
            }
            const unique = sortedUnique(surface);
            groups.push({
                name: r.name,
                kind: 'point',
                source: 'region',
                dim: r.dim,
                tag: r.tag,
                numEntries: r.entries.length,
                mapped,
                surfaceCells: null,
                surfacePoints: unique.length ? unique : null,
                note: unique.length ? undefined : 'none of its points are on the drawn surface',
            });
        } else {
            groups.push({
                name: r.name,
                kind: 'side',
                source: 'region',
                dim: r.dim,
                tag: r.tag,
                numEntries: r.entries.length / 2,
                mapped: 0,
                surfaceCells: null,
                surfacePoints: null,
                note: 'facet sets cannot be shown on the surface yet',
            });
        }
    }

    // One group per cell block, so an Exodus/Abaqus mesh with no regions still
    // has something to show, and a multi-material mesh can be told apart.
    input.blocks.forEach((b, bi) => {
        if (b.numCells === 0) return;
        const start = index.starts[bi];
        const cells = Array.from({ length: b.numCells }, (_, k) => start + k);
        groups.push(fromCells(b.type, 'block', blockDim(b.type), -1, cells, b.numCells));
    });

    return { ok: true, groups };
}

// --- categorical group labels ------------------------------------------------------ //

/**
 * One integer per surface cell: the index (into `groups`) of the first group
 * in `order` that holds it, or -1. Used to colour the surface by group; pass
 * the groups smallest first so a small group (a boundary patch) shows on top
 * of the larger one (the zone) it lies inside.
 */
export function groupLabels(groups: readonly RegionInfo[], numCells: number, order: readonly number[]): Int32Array {
    const out = new Int32Array(numCells).fill(-1);
    for (const gi of order) {
        const cells = groups[gi]?.surfaceCells;
        if (!cells) continue;
        for (let i = 0; i < cells.length; i++) if (out[cells[i]] < 0) out[cells[i]] = gi;
    }
    return out;
}

/** Group indices ordered smallest first (stable for equal sizes). */
export function smallestFirst(groups: readonly RegionInfo[]): number[] {
    return groups
        .map((g, i) => ({ i, n: g.surfaceCells?.length ?? 0 }))
        .sort((a, b) => a.n - b.n || a.i - b.i)
        .map((x) => x.i);
}

/** A group as a 0/1 mask over the surface cells. */
export function groupMask(group: RegionInfo, numCells: number): Uint8Array {
    const mask = new Uint8Array(numCells);
    if (group.surfaceCells) for (let i = 0; i < group.surfaceCells.length; i++) mask[group.surfaceCells[i]] = 1;
    return mask;
}
