/**
 * Region extraction for the Groups panel (worker side).
 *
 * The WASM surface (`surfaceBuffersOps`) carries cell and point arrays but no
 * region membership, so this reads the mesh once more with `readMesh` (which
 * does carry `Mesh.regions`) and maps the regions onto the drawn surface with
 * the pure `groups/regionMap` module. The full read is the cost of the
 * feature: it is requested separately, after the first render, and refused
 * above {@link REGION_CELL_LIMIT} cells, where a JS copy of the mesh next to
 * the MEMFS bytes and the WASM heap would not fit a 32-bit address space.
 *
 * For a solid the exact surface<->file correspondence comes from running
 * `extractSurface` on a copy of the mesh that carries a sentinel point array:
 * the skin compacts its points deterministically, so the sentinel read back
 * from the result is the surface-point -> file-point map. (Measured against
 * `surfaceBuffersOps`: same point and cell counts, identical coordinates.)
 */
import type { CellBlock, Mesh, Region } from '@meshioplusplus/wasm';

import {
    type BlockShape,
    type RawRegion,
    type RegionInfo,
    blockDim,
    cornerCount,
    directIndex,
    hasVolume,
    mapRegions,
} from '../groups/regionMap';

/** Above this many cells in the file, regions are not extracted. */
export const REGION_CELL_LIMIT = 1_500_000;

export interface RegionsMeshModule {
    readMesh(path: string, format?: string): Mesh;
    extractSurface(mesh: Mesh, recordParentIds?: boolean): Mesh;
}

export interface RegionsOutcome {
    ok: boolean;
    reason?: string;
    skipped?: string;
    numSurfaceCells: number;
    groups: RegionInfo[];
}

function blockCount(c: CellBlock): number {
    if ('rowOffsets' in c) return c.rowOffsets.length - 1;
    if ('cellOffsets' in c) return c.cellOffsets.length - 1;
    return c.nodesPerCell > 0 ? c.data.length / c.nodesPerCell : 0;
}

/** The corner nodes of local cell `k` of a rectangular block, or null. */
function cornersOf(c: CellBlock, k: number): number[] | null {
    if (!('nodesPerCell' in c)) return null;
    const n = c.nodesPerCell;
    const want = cornerCount(c.type);
    const take = want > 0 ? Math.min(want, n) : n;
    const out: number[] = [];
    for (let i = 0; i < take; i++) out.push(c.data[k * n + i]);
    return out;
}

const keyOf = (nodes: number[]): string => nodes.sort((a, b) => a - b).join(',');

export function computeRegions(
    m: RegionsMeshModule,
    path: string,
    format: string,
    totalCells: number,
    numSurfaceCells: number
): RegionsOutcome {
    if (totalCells > REGION_CELL_LIMIT) {
        return {
            ok: false,
            skipped: `groups are not read above ${REGION_CELL_LIMIT.toLocaleString('en-US')} cells`,
            numSurfaceCells,
            groups: [],
        };
    }

    const mesh = m.readMesh(path, format || undefined);
    const blocks: BlockShape[] = mesh.cells.map((c) => ({ type: c.type, numCells: blockCount(c) }));
    const regions: RawRegion[] = (mesh.regions ?? []).map((r) => ({
        name: r.name,
        kind: r.kind,
        dim: r.dim,
        tag: r.tag,
        entries: r.entries,
    }));
    const numFilePoints = Math.floor(mesh.points.length / (mesh.dim || 3));

    if (!hasVolume(blocks)) {
        const result = mapRegions({
            mode: 'direct',
            blocks,
            regions,
            numSurfaceCells,
            numFilePoints,
            numSurfacePoints: numFilePoints,
        });
        return { ok: result.ok, reason: result.reason, numSurfaceCells, groups: result.groups };
    }

    // Solid: skin it with a sentinel so the point map and the owner cells are
    // exactly the surface's own.
    const sentinel = new Float64Array(numFilePoints);
    for (let i = 0; i < numFilePoints; i++) sentinel[i] = i;
    const tagged: Mesh = { ...mesh, point_data: { ...(mesh.point_data ?? {}), viewer_orig: sentinel } };
    const skin = m.extractSurface(tagged, true);

    const skinCells = skin.cells.reduce((n, c) => n + blockCount(c), 0);
    if (skinCells !== numSurfaceCells) {
        return {
            ok: false,
            reason: `the skin has ${skinCells} cells but the drawn surface has ${numSurfaceCells}`,
            numSurfaceCells,
            groups: [],
        };
    }

    const parent = new Float64Array(skinCells);
    const parentBlocks = skin.cell_data?.['surface:parent_cell'] ?? [];
    let at = 0;
    for (const block of parentBlocks) for (let i = 0; i < block.length; i++) parent[at++] = Number(block[i]);

    const orig = skin.point_data?.viewer_orig;
    const numSurfacePoints = orig ? orig.length : 0;
    const surfacePointOf = new Int32Array(numFilePoints).fill(-1);
    if (orig) for (let i = 0; i < orig.length; i++) surfacePointOf[Number(orig[i])] = i;

    // Boundary patches (cells of lower dimension than the volume) are matched
    // to skin facets by their corner nodes. Built only when a region needs it.
    const topDim = 3;
    const starts: number[] = [];
    let acc = 0;
    for (const b of blocks) {
        starts.push(acc);
        acc += b.numCells;
    }
    const needsKeys = regions.some(
        (r) =>
            r.kind === 'cell' &&
            Array.from(r.entries).some((g) => {
                let bi = starts.length - 1;
                while (bi > 0 && starts[bi] > g) bi--;
                return blockDim(blocks[bi]?.type ?? '') < topDim;
            })
    );
    // Blocks of lower dimension are also blocks, shown as groups of their own.
    const hasLowerBlocks = blocks.some((b) => b.numCells > 0 && blockDim(b.type) < topDim);

    let facetsByKey: Map<string, number[]> | null = null;
    const lowerDimFacets = (g: number): number[] | undefined => {
        if (!facetsByKey) {
            facetsByKey = new Map();
            let surface = 0;
            for (const block of skin.cells) {
                const n = blockCount(block);
                for (let k = 0; k < n; k++, surface++) {
                    const corners = cornersOf(block, k);
                    if (!corners) continue;
                    const key = keyOf(corners.map((p) => (orig ? Number(orig[p]) : p)));
                    const list = facetsByKey.get(key);
                    if (list) list.push(surface);
                    else facetsByKey.set(key, [surface]);
                }
            }
        }
        let bi = starts.length - 1;
        while (bi > 0 && starts[bi] > g) bi--;
        const block = mesh.cells[bi];
        const corners = block ? cornersOf(block, g - starts[bi]) : null;
        return corners ? facetsByKey.get(keyOf(corners)) : undefined;
    };

    const result = mapRegions({
        mode: 'skin',
        blocks,
        regions,
        numSurfaceCells,
        parent,
        lowerDimFacets: needsKeys || hasLowerBlocks ? lowerDimFacets : undefined,
        surfacePointOf,
        numSurfacePoints,
        numFilePoints,
    });
    return { ok: result.ok, reason: result.reason, numSurfaceCells, groups: result.groups };
}

/** Every buffer a result holds, for the transfer list. */
export function regionTransfer(groups: readonly RegionInfo[]): ArrayBuffer[] {
    const out: ArrayBuffer[] = [];
    for (const g of groups) {
        if (g.surfaceCells) out.push(g.surfaceCells.buffer as ArrayBuffer);
        if (g.surfacePoints) out.push(g.surfacePoints.buffer as ArrayBuffer);
    }
    return out;
}

// --- exporting user groups with the mesh ------------------------------------------- //

export interface GroupsMeshModule extends RegionsMeshModule {
    writeMesh(path: string, mesh: Mesh, format?: string): string[];
}

export interface UserGroup {
    name: string;
    /** Cells of the drawn surface (ascending). */
    surfaceCells: ArrayLike<number>;
}

/**
 * File cells behind surface cells: a solid's facets map to the cells that own
 * them (so a group drawn on a solid's skin becomes those solid cells); a
 * surface mesh maps one to one.
 */
function fileCellsFor(m: RegionsMeshModule, mesh: Mesh, blocks: BlockShape[], groups: UserGroup[]): Int32Array[] {
    if (hasVolume(blocks)) {
        const skin = m.extractSurface(mesh, true);
        const parent: number[] = [];
        for (const block of skin.cell_data?.['surface:parent_cell'] ?? []) {
            for (let i = 0; i < block.length; i++) parent.push(Number(block[i]));
        }
        return groups.map((g) => {
            const set = new Set<number>();
            for (let i = 0; i < g.surfaceCells.length; i++) {
                const p = parent[g.surfaceCells[i]];
                if (p !== undefined) set.add(p);
            }
            return Int32Array.from([...set].sort((a, b) => a - b));
        });
    }
    const total = blocks.reduce((n, b) => n + b.numCells, 0);
    const direct = directIndex(blocks, total);
    const inverse = new Int32Array(total).fill(-1);
    if (direct) for (let g = 0; g < direct.length; g++) if (direct[g] >= 0) inverse[direct[g]] = g;
    return groups.map((g) => {
        const out: number[] = [];
        for (let i = 0; i < g.surfaceCells.length; i++) {
            const f = inverse[g.surfaceCells[i]];
            if (f >= 0) out.push(f);
        }
        return Int32Array.from(out.sort((a, b) => a - b));
    });
}

/**
 * Write the staged file in `outFormat` with `groups` added as cell regions,
 * then read the result back to report which groups the format actually kept
 * (not every format stores named regions).
 */
export function writeWithGroups(
    m: GroupsMeshModule,
    path: string,
    format: string,
    outPath: string,
    outFormat: string,
    groups: UserGroup[]
): { requested: number; kept: string[] } {
    const mesh = m.readMesh(path, format || undefined);
    const blocks: BlockShape[] = mesh.cells.map((c) => ({ type: c.type, numCells: blockCount(c) }));
    const cells = fileCellsFor(m, mesh, blocks, groups);
    const dim = Math.max(0, ...blocks.filter((b) => b.numCells > 0).map((b) => blockDim(b.type)));

    const taken = new Set((mesh.regions ?? []).map((r) => r.name));
    const added: Region[] = [];
    groups.forEach((g, i) => {
        if (cells[i].length === 0) return;
        let name = g.name;
        for (let n = 2; taken.has(name); n++) name = `${g.name}_${n}`;
        taken.add(name);
        added.push({ name, kind: 'cell', dim, tag: -1, entries: cells[i] });
    });
    if (added.length === 0) {
        m.writeMesh(outPath, mesh, outFormat);
        return { requested: groups.length, kept: [] };
    }

    m.writeMesh(outPath, { ...mesh, regions: [...(mesh.regions ?? []), ...added] }, outFormat);
    let kept: string[] = [];
    try {
        const back = m.readMesh(outPath, outFormat);
        const present = new Set((back.regions ?? []).map((r) => r.name));
        kept = added.map((r) => r.name).filter((n) => present.has(n));
    } catch {
        // A format that cannot be read back here: report nothing rather than guess.
        kept = [];
    }
    return { requested: added.length, kept };
}
