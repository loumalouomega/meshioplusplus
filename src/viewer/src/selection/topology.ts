/**
 * Cell topology of the displayed surface, built from the three PolyData cell
 * arrays (verts, lines, polys; VTK legacy layout `[n, id0 … idn-1, …]`).
 *
 * Cell numbering is the one vtk.js and the worker use: verts first, then lines,
 * then polys, each in array order — the order `cellData` values come in and
 * the order a picked cell id means.
 *
 * Pure typed-array code; no vtk.js, so it is unit-tested in Node.
 */
export interface CellTopology {
    numCells: number;
    /** `offsets[c] .. offsets[c+1]` indexes `points`. Length `numCells + 1`. */
    offsets: Uint32Array;
    points: Uint32Array;
    /** How many cells each section holds, in numbering order. */
    sections: { verts: number; lines: number; polys: number };
}

function countCells(a: ArrayLike<number>): number {
    let n = 0;
    for (let i = 0; i < a.length; i += a[i] + 1) n++;
    return n;
}

export function buildCellTopology(
    verts: ArrayLike<number>,
    lines: ArrayLike<number>,
    polys: ArrayLike<number>
): CellTopology {
    const sections = {
        verts: countCells(verts),
        lines: countCells(lines),
        polys: countCells(polys),
    };
    const numCells = sections.verts + sections.lines + sections.polys;
    const total = (a: ArrayLike<number>, n: number) => a.length - n;
    const points = new Uint32Array(
        total(verts, sections.verts) + total(lines, sections.lines) + total(polys, sections.polys)
    );
    const offsets = new Uint32Array(numCells + 1);

    let cell = 0;
    let at = 0;
    for (const array of [verts, lines, polys]) {
        for (let i = 0; i < array.length; ) {
            const n = array[i];
            offsets[cell++] = at;
            for (let k = 1; k <= n; k++) points[at++] = array[i + k];
            i += n + 1;
        }
    }
    offsets[numCells] = at;
    return { numCells, offsets, points, sections };
}

/** For each point, the cells that use it (CSR: `start[p] .. start[p+1]` into `cells`). */
export interface PointCells {
    start: Uint32Array;
    cells: Uint32Array;
}

export function buildPointCells(topo: CellTopology, numPoints: number): PointCells {
    const start = new Uint32Array(numPoints + 1);
    for (let i = 0; i < topo.points.length; i++) start[topo.points[i] + 1]++;
    for (let p = 0; p < numPoints; p++) start[p + 1] += start[p];
    const fill = start.slice(0, numPoints);
    const cells = new Uint32Array(topo.points.length);
    for (let c = 0; c < topo.numCells; c++) {
        for (let i = topo.offsets[c]; i < topo.offsets[c + 1]; i++) {
            cells[fill[topo.points[i]]++] = c;
        }
    }
    return { start, cells };
}

/** Which section a cell index belongs to. */
export function cellSection(topo: CellTopology, cell: number): 'verts' | 'lines' | 'polys' {
    if (cell < topo.sections.verts) return 'verts';
    if (cell < topo.sections.verts + topo.sections.lines) return 'lines';
    return 'polys';
}

/** Cell centroids (mean of the cell's points), `xyz` interleaved. */
export function cellCentroids(coords: ArrayLike<number>, topo: CellTopology): Float64Array {
    const out = new Float64Array(topo.numCells * 3);
    for (let c = 0; c < topo.numCells; c++) {
        const a = topo.offsets[c];
        const b = topo.offsets[c + 1];
        let x = 0;
        let y = 0;
        let z = 0;
        for (let i = a; i < b; i++) {
            const p = topo.points[i] * 3;
            x += coords[p];
            y += coords[p + 1];
            z += coords[p + 2];
        }
        const n = Math.max(1, b - a);
        out[c * 3] = x / n;
        out[c * 3 + 1] = y / n;
        out[c * 3 + 2] = z / n;
    }
    return out;
}

/**
 * Unit normals of polygon cells (Newell's method, robust for a warped quad).
 * Cells with fewer than three points get a zero normal.
 */
export function cellNormals(coords: ArrayLike<number>, topo: CellTopology): Float32Array {
    const out = new Float32Array(topo.numCells * 3);
    for (let c = 0; c < topo.numCells; c++) {
        const a = topo.offsets[c];
        const n = topo.offsets[c + 1] - a;
        if (n < 3) continue;
        let nx = 0;
        let ny = 0;
        let nz = 0;
        for (let i = 0; i < n; i++) {
            const p = topo.points[a + i] * 3;
            const q = topo.points[a + ((i + 1) % n)] * 3;
            const px = coords[p];
            const py = coords[p + 1];
            const pz = coords[p + 2];
            const qx = coords[q];
            const qy = coords[q + 1];
            const qz = coords[q + 2];
            nx += (py - qy) * (pz + qz);
            ny += (pz - qz) * (px + qx);
            nz += (px - qx) * (py + qy);
        }
        const len = Math.hypot(nx, ny, nz);
        if (len > 0) {
            out[c * 3] = nx / len;
            out[c * 3 + 1] = ny / len;
            out[c * 3 + 2] = nz / len;
        }
    }
    return out;
}
