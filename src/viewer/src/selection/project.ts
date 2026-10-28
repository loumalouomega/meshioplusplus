/**
 * Screen-space selection: which cells fall inside a rubber band or a lasso.
 *
 * Cells are tested by their centroid projected through the camera, so the work
 * is one matrix multiply per cell in typed arrays and is independent of the
 * renderer. `matrix` is the row-major world → clip matrix, which is what
 * vtk.js's `camera.getCompositeProjectionMatrix` returns.
 */
import type { Mask } from './mask';
import { emptyMask } from './mask';

export interface Projected {
    /** Display x, y in CSS pixels (y down), per cell. */
    x: Float32Array;
    y: Float32Array;
    /** 1 when the centroid is inside the view frustum's depth range. */
    inFront: Uint8Array;
}

export function projectPoints(
    xyz: ArrayLike<number>,
    matrix: ArrayLike<number>,
    width: number,
    height: number
): Projected {
    const n = Math.floor(xyz.length / 3);
    const x = new Float32Array(n);
    const y = new Float32Array(n);
    const inFront = new Uint8Array(n);
    const m = matrix;
    for (let i = 0; i < n; i++) {
        const px = xyz[i * 3];
        const py = xyz[i * 3 + 1];
        const pz = xyz[i * 3 + 2];
        const cx = m[0] * px + m[1] * py + m[2] * pz + m[3];
        const cy = m[4] * px + m[5] * py + m[6] * pz + m[7];
        const cz = m[8] * px + m[9] * py + m[10] * pz + m[11];
        const cw = m[12] * px + m[13] * py + m[14] * pz + m[15];
        if (cw <= 0) continue;
        const nx = cx / cw;
        const ny = cy / cw;
        const nz = cz / cw;
        x[i] = ((nx + 1) / 2) * width;
        y[i] = ((1 - ny) / 2) * height;
        inFront[i] = nz >= -1 && nz <= 1 ? 1 : 0;
    }
    return { x, y, inFront };
}

export interface Rect {
    x0: number;
    y0: number;
    x1: number;
    y1: number;
}

export type Shape = { kind: 'rect'; rect: Rect } | { kind: 'polygon'; points: number[] };

/** Even-odd point-in-polygon; `points` is `[x0, y0, x1, y1, …]`. */
export function pointInPolygon(px: number, py: number, points: ArrayLike<number>): boolean {
    let inside = false;
    const n = Math.floor(points.length / 2);
    for (let i = 0, j = n - 1; i < n; j = i++) {
        const xi = points[i * 2];
        const yi = points[i * 2 + 1];
        const xj = points[j * 2];
        const yj = points[j * 2 + 1];
        if (yi > py !== yj > py && px < ((xj - xi) * (py - yi)) / (yj - yi) + xi) inside = !inside;
    }
    return inside;
}

/**
 * Cells whose projected centroid lies inside `shape`.
 *
 * `eligible` (optional) restricts the test, for example to front-facing or
 * currently visible cells. A cell behind the camera is never selected.
 */
export function cellsInShape(proj: Projected, shape: Shape, eligible?: Mask): Mask {
    const n = proj.x.length;
    const out = emptyMask(n);
    if (shape.kind === 'rect') {
        const x0 = Math.min(shape.rect.x0, shape.rect.x1);
        const x1 = Math.max(shape.rect.x0, shape.rect.x1);
        const y0 = Math.min(shape.rect.y0, shape.rect.y1);
        const y1 = Math.max(shape.rect.y0, shape.rect.y1);
        for (let i = 0; i < n; i++) {
            if (!proj.inFront[i] || (eligible && !eligible[i])) continue;
            if (proj.x[i] >= x0 && proj.x[i] <= x1 && proj.y[i] >= y0 && proj.y[i] <= y1) out[i] = 1;
        }
        return out;
    }
    if (shape.points.length < 6) return out;
    let minX = Infinity;
    let minY = Infinity;
    let maxX = -Infinity;
    let maxY = -Infinity;
    for (let i = 0; i < shape.points.length; i += 2) {
        minX = Math.min(minX, shape.points[i]);
        maxX = Math.max(maxX, shape.points[i]);
        minY = Math.min(minY, shape.points[i + 1]);
        maxY = Math.max(maxY, shape.points[i + 1]);
    }
    for (let i = 0; i < n; i++) {
        if (!proj.inFront[i] || (eligible && !eligible[i])) continue;
        const x = proj.x[i];
        const y = proj.y[i];
        if (x < minX || x > maxX || y < minY || y > maxY) continue;
        if (pointInPolygon(x, y, shape.points)) out[i] = 1;
    }
    return out;
}

/**
 * Cells facing the viewer: the normal points toward the eye. Cells with no
 * normal (lines, vertices) always count, since they have no back side.
 * `eye` is the camera position for a perspective view; pass `direction` (the
 * unit view direction) for a parallel one.
 */
export function facingMask(
    centroids: ArrayLike<number>,
    normals: ArrayLike<number>,
    eye: [number, number, number] | null,
    direction: [number, number, number] | null
): Mask {
    const n = Math.floor(centroids.length / 3);
    const out = emptyMask(n);
    for (let i = 0; i < n; i++) {
        const nx = normals[i * 3];
        const ny = normals[i * 3 + 1];
        const nz = normals[i * 3 + 2];
        if (nx === 0 && ny === 0 && nz === 0) {
            out[i] = 1;
            continue;
        }
        let tx: number;
        let ty: number;
        let tz: number;
        if (eye) {
            tx = eye[0] - centroids[i * 3];
            ty = eye[1] - centroids[i * 3 + 1];
            tz = eye[2] - centroids[i * 3 + 2];
        } else {
            tx = -(direction?.[0] ?? 0);
            ty = -(direction?.[1] ?? 0);
            tz = -(direction?.[2] ?? 0);
        }
        // A face seen exactly edge-on is neither facing nor not; keep it.
        out[i] = nx * tx + ny * ty + nz * tz >= 0 ? 1 : 0;
    }
    return out;
}
