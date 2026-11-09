/**
 * Fitted elements for the measuring mode: a primitive fitted from a clicked
 * seed (`fitPrimitives` with `seeds`, roadmap 7.1.1) and the angle and
 * distance between two of them.
 *
 * Pure, so the arithmetic is unit-tested without a browser; the panel in
 * `features/measuring.ts` only formats what this returns.
 */
import type { Vector3 } from '../types';

export type ElementKind = 'plane' | 'sphere' | 'cylinder' | 'cone' | 'torus';

/** `fit:patch_type` codes, 1-based as the core writes them. */
const KINDS: readonly ElementKind[] = ['plane', 'sphere', 'cylinder', 'cone', 'torus'];

/** One fitted feature, from one row of the `fit:patch_*` field data. */
export interface FittedElement {
    /** `F1`, `F2`, ... in the order they were fitted. */
    name: string;
    kind: ElementKind;
    /** Plane: a point on it; sphere and torus: the centre; cylinder: a point on the axis; cone: the apex. */
    point: Vector3;
    /** Plane: the unit normal; cylinder, cone, torus: the unit axis; sphere: zero. */
    direction: Vector3;
    /** Sphere, cylinder: the radius; torus: the major radius; else 0. */
    radius: number;
    /** Torus: the minor radius; else 0. */
    radius2: number;
    /** Cone: the half-angle in degrees; else 0. */
    angleDeg: number;
    /** Root-mean-square residual and the residuals' standard deviation. */
    rms: number;
    sigma: number;
    /** Peak-to-peak form error over the used points (flatness, cylindricity, ...). */
    formError: number;
    /** Cells of the grown feature, and points used and rejected by the trimming. */
    numCells: number;
    numUsed: number;
    numRejected: number;
}

/** The scalars of one patch row that are not in `fit:patch_params`. */
export interface ElementStats {
    rms: number;
    sigma: number;
    formError: number;
    numCells: number;
    numUsed: number;
    numRejected: number;
}

/** The kind of a `fit:patch_type` code, or null for an unfitted patch (0) or an unknown one. */
export function elementKind(type: number): ElementKind | null {
    return KINDS[type - 1] ?? null;
}

function unit(v: Vector3): Vector3 {
    const n = Math.hypot(v[0], v[1], v[2]);
    return n > 0 ? [v[0] / n, v[1] / n, v[2] / n] : [0, 0, 0];
}

/**
 * The element of one patch row: `type` its `fit:patch_type` and `params` its
 * ten `fit:patch_params` (x, y, z, ax, ay, az, r1, r2, the cone's half-angle in degrees, 0).
 * @throws Error for an unfitted patch or fewer than nine parameters.
 */
export function elementFromRow(
    name: string,
    type: number,
    params: ArrayLike<number>,
    stats: ElementStats
): FittedElement {
    const kind = elementKind(type);
    if (!kind) throw new Error(`no primitive explains the feature (patch type ${type})`);
    if (params.length < 9) throw new Error('a fitted primitive has ten parameters');
    const point: Vector3 = [params[0], params[1], params[2]];
    const direction = kind === 'sphere' ? ([0, 0, 0] as Vector3) : unit([params[3], params[4], params[5]]);
    return {
        name,
        kind,
        point,
        direction,
        radius: kind === 'plane' || kind === 'cone' ? 0 : params[6],
        radius2: kind === 'torus' ? params[7] : 0,
        angleDeg: kind === 'cone' ? params[8] : 0,
        ...stats,
    };
}

function dot(a: Vector3, b: Vector3): number {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

function sub(a: Vector3, b: Vector3): Vector3 {
    return [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
}

function cross(a: Vector3, b: Vector3): Vector3 {
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
}

function norm(a: Vector3): number {
    return Math.hypot(a[0], a[1], a[2]);
}

/** Distance from `p` to the line through `o` along the unit `d`. */
function pointLine(p: Vector3, o: Vector3, d: Vector3): number {
    return norm(cross(sub(p, o), d));
}

/** Distance from `p` to the plane through `o` with unit normal `n`. */
function pointPlane(p: Vector3, o: Vector3, n: Vector3): number {
    return Math.abs(dot(sub(p, o), n));
}

/** The angle between two undirected lines, folded into [0, 90] degrees. */
function lineAngle(a: Vector3, b: Vector3): number {
    const c = Math.min(1, Math.abs(dot(a, b)));
    return (Math.acos(c) * 180) / Math.PI;
}

/** Angles below this many degrees count as parallel for a distance. */
export const PARALLEL_DEG = 0.5;

/** How two elements relate: the angle between them and the distance, where defined. */
export interface Relation {
    /**
     * Degrees in [0, 90]: between two planes, between two axes, or between a
     * plane and an axis (0 when the axis lies parallel to the plane). Null when
     * an element has no direction (a sphere).
     */
    angleDeg: number | null;
    /**
     * Two parallel planes: their gap; two axes: the shortest distance between
     * the lines; a plane and a parallel axis: the axis from the plane; a sphere:
     * its centre from the other element. Null when the elements meet at an angle.
     */
    distance: number | null;
    /** What `distance` measures, for the panel ("plane to plane", ...). */
    what: string;
}

type Shape = 'plane' | 'line' | 'point';

function shapeOf(e: FittedElement): Shape {
    if (e.kind === 'plane') return 'plane';
    if (e.kind === 'sphere') return 'point';
    return 'line';
}

/** The angle and distance between `a` and `b` (see {@link Relation}). */
export function relate(a: FittedElement, b: FittedElement): Relation {
    const sa = shapeOf(a);
    const sb = shapeOf(b);
    if (sa === 'point' || sb === 'point') {
        const [p, q, sq] = sa === 'point' ? [a, b, sb] : [b, a, sa];
        if (sq === 'point') return { angleDeg: null, distance: norm(sub(p.point, q.point)), what: 'centre to centre' };
        if (sq === 'plane')
            return { angleDeg: null, distance: pointPlane(p.point, q.point, q.direction), what: 'centre to plane' };
        return { angleDeg: null, distance: pointLine(p.point, q.point, q.direction), what: 'centre to axis' };
    }
    if (sa === 'plane' && sb === 'plane') {
        const angle = lineAngle(a.direction, b.direction);
        return {
            angleDeg: angle,
            distance: angle < PARALLEL_DEG ? pointPlane(b.point, a.point, a.direction) : null,
            what: 'plane to plane',
        };
    }
    if (sa === 'line' && sb === 'line') {
        const angle = lineAngle(a.direction, b.direction);
        const c = cross(a.direction, b.direction);
        const distance =
            angle < PARALLEL_DEG
                ? pointLine(b.point, a.point, a.direction)
                : Math.abs(dot(sub(b.point, a.point), c)) / norm(c);
        return { angleDeg: angle, distance, what: 'axis to axis' };
    }
    const [plane, line] = sa === 'plane' ? [a, b] : [b, a];
    const angle = 90 - lineAngle(plane.direction, line.direction);
    return {
        angleDeg: angle,
        distance: angle < PARALLEL_DEG ? pointPlane(line.point, plane.point, plane.direction) : null,
        what: 'axis to plane',
    };
}
