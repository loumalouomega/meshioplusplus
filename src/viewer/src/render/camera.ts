/**
 * Camera poses: the standard views, and saved/shared camera states.
 *
 * Pure (no vtk.js), so the geometry and the validation of pasted JSON are
 * unit-tested in Node. meshio++ meshes are Z-up, so "front" looks along +Y
 * from -Y and "top" looks down -Z from +Z.
 */
import type { Vector3 } from '../types';

export type ViewName = 'front' | 'back' | 'right' | 'left' | 'top' | 'bottom' | 'iso';

export const VIEW_NAMES: readonly ViewName[] = ['front', 'back', 'right', 'left', 'top', 'bottom', 'iso'];

/** Digit shortcuts, in the order the view buttons appear. */
export const VIEW_KEYS: Readonly<Record<string, ViewName>> = {
    '1': 'front',
    '2': 'back',
    '3': 'right',
    '4': 'left',
    '5': 'top',
    '6': 'bottom',
    '0': 'iso',
};

export interface ViewPose {
    /** From the focal point toward the eye (not normalised for `iso`). */
    direction: Vector3;
    up: Vector3;
}

export function viewPose(name: ViewName): ViewPose {
    switch (name) {
        case 'front':
            return { direction: [0, -1, 0], up: [0, 0, 1] };
        case 'back':
            return { direction: [0, 1, 0], up: [0, 0, 1] };
        case 'right':
            return { direction: [1, 0, 0], up: [0, 0, 1] };
        case 'left':
            return { direction: [-1, 0, 0], up: [0, 0, 1] };
        case 'top':
            return { direction: [0, 0, 1], up: [0, 1, 0] };
        case 'bottom':
            return { direction: [0, 0, -1], up: [0, 1, 0] };
        case 'iso':
            // A three-quarter view: an axis-aligned start renders a box as a
            // flat rectangle, which reads as a broken viewer.
            return { direction: [1, -1, 0.6], up: [0, 0, 1] };
    }
}

export interface CameraState {
    position: Vector3;
    focalPoint: Vector3;
    viewUp: Vector3;
    parallel: boolean;
    /** Half the visible height in world units (parallel projection). */
    parallelScale?: number;
    /** Vertical field of view in degrees (perspective). */
    viewAngle?: number;
}

function isVec3(v: unknown): v is Vector3 {
    return Array.isArray(v) && v.length === 3 && v.every((x) => typeof x === 'number' && Number.isFinite(x));
}

/** Whether `v` is a usable camera state. Pasted text is untrusted, so check everything. */
export function isCameraState(v: unknown): v is CameraState {
    if (typeof v !== 'object' || v === null) return false;
    const c = v as Record<string, unknown>;
    if (!isVec3(c.position) || !isVec3(c.focalPoint) || !isVec3(c.viewUp)) return false;
    if (typeof c.parallel !== 'boolean') return false;
    for (const k of ['parallelScale', 'viewAngle'] as const) {
        if (c[k] !== undefined && !(typeof c[k] === 'number' && Number.isFinite(c[k] as number) && (c[k] as number) > 0)) {
            return false;
        }
    }
    const [px, py, pz] = c.position;
    const [fx, fy, fz] = c.focalPoint;
    // The eye must not sit on the focal point, and the up vector must not vanish.
    if (Math.hypot(px - fx, py - fy, pz - fz) === 0) return false;
    return Math.hypot(...c.viewUp) > 0;
}

export function parseCameraJson(text: string): CameraState | null {
    try {
        const v: unknown = JSON.parse(text);
        return isCameraState(v) ? v : null;
    } catch {
        return null;
    }
}

export interface Bookmark {
    name: string;
    camera: CameraState;
}

/** Tolerant: keeps the valid entries of whatever was stored. */
export function parseBookmarks(raw: unknown): Bookmark[] {
    if (!Array.isArray(raw)) return [];
    const out: Bookmark[] = [];
    for (const item of raw) {
        if (typeof item !== 'object' || item === null) continue;
        const b = item as Record<string, unknown>;
        if (typeof b.name === 'string' && b.name.trim() && isCameraState(b.camera)) {
            out.push({ name: b.name.trim(), camera: b.camera });
        }
    }
    return out;
}

/** Add or replace a bookmark by name. */
export function withBookmark(list: readonly Bookmark[], bookmark: Bookmark): Bookmark[] {
    return [...list.filter((b) => b.name !== bookmark.name), bookmark];
}
