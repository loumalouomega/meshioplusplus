/**
 * The only place in `src/` that casts around a vtk.js typing gap.
 *
 * Methods that merely *exist at runtime but are undeclared* are handled by
 * declaration merging in `vtk-augmentations.d.ts`, which is preferable because
 * the result is still type-checked. This file is for the residue: methods
 * whose declared signature is wrong or absent in a way merging cannot express.
 *
 * (The augmentation file is deliberately *not* named `vtk-shims.d.ts`:
 * TypeScript would resolve it and this file to the same module name, keep the
 * `.ts`, and silently drop the declarations.)
 *
 * Every cast here is named, commented, and narrow. `as any` is banned in
 * `src/` and checked in CI; `as unknown as <explicit shape>` is what that ban
 * pushes you towards, and it still describes what is being called.
 */
import type vtkAxesActor from '@kitware/vtk.js/Rendering/Core/AxesActor';

import type { Vector3 } from './types';

/**
 * `AxesActor` inherits `setPosition`/`setScale` from `vtkProp3D` via
 * `macro.setGetArray`, but `AxesActor.d.ts` re-declares the interface without
 * them.
 */
export function placeAxes(axes: vtkAxesActor, origin: Vector3, scale: number): void {
    const prop = axes as unknown as {
        setPosition(x: number, y: number, z: number): void;
        setScale(x: number, y: number, z: number): void;
    };
    prop.setPosition(...origin);
    prop.setScale(scale, scale, scale);
}

/** The subset of `AxesActor`'s config this viewer sets. */
export interface AxesConfig {
    recenter: boolean;
    xLabel: string;
    yLabel: string;
    zLabel: string;
    tipLength: number;
    tipRadius: number;
    shaftRadius: number;
}

/**
 * `AxesActor.d.ts` types `setConfig` as taking the full internal config
 * object; in practice it merges a partial one.
 */
export function setAxesConfig(axes: vtkAxesActor, config: AxesConfig): void {
    (axes as unknown as { setConfig(c: AxesConfig): void }).setConfig(config);
}
