/**
 * Declaration merging for vtk.js methods that exist at runtime but are absent
 * from the shipped `.d.ts`.
 *
 * These are all real `macro.setGet` accessors — verified present in the `.js`
 * — so this adds no behaviour, it only tells TypeScript what is already there.
 * This file plus `vtk-shims.ts` are the only places vtk.js typing gaps are
 * papered over; `as any` is banned in `src/` and checked in CI.
 *
 * Note that `vtkMapper`'s coincident-topology methods are NOT here: they are
 * already declared in `Rendering/Core/Mapper.d.ts`.
 */

declare module '@kitware/vtk.js/Rendering/Core/ColorTransferFunction' {
    interface vtkColorTransferFunction {
        /** ColorTransferFunction.js:1164 — colour for non-finite values. */
        setNanColor(r: number, g: number, b: number, a?: number): boolean;
        /** ColorTransferFunction.js:545/562 — clamp indicators. */
        setUseAboveRangeColor(v: boolean): boolean;
        setUseBelowRangeColor(v: boolean): boolean;
        setAboveRangeColor(r: number, g: number, b: number, a?: number): boolean;
        setBelowRangeColor(r: number, g: number, b: number, a?: number): boolean;
    }
}

export {};
