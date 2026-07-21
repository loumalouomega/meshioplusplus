/**
 * Saving the current view as a PNG.
 *
 * Uses `captureNextImage`, which schedules the grab against vtk.js's own
 * render loop. Reading the canvas directly usually yields a black image: the
 * WebGL back buffer is not preserved after a frame is presented, so
 * `toDataURL` sees an empty surface unless it happens to run inside the same
 * tick as the draw. That is the classic black-screenshot bug, and it is why
 * this is not a one-liner.
 */
import type { Renderer } from './renderer';

/** Capture the current view. Returns a `data:image/png;base64,...` URL. */
export async function capturePng(renderer: Renderer): Promise<string> {
    const view = renderer.view;
    const pending = view.captureNextImage('image/png');
    // The capture only resolves once another frame is drawn.
    renderer.render();
    try {
        return await pending;
    } catch {
        // Fall back to the raw canvas. Correct whenever the back buffer does
        // happen to be preserved, and better than failing outright.
        const canvas = view.getCanvas();
        if (!canvas) throw new Error('the renderer has no canvas to capture');
        return canvas.toDataURL('image/png');
    }
}

/** Capture and hand the browser a download. */
export async function downloadPng(renderer: Renderer, filename: string): Promise<void> {
    const url = await capturePng(renderer);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    a.click();
}
