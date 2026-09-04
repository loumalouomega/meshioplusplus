/**
 * Card thumbnails: a preview render, captured through the same
 * `captureNextImage` path the viewer's screenshot uses (reading the canvas
 * directly yields a black image once the frame is presented), downscaled to
 * a small JPEG data URL. Session-only — a thumbnail is a convenience, never
 * state worth persisting.
 */

import type { Renderer } from '../render/renderer';
import { capturePng } from '../render/screenshot';

export const THUMB_WIDTH = 192;
export const THUMB_HEIGHT = 128;

function loadImage(url: string): Promise<HTMLImageElement> {
    return new Promise((resolve, reject) => {
        const img = new Image();
        img.onload = () => resolve(img);
        img.onerror = () => reject(new Error('meshio++: thumbnail image failed to load'));
        img.src = url;
    });
}

/** Capture the current view as a `THUMB_WIDTH`x`THUMB_HEIGHT` JPEG data URL
 * (cover-fit, centre-cropped). */
export async function captureThumbnail(renderer: Renderer): Promise<string> {
    const img = await loadImage(await capturePng(renderer));
    const canvas = document.createElement('canvas');
    canvas.width = THUMB_WIDTH;
    canvas.height = THUMB_HEIGHT;
    const ctx = canvas.getContext('2d');
    if (!ctx) throw new Error('meshio++: no 2D canvas for the thumbnail');
    const scale = Math.max(THUMB_WIDTH / img.width, THUMB_HEIGHT / img.height);
    const w = img.width * scale;
    const h = img.height * scale;
    ctx.drawImage(img, (THUMB_WIDTH - w) / 2, (THUMB_HEIGHT - h) / 2, w, h);
    return canvas.toDataURL('image/jpeg', 0.8);
}
