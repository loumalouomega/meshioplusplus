/**
 * Capture the documentation screenshots of the browser viewer.
 *
 * Run manually after a UI change, not in CI: the images are committed, and
 * regenerating them on every run would churn binaries for anti-aliasing noise.
 *
 *     cd viewer
 *     npm run build:web
 *     npm run preview &
 *     node tools/screenshots.mjs
 *
 * `--enable-unsafe-swiftshader` is required: headless Chromium refuses to fall
 * back to software WebGL without it, and every shot would come out blank.
 */
import { mkdir } from 'node:fs/promises';
import { fileURLToPath } from 'node:url';

import { chromium } from '@playwright/test';

const BASE = process.env.VIEWER_URL ?? 'http://localhost:4173/meshioplusplus/viewer/';
const OUT = fileURLToPath(new URL('../../../doc/public/viewer/', import.meta.url));

/** Shots to take: which sample, which array to colour by, what to call it. */
const SHOTS = [
    {
        name: 'browser-viewer',
        sample: 'Wave (point data)',
        colorBy: 'height',
        edges: false,
    },
    {
        name: 'browser-viewer-volume',
        sample: 'Block (volume)',
        colorBy: 'layer',
        edges: true,
    },
];

const browser = await chromium.launch({
    args: ['--enable-unsafe-swiftshader', '--use-gl=swiftshader'],
});
const page = await browser.newPage({
    viewport: { width: 1200, height: 750 },
    deviceScaleFactor: 2,
});

await mkdir(OUT, { recursive: true });

for (const shot of SHOTS) {
    await page.goto(BASE);
    await page.waitForFunction(() => window.__viewerState?.status === 'idle');

    await page.getByRole('button', { name: shot.sample }).click();
    await page.waitForFunction(() => window.__viewerState?.status === 'ready');

    const option = await page
        .locator('#color-by option')
        .filter({ hasText: shot.colorBy })
        .first()
        .getAttribute('value');
    if (option) await page.locator('#color-by').selectOption(option);
    if (shot.edges) await page.locator('#edges').check();

    // Let the render settle before capturing.
    await page.waitForTimeout(1200);
    await page.screenshot({ path: `${OUT}${shot.name}.png` });
    console.log(`wrote ${OUT}${shot.name}.png`);
}

await browser.close();
