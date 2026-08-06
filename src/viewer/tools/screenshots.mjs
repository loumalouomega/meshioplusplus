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

// --- the dataset-manager page ---------------------------------------------- //
// Populated with a few cases from the same sample files, via the
// webkitdirectory fallback (the same path tests/viewer/dataset.spec.mjs
// drives) -- no File System Access dialog to simulate for a screenshot.
const DATASET_BASE = BASE.replace(/\/$/, '') + '/dataset.html';
const SAMPLES_DIR = fileURLToPath(new URL('../public/samples', import.meta.url));

async function addCase(page, relPath, { split, tags } = {}) {
    await page.locator('#e-add').click();
    await page
        .locator('#add-files input[type=checkbox]')
        .and(page.locator(`[value="${relPath}"]`))
        .check();
    await page.locator('#add-commit').click();
    await page.waitForFunction(() => window.__datasetState?.status === 'ready');
    if (split) {
        await page.locator('#d-split').fill(split);
        await page.locator('#d-split').blur();
    }
    if (tags) {
        await page.locator('#d-tags').fill(tags);
        await page.locator('#d-tags').blur();
    }
}

await page.goto(DATASET_BASE);
await page.waitForFunction(() => window.__datasetState?.status === 'idle');
await page.setInputFiles('#ws-input', SAMPLES_DIR);
await page.waitForFunction(() => window.__datasetState?.backendMode === 'fallback');
await page.locator('#manifest-open').click();
await page.waitForFunction(() => window.__datasetState?.status === 'ready');

await addCase(page, 'block.vtu', { split: 'train', tags: 'raw, coarse' });
await addCase(page, 'wave.vtp', { split: 'valid' });
await addCase(page, 'series.xdmf', { split: 'test', tags: 'transient' });

// Select the multi-step case last, so the screenshot shows the scrubber.
await page.locator('#entry-list li', { hasText: 'series' }).click();
await page.locator('#d-preview').click();
await page.waitForFunction(() => window.__datasetState?.status === 'ready');

const tOption = await page
    .locator('#color-by option')
    .filter({ hasText: 'T' })
    .first()
    .getAttribute('value');
if (tOption) await page.locator('#color-by').selectOption(tOption);

// The detail editor scrolls the panel down when it opens; scroll back to the
// top so the workspace/manifest header is in frame too.
await page.locator('#panel').evaluate((el) => el.scrollTo(0, 0));

await page.waitForTimeout(1200);
await page.screenshot({ path: `${OUT}dataset-manager.png` });
console.log(`wrote ${OUT}dataset-manager.png`);

await browser.close();
