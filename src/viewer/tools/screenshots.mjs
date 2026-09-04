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

// --- the dashboard (overview depth) ------------------------------------------ //
// Two manifests over the same sample files, through a mocked File System
// Access directory (the fallback input cannot hold a manifest that was not
// on disk), then "Scan all" with thumbnails so the cards carry badges and a
// preview. The same mock shape tests/viewer/dataset.spec.mjs drives.
const dashboard = await browser.newPage({
    viewport: { width: 1200, height: 750 },
    deviceScaleFactor: 2,
});
await dashboard.addInitScript(() => {
    const manifest = (name, entries) =>
        JSON.stringify({ Version: 1, Name: name, Entries: entries }, null, 2) + '\n';
    const CAMPAIGN = manifest('heat-campaign', [
        { Id: 'block', Source: { Path: 'block.vtu' }, Split: 'train', Tags: ['raw', 'coarse'] },
        { Id: 'wave', Source: { Path: 'wave.vtp' }, Split: 'valid', Tags: ['raw'] },
        { Id: 'series', Source: { Path: 'series.xdmf' }, Split: 'test', Tags: ['transient'] },
    ]);
    const SWEEP = manifest('re-sweep', [
        { Id: 'block', Source: { Path: 'block.vtu' }, Split: 'train', Group: 'cyl/laminar' },
        { Id: 'wave', Source: { Path: 'wave.vtp' }, Group: 'cyl/turbulent' },
    ]);
    const sample = async (name) =>
        (await fetch(new URL(`samples/${name}`, location.href))).arrayBuffer();
    const fileHandle = (name, getBytes, lastModified) => ({
        kind: 'file',
        name,
        getFile: async () => new File([await getBytes()], name, { lastModified }),
    });
    const handles = {
        'campaign.json': fileHandle('campaign.json', async () => CAMPAIGN, 1756800000000),
        'sweep.json': fileHandle('sweep.json', async () => SWEEP, 1756200000000),
        'block.vtu': fileHandle('block.vtu', () => sample('block.vtu'), 1755000000000),
        'wave.vtp': fileHandle('wave.vtp', () => sample('wave.vtp'), 1755000000000),
        'series.xdmf': fileHandle('series.xdmf', () => sample('series.xdmf'), 1755000000000),
    };
    window.showDirectoryPicker = async () => ({
        kind: 'directory',
        name: 'cases',
        async *entries() {
            for (const [name, handle] of Object.entries(handles)) yield [name, handle];
        },
        queryPermission: async () => 'granted',
    });
});
await dashboard.goto(DATASET_BASE);
await dashboard.waitForFunction(() => window.__datasetState?.status === 'idle');
await dashboard.locator('#ws-pick').click();
await dashboard.waitForFunction(() => window.__datasetState?.manifests.length === 2);
await dashboard.locator('#ov-scan-all').click();
await dashboard.waitForFunction(
    () =>
        window.__datasetState?.view === 'overview' &&
        window.__datasetState?.manifests.every((c) => c.health && c.thumbnail),
    null,
    { timeout: 120_000 },
);
await dashboard.waitForTimeout(1200);
await dashboard.screenshot({ path: `${OUT}dataset-dashboard.png` });
console.log(`wrote ${OUT}dataset-dashboard.png`);

await browser.close();
