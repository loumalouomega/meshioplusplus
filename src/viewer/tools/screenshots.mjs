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
import { createHash } from 'node:crypto';
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

// --- a training run (doc/dashboard.md) --------------------------------------- //
// Against a stub companion process, routed same-origin: the shot documents the
// run panel, and standing up a real server (and a real GPU run) for a picture
// would be a much bigger dependency than the picture is worth. Service workers
// are blocked because the COOP/COEP worker would answer the fetch first.
const training = await browser.newContext({
    viewport: { width: 1200, height: 900 },
    deviceScaleFactor: 2,
    serviceWorkers: 'block',
});
const trainingPage = await training.newPage();
// A workspace holding the manifest the stub server also reports, bound by its
// content hash exactly as the real page binds them -- so the launch form is
// populated from a real open manifest rather than sitting empty.
const TRAIN_MANIFEST =
    JSON.stringify(
        {
            Version: 1,
            Name: 'heat-campaign',
            Entries: [
                { Id: 'block', Source: { Path: 'block.vtu' }, Split: 'train', Tags: ['raw'] },
            ],
        },
        null,
        2,
    ) + '\n';
const TRAIN_SHA = createHash('sha256').update(TRAIN_MANIFEST).digest('hex');
await trainingPage.addInitScript((manifestText) => {
    const fileHandle = (name, getBytes, lastModified) => ({
        kind: 'file',
        name,
        getFile: async () => new File([await getBytes()], name, { lastModified }),
    });
    const handles = {
        'dataset_manifest.json': fileHandle('dataset_manifest.json', async () => manifestText, 1756800000000),
        'block.vtu': fileHandle(
            'block.vtu',
            async () => (await fetch(new URL('samples/block.vtu', location.href))).arrayBuffer(),
            1755000000000,
        ),
    };
    window.showDirectoryPicker = async () => ({
        kind: 'directory',
        name: 'cases',
        async *entries() {
            for (const [name, handle] of Object.entries(handles)) yield [name, handle];
        },
        queryPermission: async () => 'granted',
    });
}, TRAIN_MANIFEST);
const EPOCHS = 100;
const rows = Array.from({ length: EPOCHS }, (_, epoch) => ({
    epoch,
    train_loss: 0.79 * Math.exp(-epoch / 12) + 9e-4 + 3e-4 * Math.sin(epoch / 3),
    valid_loss: 0.25 * Math.exp(-epoch / 11) + 6e-4 + 4e-4 * Math.cos(epoch / 4),
    lr: 0.001,
    elapsed: (epoch + 1) * 0.5,
    epoch_seconds: 0.5,
    timestamp: 1700000000 + epoch,
}));
const ckpt = (name, kind, epoch, loss, best) => ({
    path: `/srv/cases/runs/20260904-101500-3f2a/checkpoints/${name}`,
    name,
    kind,
    epoch,
    valid_loss: loss,
    size: 1_310_720,
    is_best: best,
});
await trainingPage.route('**/__mock-api/**', async (route) => {
    const url = new URL(route.request().url());
    const json = (body) =>
        route.fulfill({ status: 200, contentType: 'application/json', body: JSON.stringify(body) });
    const tail = url.pathname.split('/').pop();
    if (url.pathname.endsWith('/api/health')) {
        return json({
            version: '10.24.0',
            root: '/srv/cases',
            runs_dir: '/srv/cases/runs',
            tools: ['train_start'],
            mcp: '/mcp',
            transport: 'streamable-http',
            auth: 'token',
        });
    }
    if (tail === 'dataset_find') {
        return json({
            root: '/srv/cases',
            manifests: [
                {
                    path: '/srv/cases/dataset_manifest.json',
                    relpath: 'dataset_manifest.json',
                    sha256: TRAIN_SHA,
                    name: 'heat-campaign',
                    num_entries: 200,
                    splits: { train: 160, valid: 20, test: 20 },
                    mtime: 1756800000000,
                },
            ],
        });
    }
    if (tail === 'train_defaults') {
        return json({
            manifest_path: '/srv/cases/dataset_manifest.json',
            num_entries: 200,
            splits: { train: 160, valid: 20, test: 20 },
            available_fields: { point: ['q_scaled', 'T'], cell: [] },
            runs_dir: '/srv/cases/runs',
            frameworks: { torch_geometric: true, physicsnemo: true },
            spec: {},
        });
    }
    if (tail === 'train_status') {
        return json({
            job_id: '20260904-101500-3f2a',
            run_dir: '/srv/cases/runs/20260904-101500-3f2a',
            status: 'finished',
            pid: 4242,
            started: 1700000000,
            finished: 1700000050,
            exit_code: 0,
            manifest: '/srv/cases/dataset_manifest.json',
            best_checkpoint: '/srv/cases/runs/20260904-101500-3f2a/checkpoints/best.mdlus',
            epoch: EPOCHS,
            epochs: EPOCHS,
            best_epoch: 97,
            best_valid_loss: 1.648e-4,
            eta_seconds: 0,
            device: 'cuda',
            completed: true,
            num_metrics: EPOCHS,
            last: rows[EPOCHS - 1],
        });
    }
    if (tail === 'train_metrics') {
        const body = route.request().postDataJSON() ?? {};
        const scale = body.job_id === '20260903-173000-9c11' ? 1.7 : 1;
        return json({ rows: rows.map((r) => ({ ...r, valid_loss: r.valid_loss * scale })) });
    }
    if (tail === 'train_log') {
        const text = rows
            .filter((r) => r.epoch % 10 === 0 || r.epoch === EPOCHS - 1)
            .map(
                (r) =>
                    `epoch ${String(r.epoch).padStart(4)}  train ${r.train_loss.toExponential(3)}` +
                    `  valid ${r.valid_loss.toExponential(3)}  (0.5 s)\n`,
            )
            .join('');
        return json({ text, offset: 0, next_offset: text.length, size: text.length, done: true });
    }
    if (tail === 'train_checkpoints') {
        return json({
            best_checkpoint: '/srv/cases/runs/20260904-101500-3f2a/checkpoints/best.mdlus',
            checkpoints: [
                ckpt('MeshGraphNet.0.24.mdlus', 'periodic', 24, 1.9e-3, false),
                ckpt('MeshGraphNet.0.49.mdlus', 'periodic', 49, 7.4e-4, false),
                ckpt('best.mdlus', 'best', 97, 1.648e-4, true),
                ckpt('final.mdlus', 'final', 99, 1.219e-3, false),
            ],
        });
    }
    if (tail === 'train_list') {
        const base = {
            status: 'finished',
            manifest: '/srv/cases/dataset_manifest.json',
            fields: ['q_scaled'],
            target_fields: ['T'],
            train_split: 'train',
            valid_split: 'valid',
            epochs: 100,
            batch_size: 8,
            learning_rate: 0.001,
            seed: 0,
            processor_size: 8,
            tags: [],
        };
        return json({
            runs_dir: '/srv/cases/runs',
            jobs: [
                { ...base, job_id: '20260904-101500-3f2a', hidden_dim: 64, best_valid_loss: 1.648e-4, duration_seconds: 50, tags: ['baseline'] },
                { ...base, job_id: '20260903-173000-9c11', hidden_dim: 128, batch_size: 16, best_valid_loss: 2.71e-4, duration_seconds: 88, tags: ['wide'] },
                { ...base, job_id: '20260902-090012-4ab7', status: 'failed', hidden_dim: 64, best_valid_loss: null, duration_seconds: 4 },
            ],
        });
    }
    return json({ error: `unknown ${tail}`, error_type: 'KeyError' });
});
await trainingPage.goto(DATASET_BASE);
await trainingPage.waitForFunction(() => window.__datasetState?.status === 'idle');
await trainingPage.locator('#ws-pick').click();
await trainingPage.waitForFunction(() => window.__datasetState?.manifests.length === 1);
await trainingPage.locator('#srv-url').fill(new URL('__mock-api', trainingPage.url()).toString());
await trainingPage.locator('#srv-token').fill('token');
await trainingPage.locator('#srv-connect').click();
await trainingPage.waitForFunction(() => window.__datasetState?.server?.connected === true);
await trainingPage.locator('.card[data-path="dataset_manifest.json"] .card-open').click();
await trainingPage.waitForFunction(() => window.__datasetState?.view === 'manifest');
await trainingPage.locator('#t-fields').selectOption(['q_scaled']);
await trainingPage.locator('#t-targets').selectOption(['T']);
await trainingPage.locator('#t-runs').click();
await trainingPage.waitForSelector('#runs-wrap:not([hidden])');
await trainingPage.locator('.run-open', { hasText: '20260904-101500-3f2a' }).click();
await trainingPage.waitForFunction(
    () => window.__datasetState?.activeJob?.metrics.length === 100,
);
await trainingPage.waitForTimeout(800);
await trainingPage.screenshot({ path: `${OUT}training-run.png` });
console.log(`wrote ${OUT}training-run.png`);

// --- the run history and comparison ------------------------------------------ //
await trainingPage.locator('#run-close').click();
await trainingPage.locator('#t-runs').click();
await trainingPage.waitForSelector('#runs-wrap:not([hidden])');
for (const id of ['20260904-101500-3f2a', '20260903-173000-9c11']) {
    await trainingPage
        .locator('#runs-table tr', { hasText: id })
        .locator('input')
        .check();
}
await trainingPage.waitForFunction(() => window.__datasetState?.compare.length === 2);
await trainingPage.waitForTimeout(800);
await trainingPage.screenshot({ path: `${OUT}run-history.png` });
console.log(`wrote ${OUT}run-history.png`);

await browser.close();
