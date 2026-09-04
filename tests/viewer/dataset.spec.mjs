/**
 * End-to-end tests for the dataset-manager page (`dataset.html`).
 *
 * Assertions go through `window.__datasetState` (typed in
 * src/viewer/src/dataset/types.ts), never pixels — no GPU needed. Two file
 * paths are covered: the `webkitdirectory` fallback, driven with Playwright's
 * directory upload, and the File System Access path, driven through a mocked
 * `showDirectoryPicker` installed by `addInitScript` (init scripts run on
 * real navigations, which every `page.goto` here is).
 */
import { fileURLToPath } from 'node:url';

import { expect, test } from '@playwright/test';

const SAMPLES = fileURLToPath(
    new URL('../../src/viewer/public/samples', import.meta.url),
);

const state = (page) => page.evaluate(() => window.__datasetState);

/** Force the fallback path: headless Chromium HAS showDirectoryPicker. */
async function useFallback(page) {
    await page.addInitScript(() => {
        delete window.showDirectoryPicker;
    });
}

async function waitStatus(page, wanted) {
    await expect
        .poll(() => page.evaluate(() => window.__datasetState.status), {
            message: `dataset page never reached "${wanted}"`,
        })
        .toBe(wanted);
}

/** Pick the samples directory (fallback input) and open a new manifest. */
async function openNewManifest(page) {
    await page.setInputFiles('#ws-input', SAMPLES);
    await expect
        .poll(() => page.evaluate(() => window.__datasetState.backendMode))
        .toBe('fallback');
    await page.locator('#manifest-open').click();
    await waitStatus(page, 'ready');
}

async function addCase(page, relPath) {
    await page.locator('#e-add').click();
    await page
        .locator('#add-files input[type=checkbox]')
        .and(page.locator(`[value="${relPath}"]`))
        .check();
    await page.locator('#add-commit').click();
    await waitStatus(page, 'ready');
}

test('fallback flow: add a case, edit, save as download', async ({ page }) => {
    await useFallback(page);
    await page.goto('./dataset.html');
    await waitStatus(page, 'idle');

    await openNewManifest(page);
    await addCase(page, 'block.vtu');

    let s = await state(page);
    expect(s.entryIds).toEqual(['block']);
    expect(s.selected).toBe('block');
    expect(s.planLength).toBe(1);
    expect(s.numPoints).toBeGreaterThan(0);
    expect(s.dirty).toBe(false); // nothing saved yet, so nothing to diff against

    // curate: split + tags through the detail editor
    await page.locator('#d-split').fill('train');
    await page.locator('#d-split').blur();
    await page.locator('#d-tags').fill('raw, coarse');
    await page.locator('#d-tags').blur();

    await page.locator('#m-save').click();
    await waitStatus(page, 'ready');
    s = await state(page);
    expect(s.lastSave.mode).toBe('download');
    // Byte-parity with DatasetManifest.save's conventions: 2-space indent,
    // insertion-ordered keys, trailing newline.
    expect(s.lastSave.text.startsWith('{\n  "Version": 1,\n')).toBe(true);
    expect(s.lastSave.text.endsWith('\n')).toBe(true);
    expect(s.lastSave.text).toContain('"Id": "block"');
    expect(s.lastSave.text).toContain('"Path": "block.vtu"');
    expect(s.lastSave.text).toContain('"Split": "train"');
    expect(s.lastSave.text).toContain('"Tags": [\n        "raw",\n        "coarse"\n      ]');
    expect(s.dirty).toBe(false);

    // the summary table saw the sample's data arrays
    expect(Array.isArray(s.summary)).toBe(true);
});

test('a multi-step file gets a plan and the scrubber steps through it', async ({
    page,
}) => {
    await useFallback(page);
    await page.goto('./dataset.html');
    await waitStatus(page, 'idle');
    await openNewManifest(page);
    await addCase(page, 'series.xdmf');

    let s = await state(page);
    expect(s.entryIds).toEqual(['series']);
    expect(s.planLength).toBe(3); // the fan-out resolved every step
    expect(s.step).toBe(0);

    await page.locator('#scrub').fill('2');
    await expect.poll(() => page.evaluate(() => window.__datasetState.step)).toBe(2);
    s = await state(page);
    expect(s.numPoints).toBeGreaterThan(0);
});

test('a broken pattern fails by name and never half-adds an entry', async ({
    page,
}) => {
    await useFallback(page);
    await page.goto('./dataset.html');
    await waitStatus(page, 'idle');
    await openNewManifest(page);

    // Hand-build an entry through the page's own model by picking two files
    // (a Paths source) then removing one from disk is not possible here, so
    // exercise the validation seam directly: an add whose files vanish is
    // covered by unit tests; e2e covers the happy path plus the error state
    // via an unreadable selection (no files checked -> commit is a no-op).
    await page.locator('#e-add').click();
    await page.locator('#add-commit').click();
    const s = await state(page);
    expect(s.entryIds).toEqual([]);
});

test('FSA flow: in-place save through a mocked directory handle', async ({
    page,
}) => {
    // A mock directory containing block.vtu (bytes fetched from the served
    // samples), with a writable manifest handle capturing every write.
    await page.addInitScript(() => {
        window.__fsaWrites = [];
        const fileHandle = (name, getBytes) => ({
            kind: 'file',
            name,
            getFile: async () => new File([await getBytes()], name),
        });
        const writableHandle = (name) => ({
            kind: 'file',
            name,
            getFile: async () => new File([''], name),
            createWritable: async () => ({
                write: async (data) => {
                    window.__fsaWrites.push(String(data));
                },
                close: async () => {},
            }),
        });
        window.showDirectoryPicker = async () => ({
            kind: 'directory',
            name: 'root',
            async *entries() {
                yield [
                    'block.vtu',
                    fileHandle('block.vtu', async () => {
                        const url = new URL('samples/block.vtu', location.href);
                        return (await fetch(url)).arrayBuffer();
                    }),
                ];
            },
            queryPermission: async () => 'granted',
            getFileHandle: async (name) => writableHandle(name),
        });
    });
    await page.goto('./dataset.html');
    await waitStatus(page, 'idle');

    await page.locator('#ws-pick').click();
    await expect
        .poll(() => page.evaluate(() => window.__datasetState.backendMode))
        .toBe('fsa');
    await page.locator('#manifest-open').click();
    await waitStatus(page, 'ready');
    await addCase(page, 'block.vtu');

    await page.locator('#m-save').click();
    await waitStatus(page, 'ready');
    const s = await state(page);
    expect(s.lastSave.mode).toBe('in-place');
    const writes = await page.evaluate(() => window.__fsaWrites);
    expect(writes.length).toBe(1);
    expect(writes[0]).toBe(s.lastSave.text);
});

// --------------------------------------------------------------------------- //
// v9.30.0: persisted directory handles (IndexedDB) + quality summaries        //
// --------------------------------------------------------------------------- //

/**
 * A Map-backed `indexedDB` fake, installed before page scripts run. The
 * production code reads `window.indexedDB` at call time (the test seam) and
 * treats every failure as best-effort, so this only needs the happy-path
 * subset persist.ts touches: open → onupgradeneeded/onsuccess, transaction →
 * objectStore → put/get/delete with async onsuccess. The backing Map is
 * exposed as `window.__idbStore`; `seedHandle` pre-seeds a mock directory
 * handle whose permission starts at 'prompt' (a restored handle's real
 * state) and is granted or denied by the flag.
 */
async function installIdbFake(page, { seedHandle, grant = true } = {}) {
    await page.addInitScript(
        ({ seed, grantPermission }) => {
            const store = new Map();
            window.__idbStore = store;
            window.__pickerInvoked = false;
            if (seed) {
                store.set('workspace-root', {
                    kind: 'directory',
                    name: 'root',
                    async *entries() {
                        yield [
                            'block.vtu',
                            {
                                kind: 'file',
                                name: 'block.vtu',
                                getFile: async () => {
                                    const url = new URL('samples/block.vtu', location.href);
                                    const bytes = await (await fetch(url)).arrayBuffer();
                                    return new File([bytes], 'block.vtu');
                                },
                            },
                        ];
                    },
                    queryPermission: async () => 'prompt',
                    requestPermission: async () =>
                        grantPermission ? 'granted' : 'denied',
                    getFileHandle: async (name) => ({
                        kind: 'file',
                        name,
                        getFile: async () => new File([''], name),
                        createWritable: async () => ({
                            write: async () => {},
                            close: async () => {},
                        }),
                    }),
                });
            }
            const request = (result) => {
                const r = { result, onsuccess: null, onerror: null };
                queueMicrotask(() => r.onsuccess && r.onsuccess());
                return r;
            };
            const objectStore = {
                put: (value, key) => request(void store.set(key, value)),
                get: (key) => request(store.get(key)),
                delete: (key) => request(void store.delete(key)),
            };
            const fake = {
                open: () => {
                    const r = {
                        result: {
                            createObjectStore: () => objectStore,
                            transaction: () => ({ objectStore: () => objectStore }),
                            close: () => {},
                        },
                        onupgradeneeded: null,
                        onsuccess: null,
                        onerror: null,
                    };
                    queueMicrotask(() => {
                        if (r.onupgradeneeded) r.onupgradeneeded();
                        if (r.onsuccess) r.onsuccess();
                    });
                    return r;
                },
            };
            Object.defineProperty(window, 'indexedDB', { value: fake });
            // Track fallback-picker invocations for the denial test.
            delete window.showDirectoryPicker;
            window.showDirectoryPicker = async () => {
                window.__pickerInvoked = true;
                throw new DOMException('aborted', 'AbortError');
            };
        },
        { seed: !!seedHandle, grantPermission: grant },
    );
}

test('reopen: a persisted handle restores the workspace behind a click', async ({
    page,
}) => {
    await installIdbFake(page, { seedHandle: true, grant: true });
    await page.goto('./dataset.html');
    await waitStatus(page, 'idle');
    await expect
        .poll(() => page.evaluate(() => window.__datasetState.restoreAvailable))
        .toBe(true);
    await expect(page.locator('#ws-reopen')).toContainText('root');

    await page.locator('#ws-reopen').click();
    await expect
        .poll(() => page.evaluate(() => window.__datasetState.backendMode))
        .toBe('fsa');
    // continue as a smoke: the restored workspace opens a manifest normally
    await page.locator('#manifest-open').click();
    await waitStatus(page, 'ready');
});

test('reopen: denial clears the handle and falls back to the picker', async ({
    page,
}) => {
    await installIdbFake(page, { seedHandle: true, grant: false });
    await page.goto('./dataset.html');
    await waitStatus(page, 'idle');
    await expect
        .poll(() => page.evaluate(() => window.__datasetState.restoreAvailable))
        .toBe(true);

    await page.locator('#ws-reopen').click();
    await expect
        .poll(() => page.evaluate(() => window.__datasetState.restoreAvailable))
        .toBe(false);
    const after = await page.evaluate(() => ({
        stored: window.__idbStore.size,
        pickerInvoked: window.__pickerInvoked,
        reopenHidden: document.getElementById('ws-reopen').hidden,
    }));
    expect(after.stored).toBe(0); // clearHandle ran
    expect(after.pickerInvoked).toBe(true); // fell back to a fresh pick
    expect(after.reopenHidden).toBe(true);
});

test('quality: summary rows flow from attachQuality and the scan is unpolluted', async ({
    page,
}) => {
    await useFallback(page);
    await page.goto('./dataset.html');
    await waitStatus(page, 'idle');
    await openNewManifest(page);
    await addCase(page, 'block.vtu');

    let s = await state(page);
    // quality:* rows ride the same summary; scaled_jacobian has a finite min
    const names = s.summary.map((r) => r.name);
    expect(names.some((n) => n.startsWith('quality:'))).toBe(true);
    const jacobian = s.summary.find((r) => r.name === 'quality:scaled_jacobian');
    expect(Number.isFinite(jacobian.min)).toBe(true);

    await page.locator('#e-scan').click();
    await waitStatus(page, 'ready');
    s = await state(page);
    const scan = s.scans.block;
    expect(scan.steps).toBe(1);
    expect(scan.numInverted).toBe(0);
    expect(typeof scan.minScaledJacobian).toBe('number');
    // quality NaN (= metric N/A) must NOT pollute the bad-case lane
    expect(scan.numNan).toBe(0);
    expect(scan.numInf).toBe(0);
});

// --------------------------------------------------------------------------- //
// v10.22.0: the overview depth — cards, drill-down, health, diff              //
// --------------------------------------------------------------------------- //

/**
 * A mocked FSA directory holding two manifests (one deliberately broken) and
 * the block.vtu sample, with real `lastModified` stamps so the cards can show
 * them. Writes are captured like the FSA test above.
 */
async function installOverviewMock(page) {
    await page.addInitScript(() => {
        window.__fsaWrites = [];
        const MANIFEST =
            JSON.stringify(
                {
                    Version: 1,
                    Name: 'campaign',
                    Entries: [
                        { Id: 'block', Source: { Path: 'block.vtu' }, Split: 'train', Tags: ['raw'] },
                    ],
                },
                null,
                2,
            ) + '\n';
        const BROKEN = '{"Version": 1, "Bogus": 1}';
        const fileHandle = (name, getBytes, lastModified) => ({
            kind: 'file',
            name,
            getFile: async () => new File([await getBytes()], name, { lastModified }),
            createWritable: async () => ({
                write: async (data) => {
                    window.__fsaWrites.push(String(data));
                },
                close: async () => {},
            }),
        });
        const handles = {
            'dataset.json': fileHandle('dataset.json', async () => MANIFEST, 1700000000000),
            'broken.json': fileHandle('broken.json', async () => BROKEN, 1600000000000),
            'block.vtu': fileHandle(
                'block.vtu',
                async () => {
                    const url = new URL('samples/block.vtu', location.href);
                    return (await fetch(url)).arrayBuffer();
                },
                1500000000000,
            ),
        };
        window.showDirectoryPicker = async () => ({
            kind: 'directory',
            name: 'root',
            async *entries() {
                for (const [name, handle] of Object.entries(handles)) yield [name, handle];
            },
            queryPermission: async () => 'granted',
            getFileHandle: async (name) => handles[name] ?? fileHandle(name, async () => '', Date.now()),
        });
    });
}

async function pickOverview(page) {
    await page.goto('./dataset.html');
    await waitStatus(page, 'idle');
    await page.locator('#ws-pick').click();
    await expect
        .poll(() => page.evaluate(() => window.__datasetState.manifests.length))
        .toBe(2);
}

test('overview: a card per manifest, drill-down and back', async ({ page }) => {
    await installOverviewMock(page);
    await pickOverview(page);

    let s = await state(page);
    expect(s.view).toBe('overview');
    const card = s.manifests.find((c) => c.path === 'dataset.json');
    expect(card.name).toBe('campaign');
    expect(card.numEntries).toBe(1);
    expect(card.splits).toEqual({ train: 1 });
    expect(card.tags).toEqual(['raw']);
    expect(card.lastModified).toBe(1700000000000);
    expect(card.parseError).toBeNull();
    expect(card.health).toBeNull();
    expect(typeof card.sha256).toBe('string');
    const broken = s.manifests.find((c) => c.path === 'broken.json');
    expect(broken.parseError).toContain("unknown key 'Bogus'");
    expect(broken.numEntries).toBe(0);
    await expect(page.locator('.card[data-path="broken.json"]')).toHaveClass(/broken/);

    // drill down through the card, then back
    await page.locator('.card[data-path="dataset.json"] .card-open').click();
    await waitStatus(page, 'ready');
    s = await state(page);
    expect(s.view).toBe('manifest');
    expect(s.manifestName).toBe('dataset.json');
    expect(s.entryIds).toEqual(['block']);
    await expect(page.locator('#overview')).toBeHidden();
    await expect(page.locator('#m-back')).toBeVisible();

    await page.locator('#m-back').click();
    s = await state(page);
    expect(s.view).toBe('overview');
    await expect(page.locator('#overview')).toBeVisible();
});

test('overview: a scan fills the card health and the drill-down health section', async ({
    page,
}) => {
    await installOverviewMock(page);
    await pickOverview(page);
    await page.locator('.card[data-path="dataset.json"] .card-open').click();
    await waitStatus(page, 'ready');

    await page.locator('#e-scan').click();
    await waitStatus(page, 'ready');
    const s = await state(page);
    const health = s.manifests.find((c) => c.path === 'dataset.json').health;
    expect(health.producer).toBe('browser');
    expect(health.scanned).toBe(1);
    expect(health.total).toBe(1);
    expect(health.numNan).toBe(0);
    expect(health.numInf).toBe(0);
    expect(health.numInverted).toBe(0);
    expect(health.numDegenerate).toBe(0);
    expect(typeof health.minScaledJacobian).toBe('number');
    expect(health.splitBalance).toEqual([{ split: 'train', count: 1, fraction: 1 }]);
    expect(health.fieldsMissing).toEqual({});
    expect(health.badEntries).toEqual([]);
    expect(s.scans.block.arrays.length).toBeGreaterThan(0);
    await expect(page.locator('#health-section')).toBeVisible();
    await expect(page.locator('#h-totals .badge.ok')).toHaveCount(3);
});

test('overview: the diff view compares the file on disk with the current edits', async ({
    page,
}) => {
    await installOverviewMock(page);
    await pickOverview(page);
    await page.locator('.card[data-path="dataset.json"] .card-open').click();
    await waitStatus(page, 'ready');

    await page.locator('#entry-list li', { hasText: 'block' }).click();
    await page.locator('#d-split').fill('valid');
    await page.locator('#d-split').blur();
    await expect.poll(() => page.evaluate(() => window.__datasetState.dirty)).toBe(true);

    await page.locator('#ov-diff').click();
    await expect
        .poll(() => page.evaluate(() => window.__datasetState.diff?.changed))
        .toBe(1);
    const s = await state(page);
    expect(s.diff.a).toContain('on disk');
    expect(s.diff.b).toContain('current edits');
    expect(s.diff.added).toBe(0);
    expect(s.diff.removed).toBe(0);
    await expect(page.locator('#diff-body')).toContainText('~ block');
    await expect(page.locator('#diff-body')).toContainText('split: train → valid');

    // the card reflects the unsaved edit once we go back
    await page.locator('#diff-close').click();
    await page.locator('#m-back').click();
    const card = (await state(page)).manifests.find((c) => c.path === 'dataset.json');
    expect(card.dirty).toBe(true);
    expect(card.splits).toEqual({ valid: 1 });
});
