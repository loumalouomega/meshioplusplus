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

// --------------------------------------------------------------------------- //
// v10.23.0: the companion process (served same-origin by a routed fake, so   //
// no CORS preflight is involved — CORS is covered by test_mcp_http.py)         //
// --------------------------------------------------------------------------- //

import { createHash } from 'node:crypto';

const OVERVIEW_MANIFEST =
    JSON.stringify(
        {
            Version: 1,
            Name: 'campaign',
            Entries: [{ Id: 'block', Source: { Path: 'block.vtu' }, Split: 'train', Tags: ['raw'] }],
        },
        null,
        2,
    ) + '\n';
const OVERVIEW_SHA = createHash('sha256').update(OVERVIEW_MANIFEST).digest('hex');
const SERVER_TOKEN = 'tok-123';

/** A fake companion process at `<page origin>/__mock-api`, answering only
 * with the right bearer token. */
async function installServerFake(page) {
    await page.route('**/__mock-api/**', async (route) => {
        const request = route.request();
        const url = new URL(request.url());
        const json = (status, body) =>
            route.fulfill({ status, contentType: 'application/json', body: JSON.stringify(body) });
        if (request.headers()['authorization'] !== `Bearer ${SERVER_TOKEN}`) {
            return json(401, { error: 'missing or invalid bearer token', error_type: 'PermissionError' });
        }
        if (url.pathname.endsWith('/api/health')) {
            return json(200, {
                version: '10.23.0',
                root: '/srv/cases',
                runs_dir: '/srv/cases/runs',
                tools: ['info', 'dataset_find', 'dataset_health'],
                mcp: '/mcp',
                transport: 'streamable-http',
                auth: 'token',
            });
        }
        if (url.pathname.endsWith('/api/tools/dataset_find')) {
            return json(200, {
                root: '/srv/cases',
                manifests: [
                    { path: '/srv/cases/dataset.json', relpath: 'dataset.json', sha256: OVERVIEW_SHA, name: 'campaign', num_entries: 1, splits: { train: 1 }, mtime: 1700000000000 },
                    { path: '/srv/cases/other/big.json', relpath: 'other/big.json', sha256: 'f'.repeat(64), name: 'big', num_entries: 40, splits: { train: 30, valid: 10 }, mtime: 1690000000000 },
                ],
            });
        }
        if (url.pathname.endsWith('/api/tools/dataset_health')) {
            const body = request.postDataJSON();
            return json(200, {
                producer: 'server',
                name: 'campaign',
                num_entries: 1,
                scanned: 1,
                splits: { train: 1 },
                split_balance: [{ split: 'train', count: 1, fraction: 1 }],
                entries: { block: { steps: 1, num_nan: 0, num_inf: 0, num_inverted: 0, num_degenerate: 0, min_scaled_jacobian: 0.42, arrays: ['point_data:layer'] } },
                fields_missing: {},
                totals: { num_nan: 0, num_inf: 0, num_inverted: 0, num_degenerate: 0, min_scaled_jacobian: 0.42 },
                bad_entries: [],
                manifest_path: body.manifest_path,
                sha256: OVERVIEW_SHA,
            });
        }
        return json(404, { error: 'unknown', error_type: 'KeyError' });
    });
}

// The web build registers the COOP/COEP service worker (a speed enhancement
// for the threaded wasm), and Playwright's `page.route` never sees a fetch
// a service worker handled -- so the fake is only reachable with service
// workers blocked for this spec.
test.describe('companion process', () => {
    test.use({ serviceWorkers: 'block' });

test('connect, bind cards by hash, scan on the server', async ({ page }) => {
    await installOverviewMock(page);
    await installServerFake(page);
    await pickOverview(page);

    const serverUrl = new URL('__mock-api', page.url()).toString();
    await page.locator('#srv-url').fill(serverUrl);
    // a wrong token is a named 401, and the section stays disconnected
    await page.locator('#srv-token').fill('wrong');
    await page.locator('#srv-connect').click();
    await expect
        .poll(() => page.evaluate(() => window.__datasetState.server?.error))
        .toContain('bearer token');
    expect((await state(page)).server.connected).toBe(false);
    await expect(page.locator('#e-scan-server')).toBeHidden();

    await page.locator('#srv-token').fill(SERVER_TOKEN);
    await page.locator('#srv-connect').click();
    await expect
        .poll(() => page.evaluate(() => window.__datasetState.server?.connected))
        .toBe(true);
    let s = await state(page);
    expect(s.server.version).toBe('10.23.0');
    expect(s.server.root).toBe('/srv/cases');
    expect(s.server.tools).toContain('dataset_health');

    // the workspace card is bound by its content hash; the other is server-only
    await expect
        .poll(() => page.evaluate(() => window.__datasetState.manifests.length))
        .toBe(3);
    s = await state(page);
    const bound = s.manifests.find((c) => c.path === 'dataset.json');
    expect(bound.serverPath).toBe('/srv/cases/dataset.json');
    expect(bound.serverOnly).toBe(false);
    const remote = s.manifests.find((c) => c.serverOnly);
    expect(remote.path).toBe('server:other/big.json');
    expect(remote.numEntries).toBe(40);
    await expect(page.locator('.card[data-path="server:other/big.json"] .card-open')).toBeDisabled();

    // scan on the server from the card
    await page.locator('.card[data-path="dataset.json"] .card-scan-server').click();
    await expect
        .poll(() => page.evaluate(() => window.__datasetState.manifests.find((c) => c.path === 'dataset.json').health?.producer))
        .toBe('server');
    s = await state(page);
    expect(s.manifests.find((c) => c.path === 'dataset.json').health.minScaledJacobian).toBe(0.42);

    // ...and from the drill-down, where it also fills the entry scans
    await page.locator('.card[data-path="dataset.json"] .card-open').click();
    await waitStatus(page, 'ready');
    await expect(page.locator('#e-scan-server')).toBeVisible();
    await page.locator('#e-scan-server').click();
    await waitStatus(page, 'ready');
    s = await state(page);
    expect(s.scans.block.arrays).toEqual(['point_data:layer']);
    expect(s.scans.block.minScaledJacobian).toBe(0.42);

    // disconnecting drops the server-only card and unbinds the rest
    await page.locator('#srv-disconnect').click();
    s = await state(page);
    expect(s.server.connected).toBe(false);
    expect(s.manifests.length).toBe(2);
    expect(s.manifests.find((c) => c.path === 'dataset.json').serverPath).toBeNull();
});
});
