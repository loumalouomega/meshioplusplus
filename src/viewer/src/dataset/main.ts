/**
 * The dataset-manager page (`dataset.html`) — build and curate
 * `DatasetManifest` JSONs against a local case directory, previewing entries
 * through the same worker/renderer the viewer uses.
 *
 * Two depths of one page: the **overview** (a card per manifest found at the
 * workspace root — counts, split balance, health, thumbnail) and the
 * **manifest** view (one manifest's curation: entries, detail editor,
 * preview, scan). Structure mirrors `src/main.ts`: static markup in the HTML,
 * `$()` lookups that throw on a missing id, one state object
 * (`window.__datasetState`) that the Playwright spec asserts on. All document
 * logic lives in the pure modules (`manifest.ts`, `glob.ts`, `suggest.ts`,
 * `health.ts`, `overview.ts`, `diff.ts`); all file access in `fs.ts`;
 * everything wasm-side behind the shared worker.
 */

import { Legend } from '../render/legend';
import { Renderer, SOLID_COLOR, colorKey } from '../render/renderer';
import { $, setOptions, show } from '../ui/dom';
import type { ArrayEntry } from '../types';
import * as client from '../worker/client';
import type { ArraySummary, PlanEntry } from '../worker/protocol';
import { ApiClient, isToolError, type ServerManifest } from './api';
import { diffLines, diffManifests, renderDiff, summarizeDiff } from './diff';
import {
    hasFsAccess,
    pickWorkspaceFsa,
    workspaceFromFileList,
    workspaceFromHandle,
    type Workspace,
} from './fs';
import { globMatch } from './glob';
import { sha256Hex } from './hash';
import {
    emptyScan,
    entryScanFromSummaries,
    fromServerHealth,
    healthBadges,
    healthFromScans,
    isQualityRow,
    qualityRowIsBad,
    type ServerHealthReport,
} from './health';
import {
    addEntry,
    deriveId,
    emptyManifest,
    getEntry,
    parseManifest,
    removeEntry,
    stringifyManifest,
    type Manifest,
    type ManifestEntry,
    type SourceSpec,
} from './manifest';
import {
    buildCard,
    cardFromManifest,
    cardFromServerManifest,
    renderCards,
    sortCards,
    splitBar,
    splitPaletteOf,
    type CardSort,
} from './overview';
import { notificationState } from './notify';
import { clearHandle, getValue, loadHandle, putValue, saveHandle } from './persist';
import { suggestSource } from './suggest';
import { TrainPanel } from './train';
import { captureThumbnail } from './thumbs';
import type {
    DatasetState,
    DatasetView,
    EntryScan,
    ManifestCard,
    ManifestHealth,
    ServerState,
} from './types';

// --- state ----------------------------------------------------------------- //

window.__datasetState = {
    status: 'booting',
    backendMode: null,
    restoreAvailable: false,
    view: 'overview',
    manifests: [],
    server: null,
    jobs: [],
    activeJob: null,
    notifications: notificationState(),
    manifestName: null,
    entryIds: [],
    selected: null,
    planLength: 0,
    step: 0,
    numPoints: 0,
    numCells: 0,
    summary: null,
    scans: {},
    diff: null,
    dirty: false,
    lastSave: null,
    error: null,
} satisfies DatasetState;

let workspace: Workspace | null = null;
let manifest: Manifest = emptyManifest();
let manifestPath: string | null = null;
/** The serialization on disk (or last saved) — dirty = current !== saved. */
let savedText: string | null = null;
let selected: string | null = null;
/** The staged entry (id) and its plan, for the scrubber. */
let staged: { id: string; plan: PlanEntry[] } | null = null;
const scans = new Map<string, EntryScan>();
let readableExtensions: Set<string> = new Set();
/** The overview's cards — one per manifest at the workspace root (plus an
 * unsaved new manifest while it is being edited). */
let cards: ManifestCard[] = [];
/** The fixed split -> colour assignment over every split in the workspace. */
let splitPalette = new Map<string, string>();
/** Set while "Scan all" drives the manifest view, so it does not flip the
 * page between depths on every manifest it opens. */
let bulkScan = false;
/** The companion-process client while connected (`api.ts`). */
let serverApi: ApiClient | null = null;
/** What the server found (`dataset_find`), merged into the cards by hash. */
let serverManifests: ServerManifest[] = [];

function setState(patch: Partial<DatasetState>): void {
    Object.assign(window.__datasetState, patch);
}

function setStatus(status: DatasetState['status'], message = ''): void {
    setState({ status, error: status === 'error' ? message : null });
    $('status').textContent =
        status === 'error' ? `Error: ${message}` : message || status;
}

function fail(e: unknown): void {
    setStatus('error', e instanceof Error ? e.message : String(e));
}

function currentText(): string {
    return stringifyManifest(manifest);
}

function currentCard(): ManifestCard | undefined {
    return manifestPath ? cards.find((c) => c.path === manifestPath) : undefined;
}

function refreshDirty(): void {
    const dirty = savedText !== null && currentText() !== savedText;
    setState({ dirty });
    $('m-dirty').hidden = !dirty;
    const card = currentCard();
    if (card) card.dirty = dirty;
}

// --- the two depths -------------------------------------------------------- //

function setView(view: DatasetView): void {
    setState({ view });
    // The CSS hides the curation sections (editor, entries, health, detail)
    // while the overview shows, without touching their `hidden` attributes —
    // those still mean "no manifest open"/"nothing selected".
    document.body.dataset.view = view;
    show($('overview'), view === 'overview');
    show($('m-back'), view === 'manifest');
    if (view === 'overview') {
        $('diff-wrap').hidden = true;
        // The run panel is a stage overlay like the diff: hidden here, but
        // its poller keeps running, so a run still finishes and announces
        // itself while the overview is up.
        $('run-wrap').hidden = true;
    }
}

// --- rendering ------------------------------------------------------------- //

const renderer = new Renderer($('render'));
let arrays: ArrayEntry[] = [];
const legend = new Legend($('legend'), {
    onRange: (min, max) => {
        renderer.applyRange(min, max);
        const range = renderer.range;
        if (range) legend.update({ ...range, min, max });
    },
    onRescaleToData: () => {
        const range = renderer.dataRange();
        if (!range) return;
        renderer.applyRange(range.min, range.max);
        legend.update(range);
    },
});

function show3d(vtp: ArrayBuffer): void {
    const info = renderer.load(vtp);
    arrays = info.arrays;
    setState({ numPoints: info.numPoints, numCells: info.numCells });
    setOptions(
        $<HTMLSelectElement>('color-by'),
        [
            { value: SOLID_COLOR, label: 'Solid colour' },
            ...arrays.map((a) => ({ value: colorKey(a), label: a.label })),
        ],
        SOLID_COLOR,
    );
    legend.hide();
}

function applyColorBy(key: string): void {
    const range = renderer.setColorBy(key);
    const entry = arrays.find((a) => colorKey(a) === key);
    if (!range || !entry) {
        legend.hide();
        return;
    }
    legend.show(entry.label, renderer.colormap, range);
}

// --- overview -------------------------------------------------------------- //

function renderOverview(): void {
    splitPalette = splitPaletteOf(cards.flatMap((c) => Object.keys(c.splits)));
    const key = $<HTMLSelectElement>('ov-sort').value as CardSort;
    renderCards($('ov-cards'), sortCards(cards, key), splitPalette, {
        onOpen: (path) => void openManifestPath(path).catch(fail),
        onScan: (path) => void scanManifestPath(path).catch(fail),
        onScanServer: serverApi
            ? (path) => void scanOnServer(path).catch(fail)
            : undefined,
    });
    show($('ov-empty'), cards.length === 0);
    $('ov-counts').textContent = `(${cards.length})`;
    setState({ manifests: cards });
}

/** Rebuild the cards from the workspace's root-level manifests, keeping the
 * health/thumbnail of any card whose file is unchanged. */
async function buildOverview(): Promise<void> {
    if (!workspace) return;
    const byPath = new Map(workspace.files.map((f) => [f.relPath, f]));
    const built: ManifestCard[] = [];
    for (const path of workspace.manifestCandidates()) {
        const file = byPath.get(path);
        if (!file) continue;
        const [text, mtime] = await Promise.all([
            workspace.readText(path),
            file.lastModified(),
        ]);
        const card = buildCard(path, text, mtime);
        card.sha256 = await sha256Hex(text);
        const previous = cards.find((c) => c.path === path);
        if (previous && previous.sha256 === card.sha256) {
            card.health = previous.health;
            card.thumbnail = previous.thumbnail;
        }
        built.push(card);
    }
    cards = built;
    mergeServerCards();
    renderOverview();
}

/** Bring the open manifest's card in line with the in-memory document
 * (counts/splits reflect unsaved edits, flagged dirty). A brand-new manifest
 * gets a card of its own until it is saved. */
function refreshCurrentCard(): void {
    if (!manifestPath) return;
    let card = currentCard();
    const fresh = cardFromManifest(manifestPath, manifest, card?.lastModified ?? null);
    if (card) {
        card.name = fresh.name;
        card.description = fresh.description;
        card.numEntries = fresh.numEntries;
        card.splits = fresh.splits;
        card.tags = fresh.tags;
        card.groups = fresh.groups;
    } else {
        card = fresh;
        cards.push(card);
    }
    card.dirty = window.__datasetState.dirty;
    if (scans.size) card.health = healthFromScans(manifest, scans);
    renderOverview();
}

// --- entry list / detail --------------------------------------------------- //

function chip(text: string, cls: string): HTMLSpanElement {
    const span = document.createElement('span');
    span.className = `chip ${cls}`;
    span.textContent = text;
    return span;
}

function matchesFilter(entry: ManifestEntry, needle: string): boolean {
    if (!needle) return true;
    const hay = [entry.id, entry.split ?? '', entry.group ?? '', ...entry.tags]
        .join(' ')
        .toLowerCase();
    return hay.includes(needle.toLowerCase());
}

function renderEntryList(): void {
    const needle = $<HTMLInputElement>('e-filter').value.trim();
    const items = manifest.entries
        .filter((entry) => matchesFilter(entry, needle))
        .map((entry) => {
            const li = document.createElement('li');
            if (entry.id === selected) li.classList.add('selected');
            const id = document.createElement('span');
            id.className = 'eid';
            id.textContent = entry.id;
            li.append(id);
            if (entry.split) li.append(chip(entry.split, 'split'));
            for (const tag of entry.tags) li.append(chip(tag, 'tag'));
            const scan = scans.get(entry.id);
            if (scan && scan.steps > 1) li.append(chip(`${scan.steps} steps`, 'tag'));
            if (scan && scan.numNan + scan.numInf > 0) {
                li.append(chip(`NaN/Inf ${scan.numNan + scan.numInf}`, 'warn'));
            }
            if (scan && scan.numInverted > 0) {
                li.append(chip(`inverted ${scan.numInverted}`, 'warn'));
            }
            if (scan && scan.numDegenerate > 0) {
                li.append(chip(`degenerate ${scan.numDegenerate}`, 'warn'));
            }
            if (scan && scan.steps === 0) li.append(chip('unreadable', 'warn'));
            li.addEventListener('click', () => selectEntry(entry.id));
            return li;
        });
    $('entry-list').replaceChildren(...items);
    $('e-counts').textContent = `(${manifest.entries.length})`;
    setState({ entryIds: manifest.entries.map((e) => e.id) });
    refreshDirty();
}

function describeSource(source: SourceSpec): string {
    if (source.Pattern) return `pattern ${source.Pattern}`;
    if (source.Path) return source.Path;
    return `${source.Paths?.length ?? 0} files: ${source.Paths?.join(', ') ?? ''}`;
}

function renderDetail(): void {
    const section = $('detail-section');
    if (!selected) {
        section.hidden = true;
        return;
    }
    const entry = getEntry(manifest, selected);
    $('d-id').textContent = entry.id;
    $('d-source').textContent = describeSource(entry.source);
    $<HTMLInputElement>('d-split').value = entry.split ?? '';
    $<HTMLInputElement>('d-tags').value = entry.tags.join(', ');
    $<HTMLInputElement>('d-group').value = entry.group ?? '';
    $<HTMLTextAreaElement>('d-notes').value = entry.notes ?? '';
    $<HTMLTextAreaElement>('d-meta').value = Object.keys(entry.metadata).length
        ? JSON.stringify(entry.metadata, null, 1)
        : '';
    $('d-meta-error').hidden = true;
    const options = [...new Set(manifest.entries.map((e) => e.split).filter(Boolean))];
    $('split-options').replaceChildren(
        ...options.map((value) => {
            const option = document.createElement('option');
            option.value = value ?? '';
            return option;
        }),
    );
    show(section);
}

function selectEntry(id: string | null): void {
    selected = id;
    setState({ selected: id });
    renderEntryList();
    renderDetail();
}

function commitDetail(): void {
    if (!selected) return;
    const entry = getEntry(manifest, selected);
    const split = $<HTMLInputElement>('d-split').value.trim();
    const tags = $<HTMLInputElement>('d-tags')
        .value.split(',')
        .map((t) => t.trim())
        .filter(Boolean);
    const group = $<HTMLInputElement>('d-group').value.trim();
    const notes = $<HTMLTextAreaElement>('d-notes').value;
    const metaText = $<HTMLTextAreaElement>('d-meta').value.trim();
    let metadata: Record<string, unknown> = {};
    if (metaText) {
        try {
            const parsed: unknown = JSON.parse(metaText);
            if (typeof parsed !== 'object' || parsed === null || Array.isArray(parsed)) {
                throw new Error('metadata must be a JSON object');
            }
            metadata = parsed as Record<string, unknown>;
        } catch (e) {
            const box = $('d-meta-error');
            box.textContent = e instanceof Error ? e.message : String(e);
            box.hidden = false;
            return;
        }
    }
    $('d-meta-error').hidden = true;
    const index = manifest.entries.indexOf(entry);
    manifest.entries[index] = {
        ...entry,
        split: split || null,
        tags,
        group: group || null,
        notes: notes || null,
        metadata,
    };
    renderEntryList();
    updateHealth();
}

// --- health ---------------------------------------------------------------- //

function renderHealthSection(health: ManifestHealth | null): void {
    const section = $('health-section');
    if (!health) {
        section.hidden = true;
        return;
    }
    $('h-splits').replaceChildren(splitBar(health.splitBalance, splitPalette));
    $('h-totals').replaceChildren(
        ...healthBadges(health).map((b) => {
            const span = document.createElement('span');
            span.className = `badge ${b.level}`;
            span.textContent = b.label;
            return span;
        }),
    );
    const missing = Object.entries(health.fieldsMissing);
    const box = $('h-missing');
    if (missing.length) {
        const ul = document.createElement('ul');
        for (const [id, fields] of missing) {
            const li = document.createElement('li');
            li.textContent = `${id} lacks `;
            const span = document.createElement('span');
            span.textContent = fields.join(', ');
            li.append(span);
            ul.append(li);
        }
        box.replaceChildren(ul);
    } else {
        box.replaceChildren();
    }
    show(section);
}

/** Recompute the open manifest's health from its scans and push it to the
 * drill-down section and the card. */
function updateHealth(): void {
    const health = scans.size ? healthFromScans(manifest, scans) : null;
    const card = currentCard();
    if (card) card.health = health;
    renderHealthSection(health);
    renderOverview();
}

// --- staging / preview / summary ------------------------------------------- //

function entryFiles(source: SourceSpec): { relPath: string }[] {
    if (!workspace) throw new Error('no workspace');
    const all = workspace.files.map((f) => f.relPath);
    let wanted: string[];
    if (source.Pattern) {
        wanted = all.filter((p) => globMatch(source.Pattern ?? '', p));
        if (!wanted.length) {
            throw new Error(
                `meshio++: dataset: the pattern '${source.Pattern}' matches nothing in the workspace`,
            );
        }
    } else if (source.Path) {
        if (!all.includes(source.Path)) {
            throw new Error(`meshio++: dataset: no workspace file '${source.Path}'`);
        }
        wanted = [source.Path];
    } else {
        wanted = source.Paths ?? [];
        for (const p of wanted) {
            if (!all.includes(p)) {
                throw new Error(`meshio++: dataset: no workspace file '${p}'`);
            }
        }
    }
    return wanted.map((relPath) => ({ relPath }));
}

async function stageSelected(): Promise<PlanEntry[]> {
    if (!workspace || !selected) throw new Error('nothing selected');
    const entry = getEntry(manifest, selected);
    const byPath = new Map(workspace.files.map((f) => [f.relPath, f]));
    const files = await Promise.all(
        entryFiles(entry.source).map(async ({ relPath }) => ({
            relPath,
            bytes: await byPath.get(relPath)!.bytes(),
        })),
    );
    const result = await client.stageEntry(files, entry.source);
    staged = { id: entry.id, plan: result.plan };
    return result.plan;
}

async function previewStep(step: number): Promise<void> {
    if (!staged) return;
    setStatus('loading', `step ${step + 1}/${staged.plan.length}…`);
    const render = await client.previewStep(step);
    show3d(render.vtp);
    const plan = staged.plan[step];
    $('scrub-label').textContent =
        `step ${step + 1}/${staged.plan.length} · t = ${plan?.time ?? 0} (${plan?.timeSource ?? '?'})`;
    setState({ step });
    setStatus('ready');
    void refreshSummary(step);
}

async function refreshSummary(step: number): Promise<void> {
    try {
        const { arrays: summaries } = await client.summarizeStep(step);
        renderSummary(summaries);
        setState({
            summary: summaries.map((s) => ({
                name: s.name,
                min: s.min,
                mean: s.mean,
                numNan: s.numNan,
                numInf: s.numInf,
            })),
        });
    } catch {
        // A format the flat Mesh cannot express is a missing summary, never
        // a failed preview.
        $('summary-wrap').hidden = true;
        setState({ summary: null });
    }
}

function renderSummary(summaries: ArraySummary[]): void {
    const wrap = $('summary-wrap');
    if (!summaries.length) {
        wrap.hidden = true;
        return;
    }
    const table = $('summary');
    const header = document.createElement('tr');
    for (const label of ['array', 'dtype', 'min', 'max', 'mean', 'NaN', 'Inf']) {
        const th = document.createElement('th');
        th.textContent = label;
        header.append(th);
    }
    const fmt = (v: number) =>
        Number.isFinite(v) ? Number(v.toPrecision(4)).toString() : '—';
    // Data rows first, quality metrics grouped at the end.
    const ordered = [...summaries].sort(
        (a, b) => Number(isQualityRow(a)) - Number(isQualityRow(b)),
    );
    const rows = ordered.map((s) => {
        const tr = document.createElement('tr');
        const quality = isQualityRow(s);
        if (quality) tr.classList.add('quality');
        const bad = quality ? qualityRowIsBad(s) : s.numNan + s.numInf > 0;
        if (bad) tr.classList.add('bad');
        const cells = [
            `${s.name} (${s.location})`,
            s.dtype,
            fmt(s.min),
            fmt(s.max),
            fmt(s.mean),
            String(s.numNan),
            String(s.numInf),
        ];
        for (const text of cells) {
            const td = document.createElement('td');
            td.textContent = text;
            tr.append(td);
        }
        return tr;
    });
    table.replaceChildren(header, ...rows);
    show(wrap);
}

/** Give the open manifest's card a thumbnail from the view just rendered,
 * best-effort (a card only ever lacks one, never shows a wrong one). */
async function thumbnailFromView(): Promise<void> {
    const card = currentCard();
    if (!card) return;
    try {
        card.thumbnail = await captureThumbnail(renderer);
        renderOverview();
    } catch {
        // WebGL capture unavailable (e.g. a context without a back buffer)
    }
}

async function previewSelected(): Promise<void> {
    try {
        setStatus('loading', 'staging…');
        const plan = await stageSelected();
        setState({ planLength: plan.length });
        const scrub = $<HTMLInputElement>('scrub');
        scrub.max = String(Math.max(plan.length - 1, 0));
        scrub.value = '0';
        $('scrub-row').hidden = plan.length < 2;
        await previewStep(0);
        if (!currentCard()?.thumbnail) await thumbnailFromView();
    } catch (e) {
        fail(e);
    }
}

/** The scan cache key: the manifest, the entry and the newest modification
 * time of its source files — a source that changes invalidates the scan.
 * Null when a backend cannot report modification times (no caching then). */
async function scanCacheKey(entry: ManifestEntry): Promise<string | null> {
    if (!workspace || !manifestPath) return null;
    const byPath = new Map(workspace.files.map((f) => [f.relPath, f]));
    let newest = 0;
    for (const { relPath } of entryFiles(entry.source)) {
        const mtime = await byPath.get(relPath)?.lastModified();
        if (mtime === null || mtime === undefined) return null;
        newest = Math.max(newest, mtime);
    }
    return `scan:${manifestPath}:${entry.id}:${newest}`;
}

function isEntryScan(value: unknown): value is EntryScan {
    return (
        typeof value === 'object' &&
        value !== null &&
        Array.isArray((value as EntryScan).arrays) &&
        typeof (value as EntryScan).steps === 'number'
    );
}

/** Pick up cached scans for the open manifest (from an earlier visit) so its
 * health shows without re-staging every entry. */
async function restoreScansFromCache(): Promise<void> {
    for (const entry of manifest.entries) {
        try {
            const key = await scanCacheKey(entry);
            const cached = key ? await getValue<unknown>(key) : null;
            if (isEntryScan(cached)) scans.set(entry.id, cached);
        } catch {
            // a source that does not resolve is simply not cached
        }
    }
    if (scans.size) {
        setState({ scans: Object.fromEntries(scans) });
        updateHealth();
        renderEntryList();
    }
}

async function scanEntries(): Promise<void> {
    if (!workspace) return;
    const previous = selected;
    let stagedAny = false;
    for (const entry of manifest.entries) {
        try {
            setStatus('loading', `scanning ${entry.id}…`);
            const key = await scanCacheKey(entry);
            const cached = key ? await getValue<unknown>(key) : null;
            if (isEntryScan(cached)) {
                scans.set(entry.id, cached);
                continue;
            }
            selected = entry.id;
            const plan = await stageSelected();
            stagedAny = true;
            const { arrays: summaries } = await client.summarizeStep(0);
            const scan = entryScanFromSummaries(plan.length, summaries);
            scans.set(entry.id, scan);
            if (key) await putValue(key, scan);
        } catch {
            scans.set(entry.id, emptyScan());
        }
    }
    if (stagedAny) await client.evictEntry();
    staged = null;
    selected = previous;
    setState({ scans: Object.fromEntries(scans) });
    updateHealth();
    setStatus('ready', 'scan complete');
    renderEntryList();
}

/** Scan one card's manifest from the overview (opening it first). */
async function scanManifestPath(path: string): Promise<void> {
    await openManifestPath(path);
    await scanEntries();
}

/** Open, scan and (optionally) thumbnail every manifest in turn — strictly
 * one staged entry at a time, the MEMFS rule — then return to the overview. */
async function scanAllManifests(): Promise<void> {
    if (!workspace) return;
    if (window.__datasetState.dirty) {
        setStatus('error', 'save or discard the unsaved manifest edits before scanning all');
        return;
    }
    const withThumbs = $<HTMLInputElement>('ov-thumbs').checked;
    bulkScan = true;
    try {
        for (const card of [...cards]) {
            if (card.parseError) continue;
            await openManifestPath(card.path);
            await scanEntries();
            if (withThumbs && !card.thumbnail && manifest.entries.length) {
                selectEntry(manifest.entries[0]?.id ?? null);
                await previewSelected();
                await thumbnailFromView();
                await client.evictEntry();
                staged = null;
            }
        }
    } finally {
        bulkScan = false;
    }
    selectEntry(null);
    setView('overview');
    renderOverview();
    setStatus('ready', `scanned ${cards.filter((c) => !c.parseError).length} manifest(s)`);
}

// --- companion process ------------------------------------------------------ //

const SERVER_KEY = 'server';

function setServer(patch: Partial<ServerState> & { url: string }): void {
    const previous = window.__datasetState.server;
    setState({
        server: {
            url: patch.url,
            connected: patch.connected ?? previous?.connected ?? false,
            version: patch.version ?? previous?.version ?? null,
            tools: patch.tools ?? previous?.tools ?? [],
            root: patch.root ?? previous?.root ?? null,
            runsDir: patch.runsDir ?? previous?.runsDir ?? null,
            error: patch.error ?? null,
        },
    });
}

/** Show/hide every `data-needs-server` control and the status line. */
function updateServerUi(): void {
    const server = window.__datasetState.server;
    const connected = !!server?.connected;
    for (const el of document.querySelectorAll<HTMLElement>('[data-needs-server]')) {
        el.hidden = !connected;
    }
    show($('srv-disconnect'), connected);
    show($('srv-connect'), !connected);
    const status = $('srv-status');
    status.className = 'hint';
    if (!server) {
        status.textContent = '';
    } else if (connected) {
        status.classList.add('connected');
        status.textContent =
            `connected — meshio++ ${server.version ?? '?'}` +
            (server.root ? `, root ${server.root}` : ', unrestricted paths');
    } else {
        status.classList.add('failed');
        status.textContent = server.error ?? 'not connected';
    }
    $<HTMLButtonElement>('e-scan-server').disabled = !currentCard()?.serverPath;
    trainPanel.refreshAvailability();
}

/** Bind workspace cards to the server's manifests by content hash (FSA
 * never reveals an absolute path, so bytes are the only shared identity),
 * and add a card for every manifest only the server knows. */
function mergeServerCards(): void {
    cards = cards.filter((c) => !c.serverOnly);
    for (const card of cards) card.serverPath = null;
    for (const m of serverManifests) {
        const bound = cards.find((c) => c.sha256 !== null && c.sha256 === m.sha256);
        if (bound) {
            bound.serverPath = m.path;
        } else {
            cards.push(cardFromServerManifest(m));
        }
    }
}

async function refreshServerManifests(): Promise<void> {
    if (!serverApi) return;
    const found = await serverApi.tool<{ manifests: ServerManifest[] }>('dataset_find', {
        root_dir: '.',
        max_depth: 2,
    });
    if (isToolError(found)) throw new Error(found.error);
    serverManifests = found.manifests;
    mergeServerCards();
    renderOverview();
    updateServerUi();
}

async function connectServer(url: string, token: string): Promise<void> {
    const api = new ApiClient(url, token || null);
    setStatus('loading', 'connecting…');
    try {
        const health = await api.health();
        serverApi = api;
        setServer({
            url: api.baseUrl,
            connected: true,
            version: health.version,
            tools: health.tools,
            root: health.root,
            runsDir: health.runs_dir,
            error: null,
        });
        await putValue(SERVER_KEY, { url: api.baseUrl, token: token || '' });
        updateServerUi();
        setStatus('ready', `connected to ${api.baseUrl}`);
        await refreshServerManifests();
    } catch (e) {
        serverApi = null;
        setServer({
            url,
            connected: false,
            error: e instanceof Error ? e.message : String(e),
        });
        updateServerUi();
        setStatus('ready');
    }
}

function disconnectServer(): void {
    trainPanel.dispose();
    serverApi = null;
    serverManifests = [];
    setServer({ url: $<HTMLInputElement>('srv-url').value, connected: false, error: null });
    mergeServerCards();
    renderOverview();
    updateServerUi();
}

/** Health from the server producer, onto a card and — when that card is
 * the open manifest — onto its entries, replacing the browser scans. */
async function scanOnServer(path: string): Promise<void> {
    const card = cards.find((c) => c.path === path);
    if (!serverApi || !card?.serverPath) return;
    setStatus('loading', `scanning ${card.name ?? card.path} on the server…`);
    const report = await serverApi.tool<ServerHealthReport>('dataset_health', {
        manifest_path: card.serverPath,
    });
    if (isToolError(report)) throw new Error(report.error);
    const { health, scans: serverScans } = fromServerHealth(report);
    card.health = health;
    if (!card.serverOnly && card.path === manifestPath) {
        scans.clear();
        for (const [id, scan] of Object.entries(serverScans)) scans.set(id, scan);
        setState({ scans: Object.fromEntries(scans) });
        renderHealthSection(health);
        renderEntryList();
    }
    renderOverview();
    setStatus('ready', `server scan of ${card.name ?? card.path} complete`);
}

/** Try the remembered server on boot, silently: a stale entry just leaves
 * the section disconnected with its fields prefilled. */
async function restoreServer(): Promise<void> {
    const saved = await getValue<{ url: string; token: string }>(SERVER_KEY);
    if (!saved || typeof saved.url !== 'string') return;
    $<HTMLInputElement>('srv-url').value = saved.url;
    $<HTMLInputElement>('srv-token').value = saved.token ?? '';
    await connectServer(saved.url, saved.token ?? '');
}

// --- training --------------------------------------------------------------- //

const trainPanel = new TrainPanel({
    getApi: () => serverApi,
    // Training runs on the server, against the file the server can see: the
    // manifest's own path there, which a card only has once its bytes match
    // (so an unsaved edit correctly disables the form).
    getManifestPath: () => currentCard()?.serverPath ?? null,
    setState,
    setStatus: (message) => setStatus('ready', message),
    fail,
});

// --- add-case flow ---------------------------------------------------------- //

function openAddFlow(): void {
    if (!workspace) return;
    const candidates = workspace.files
        .map((f) => f.relPath)
        .filter((p) => {
            const ext = p.slice(p.lastIndexOf('.') + 1).toLowerCase();
            return readableExtensions.size === 0 || readableExtensions.has(ext);
        });
    $('add-files').replaceChildren(
        ...candidates.map((relPath) => {
            const li = document.createElement('li');
            const label = document.createElement('label');
            const box = document.createElement('input');
            box.type = 'checkbox';
            box.value = relPath;
            box.addEventListener('change', updateSuggestion);
            const span = document.createElement('span');
            span.textContent = relPath;
            label.append(box, span);
            li.append(label);
            return li;
        }),
    );
    $('add-error').hidden = true;
    $('add-suggest').textContent = '';
    $<HTMLInputElement>('add-id').value = '';
    show($('add-flow'));
}

function selectedAddFiles(): string[] {
    return [...$('add-files').querySelectorAll<HTMLInputElement>('input:checked')].map(
        (box) => box.value,
    );
}

let addSuggestion: SourceSpec | null = null;

function updateSuggestion(): void {
    const picked = selectedAddFiles();
    if (!picked.length || !workspace) {
        addSuggestion = null;
        $('add-suggest').textContent = '';
        return;
    }
    addSuggestion = suggestSource(
        picked,
        workspace.files.map((f) => f.relPath),
    );
    $('add-suggest').textContent = `source: ${describeSource(addSuggestion)}`;
    try {
        $<HTMLInputElement>('add-id').placeholder = deriveId(
            addSuggestion,
            manifest.entries.map((e) => e.id),
        );
    } catch {
        $<HTMLInputElement>('add-id').placeholder = 'id required';
    }
}

async function commitAdd(): Promise<void> {
    if (!addSuggestion) return;
    const idText = $<HTMLInputElement>('add-id').value.trim();
    try {
        const entry = addEntry(manifest, addSuggestion, {
            id: idText || undefined,
        });
        // Validate now, like `add(validate_source=True)`: stage and resolve
        // the plan so a bad case fails here, by name, not at training time.
        selectEntry(entry.id);
        await previewSelected();
        $('add-flow').hidden = true;
        refreshCurrentCard();
    } catch (e) {
        // Roll back a half-added entry so the manifest never holds an entry
        // whose source did not resolve.
        const idAdded = manifest.entries.at(-1)?.id;
        if (idAdded && !savedTextHas(idAdded)) {
            try {
                removeEntry(manifest, idAdded);
            } catch {
                // already absent
            }
        }
        selectEntry(null);
        const box = $('add-error');
        box.textContent = e instanceof Error ? e.message : String(e);
        box.hidden = false;
    }
}

function savedTextHas(id: string): boolean {
    return savedText !== null && savedText.includes(`"Id": ${JSON.stringify(id)}`);
}

// --- manifest open/save ----------------------------------------------------- //

// Sentinel select-value for "(new) dataset.json" — cannot collide with a
// candidate, since every candidate relPath ends in `.json`.
const NEW_MANIFEST = '::new::';

function download(text: string, filename: string): void {
    const url = URL.createObjectURL(new Blob([text], { type: 'application/json' }));
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    a.click();
    URL.revokeObjectURL(url);
}

function populateManifestSelect(): void {
    if (!workspace) return;
    const candidates = workspace.manifestCandidates();
    setOptions(
        $<HTMLSelectElement>('manifest-file'),
        [
            ...candidates.map((p) => ({ value: p, label: p })),
            { value: NEW_MANIFEST, label: '(new) dataset.json' },
        ],
        candidates[0] ?? NEW_MANIFEST,
    );
    show($('manifest-section'));
}

async function openManifest(): Promise<void> {
    if (!workspace) return;
    const choice = $<HTMLSelectElement>('manifest-file').value;
    try {
        if (choice === NEW_MANIFEST) {
            manifest = emptyManifest();
            manifestPath = 'dataset.json';
            savedText = null;
        } else {
            const text = await workspace.readText(choice);
            manifest = parseManifest(text);
            manifestPath = choice;
            savedText = text;
        }
        scans.clear();
        setState({ scans: {} });
        selectEntry(null);
        $<HTMLInputElement>('m-name').value = manifest.name ?? '';
        $<HTMLInputElement>('m-desc').value = manifest.description ?? '';
        setState({ manifestName: manifestPath });
        show($('manifest-editor'));
        show($('entries-section'));
        renderHealthSection(currentCard()?.health ?? null);
        renderEntryList();
        updateServerUi();
        if (!bulkScan) setView('manifest');
        setStatus('ready', `${manifest.entries.length} entr${manifest.entries.length === 1 ? 'y' : 'ies'}`);
        await restoreScansFromCache();
    } catch (e) {
        fail(e);
    }
}

/** Open a manifest by path (a card click). Re-entering the manifest that is
 * already open with unsaved edits keeps those edits rather than re-reading
 * the file. */
async function openManifestPath(path: string): Promise<void> {
    if (path === manifestPath && window.__datasetState.dirty) {
        if (!bulkScan) setView('manifest');
        return;
    }
    const select = $<HTMLSelectElement>('manifest-file');
    if ([...select.options].some((o) => o.value === path)) select.value = path;
    else select.value = NEW_MANIFEST;
    await openManifest();
}

async function saveManifest(): Promise<void> {
    if (!workspace || !manifestPath) return;
    try {
        const text = currentText();
        const mode = await workspace.saveManifest(manifestPath, text);
        savedText = text;
        setState({ lastSave: { mode, text } });
        refreshDirty();
        if (mode === 'in-place') {
            // The file list may have gained the manifest; refresh the select
            // and the card set from disk (health/thumbnail survive by hash).
            populateManifestSelect();
            $<HTMLSelectElement>('manifest-file').value = manifestPath;
            const card = currentCard();
            const hash = await sha256Hex(text);
            if (card) {
                card.lastModified = Date.now();
                card.sha256 = hash;
            }
        }
        refreshCurrentCard();
        setStatus(
            'ready',
            mode === 'in-place' ? `saved ${manifestPath}` : `downloaded ${manifestPath}`,
        );
    } catch (e) {
        fail(e);
    }
}

// --- diff ------------------------------------------------------------------- //

const DIFF_DISK = '::disk::';
const DIFF_CURRENT = '::current::';
const DIFF_PASTE = '::paste::';

function diffChoices(): { value: string; label: string }[] {
    const choices: { value: string; label: string }[] = [];
    if (manifestPath && savedText !== null) {
        choices.push({ value: DIFF_DISK, label: `${manifestPath} (on disk)` });
    }
    if (manifestPath) {
        choices.push({ value: DIFF_CURRENT, label: `${manifestPath} (current edits)` });
    }
    for (const p of workspace?.manifestCandidates() ?? []) {
        choices.push({ value: p, label: p });
    }
    choices.push({ value: DIFF_PASTE, label: 'pasted JSON' });
    return choices;
}

async function diffText(choice: string): Promise<{ label: string; text: string }> {
    if (choice === DIFF_DISK) return { label: `${manifestPath} (on disk)`, text: savedText ?? '' };
    if (choice === DIFF_CURRENT) {
        return { label: `${manifestPath} (current edits)`, text: currentText() };
    }
    if (choice === DIFF_PASTE) {
        return { label: 'pasted JSON', text: $<HTMLTextAreaElement>('diff-paste').value };
    }
    if (!workspace) throw new Error('no workspace');
    return { label: choice, text: await workspace.readText(choice) };
}

async function runDiff(): Promise<void> {
    try {
        const [a, b] = await Promise.all([
            diffText($<HTMLSelectElement>('diff-a').value),
            diffText($<HTMLSelectElement>('diff-b').value),
        ]);
        const diff = diffManifests(parseManifest(a.text), parseManifest(b.text));
        const lines = diffLines(a.text, b.text);
        renderDiff($('diff-body'), diff, lines);
        const summary = summarizeDiff(a.label, b.label, diff);
        $('diff-summary').textContent =
            `${summary.a} → ${summary.b}: +${summary.added} added, ` +
            `−${summary.removed} removed, ~${summary.changed} changed` +
            (summary.headerChanged ? `, ${summary.headerChanged} document field(s)` : '');
        setState({ diff: summary });
    } catch (e) {
        $('diff-body').replaceChildren();
        $('diff-summary').textContent = e instanceof Error ? e.message : String(e);
        setState({ diff: null });
    }
}

function openDiff(): void {
    const choices = diffChoices();
    const a = savedText !== null && manifestPath ? DIFF_DISK : (choices[0]?.value ?? DIFF_PASTE);
    const b = manifestPath
        ? DIFF_CURRENT
        : (choices.find((c) => c.value !== a && c.value !== DIFF_PASTE)?.value ?? DIFF_PASTE);
    setOptions($<HTMLSelectElement>('diff-a'), choices, a);
    setOptions($<HTMLSelectElement>('diff-b'), choices, b);
    show($('diff-wrap'));
    void runDiff();
}

// --- workspace pick --------------------------------------------------------- //

function adoptWorkspace(ws: Workspace): void {
    workspace = ws;
    setState({ backendMode: ws.kind });
    if (ws.kind === 'fsa' && ws.root) {
        // Best-effort: a failure (mock handle, no IndexedDB, private mode)
        // degrades to picking again next visit, never to a broken adopt.
        void saveHandle(ws.root);
    }
    $('ws-info').textContent =
        `${ws.files.length} file(s) — ` +
        (ws.kind === 'fsa'
            ? 'the manifest saves back in place.'
            : 'this browser cannot write back; saving downloads the manifest.');
    show($('ws-info'));
    populateManifestSelect();
    show($('overview-tools'));
    setView('overview');
    setStatus('idle', 'workspace ready — open or create a manifest');
    void buildOverview().catch(fail);
}

function pickFlow(): void {
    if (hasFsAccess()) {
        pickWorkspaceFsa().then(adoptWorkspace).catch(fail);
    } else {
        $<HTMLInputElement>('ws-input').click();
    }
}

/** Offer "Reopen" when a handle survived in IndexedDB. Restore itself runs
 * behind the click — the persisted handle's permission is 'prompt' and
 * requestPermission demands a user gesture, so nothing restores on boot. */
function offerRestore(): void {
    void loadHandle().then((handle) => {
        if (!handle) return;
        const button = $('ws-reopen');
        button.textContent = `Reopen “${handle.name}”`;
        show(button);
        setState({ restoreAvailable: true });
        button.addEventListener('click', () => {
            workspaceFromHandle(handle)
                .then(adoptWorkspace)
                .catch(async () => {
                    // Denied or stale: forget it and fall back to a pick.
                    await clearHandle();
                    button.hidden = true;
                    setState({ restoreAvailable: false });
                    pickFlow();
                });
        });
    });
}

// --- boot ------------------------------------------------------------------- //

async function boot(): Promise<void> {
    try {
        const { formats } = await client.init();
        readableExtensions = new Set<string>();
        // The add-case list filters by readable extension; format keys mostly
        // ARE extensions in this registry, which is close enough for a filter
        // that the user can always override by picking anyway.
        for (const f of formats.readers) readableExtensions.add(f.toLowerCase());
        for (const extra of ['msh', 'e', 'exo', 'h5m', 'foam']) {
            readableExtensions.add(extra);
        }
        setStatus('idle', 'choose a case directory to begin');
    } catch (e) {
        fail(e);
        return;
    }

    $('ws-pick').addEventListener('click', pickFlow);
    offerRestore();
    $<HTMLInputElement>('ws-input').addEventListener('change', (e) => {
        const input = e.target as HTMLInputElement;
        if (input.files?.length) {
            adoptWorkspace(workspaceFromFileList(input.files, download));
        }
    });

    $('ov-sort').addEventListener('change', renderOverview);
    $('ov-scan-all').addEventListener('click', () => void scanAllManifests().catch(fail));
    $('ov-diff').addEventListener('click', openDiff);
    $('ov-new').addEventListener('click', () => {
        $<HTMLSelectElement>('manifest-file').value = NEW_MANIFEST;
        void openManifest();
    });
    $('m-back').addEventListener('click', () => {
        refreshCurrentCard();
        setView('overview');
    });

    $('srv-connect').addEventListener('click', () => {
        void connectServer(
            $<HTMLInputElement>('srv-url').value.trim(),
            $<HTMLInputElement>('srv-token').value.trim(),
        );
    });
    $('srv-disconnect').addEventListener('click', disconnectServer);
    $('e-scan-server').addEventListener('click', () => {
        const path = currentCard()?.path;
        if (path) void scanOnServer(path).catch(fail);
    });
    void restoreServer();

    $('manifest-open').addEventListener('click', () => void openManifest());
    $('m-save').addEventListener('click', () => void saveManifest());
    $<HTMLInputElement>('m-name').addEventListener('change', (e) => {
        manifest.name = (e.target as HTMLInputElement).value.trim() || null;
        refreshDirty();
    });
    $<HTMLInputElement>('m-desc').addEventListener('change', (e) => {
        manifest.description = (e.target as HTMLInputElement).value.trim() || null;
        refreshDirty();
    });

    $<HTMLInputElement>('e-filter').addEventListener('input', renderEntryList);
    $('e-add').addEventListener('click', openAddFlow);
    $('e-scan').addEventListener('click', () => void scanEntries());
    $('add-commit').addEventListener('click', () => void commitAdd());
    $('add-cancel').addEventListener('click', () => {
        $('add-flow').hidden = true;
    });

    for (const id of ['d-split', 'd-tags', 'd-group', 'd-notes', 'd-meta']) {
        $(id).addEventListener('change', commitDetail);
    }
    $('d-preview').addEventListener('click', () => void previewSelected());
    $('d-remove').addEventListener('click', () => {
        if (!selected) return;
        removeEntry(manifest, selected);
        scans.delete(selected);
        selectEntry(null);
        updateHealth();
    });

    $('diff-a').addEventListener('change', () => void runDiff());
    $('diff-b').addEventListener('change', () => void runDiff());
    $('diff-paste').addEventListener('input', () => void runDiff());
    $('diff-close').addEventListener('click', () => {
        $('diff-wrap').hidden = true;
    });

    $<HTMLInputElement>('scrub').addEventListener('input', (e) => {
        const step = Number((e.target as HTMLInputElement).value);
        clearTimeout(scrubTimer);
        scrubTimer = window.setTimeout(() => void previewStep(step).catch(fail), 120);
    });
    $<HTMLSelectElement>('color-by').addEventListener('change', (e) => {
        applyColorBy((e.target as HTMLSelectElement).value);
    });
    $('tool-reset').addEventListener('click', () => renderer.resetCamera());
}

let scrubTimer = 0;

void boot();
