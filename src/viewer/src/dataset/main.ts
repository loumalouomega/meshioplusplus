/**
 * The dataset-manager page (`dataset.html`) — build and curate a
 * `DatasetManifest` JSON against a local case directory, previewing entries
 * through the same worker/renderer the viewer uses.
 *
 * Structure mirrors `src/main.ts`: static markup in the HTML, `$()` lookups
 * that throw on a missing id, one state object (`window.__datasetState`)
 * that the Playwright spec asserts on. All document logic lives in the pure
 * modules (`manifest.ts`, `glob.ts`, `suggest.ts`); all file access in
 * `fs.ts`; everything wasm-side behind the shared worker.
 */

import { Legend } from '../render/legend';
import { Renderer, SOLID_COLOR, colorKey } from '../render/renderer';
import { $, setOptions, show } from '../ui/dom';
import type { ArrayEntry } from '../types';
import * as client from '../worker/client';
import type { ArraySummary, PlanEntry } from '../worker/protocol';
import {
    addEntry,
    annotateEntry,
    emptyManifest,
    getEntry,
    parseManifest,
    removeEntry,
    stringifyManifest,
    type Manifest,
    type ManifestEntry,
    type SourceSpec,
} from './manifest';
import { globMatch } from './glob';
import { suggestSource } from './suggest';
import { deriveId } from './manifest';
import {
    hasFsAccess,
    pickWorkspaceFsa,
    workspaceFromFileList,
    workspaceFromHandle,
    type Workspace,
} from './fs';
import { clearHandle, loadHandle, saveHandle } from './persist';
import type { DatasetState, EntryScan } from './types';

// --- state ----------------------------------------------------------------- //

window.__datasetState = {
    status: 'booting',
    backendMode: null,
    restoreAvailable: false,
    manifestName: null,
    entryIds: [],
    selected: null,
    planLength: 0,
    step: 0,
    numPoints: 0,
    numCells: 0,
    summary: null,
    scans: {},
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

function refreshDirty(): void {
    setState({ dirty: savedText !== null && currentText() !== savedText });
    $('m-dirty').hidden = !window.__datasetState.dirty;
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

/** Is a summary row a `quality:*` metric? Their NaN means "metric N/A for
 * this cell type" BY DESIGN, so quality rows are excluded from every
 * NaN-based bad-case rule — here and in `scanEntries` — or every entry
 * would badge spuriously. */
function isQualityRow(s: ArraySummary): boolean {
    return s.name.startsWith('quality:');
}

/** The quality-specific bad rules: inverted/degenerate cells present, or a
 * negative worst scaled Jacobian. */
function qualityRowIsBad(s: ArraySummary): boolean {
    if (s.name === 'quality:inverted' || s.name === 'quality:degenerate') {
        return s.mean * s.numValues > 0;
    }
    return s.name === 'quality:scaled_jacobian' && s.min < 0;
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
    } catch (e) {
        fail(e);
    }
}

async function scanEntries(): Promise<void> {
    if (!workspace) return;
    const previous = selected;
    for (const entry of manifest.entries) {
        try {
            setStatus('loading', `scanning ${entry.id}…`);
            selected = entry.id;
            const plan = await stageSelected();
            const { arrays: summaries } = await client.summarizeStep(0);
            // The NaN/Inf lane counts DATA arrays only — a quality metric's
            // NaN is "N/A for this cell type" by design, not a bad value.
            const dataRows = summaries.filter((s) => !isQualityRow(s));
            const inverted = summaries.find((s) => s.name === 'quality:inverted');
            const jacobian = summaries.find(
                (s) => s.name === 'quality:scaled_jacobian',
            );
            scans.set(entry.id, {
                steps: plan.length,
                numNan: dataRows.reduce((n, s) => n + s.numNan, 0),
                numInf: dataRows.reduce((n, s) => n + s.numInf, 0),
                numInverted: inverted
                    ? Math.round(inverted.mean * inverted.numValues)
                    : 0,
                minScaledJacobian:
                    jacobian && Number.isFinite(jacobian.min) ? jacobian.min : null,
            });
        } catch {
            scans.set(entry.id, {
                steps: 0,
                numNan: 0,
                numInf: 0,
                numInverted: 0,
                minScaledJacobian: null,
            });
        }
    }
    await client.evictEntry();
    staged = null;
    selected = previous;
    setState({ scans: Object.fromEntries(scans) });
    setStatus('ready', 'scan complete');
    renderEntryList();
}

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
        selectEntry(null);
        $<HTMLInputElement>('m-name').value = manifest.name ?? '';
        $<HTMLInputElement>('m-desc').value = manifest.description ?? '';
        setState({ manifestName: manifestPath });
        show($('manifest-editor'));
        show($('entries-section'));
        renderEntryList();
        setStatus('ready', `${manifest.entries.length} entr${manifest.entries.length === 1 ? 'y' : 'ies'}`);
    } catch (e) {
        fail(e);
    }
}

async function saveManifest(): Promise<void> {
    if (!workspace || !manifestPath) return;
    try {
        const text = currentText();
        const mode = await workspace.saveManifest(manifestPath, text);
        savedText = text;
        setState({ lastSave: { mode, text } });
        refreshDirty();
        setStatus(
            'ready',
            mode === 'in-place' ? `saved ${manifestPath}` : `downloaded ${manifestPath}`,
        );
    } catch (e) {
        fail(e);
    }
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
    setStatus('idle', 'workspace ready — open or create a manifest');
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
