/**
 * Viewer entry point.
 *
 * Resolves where the mesh comes from, wires the UI to the renderer, and — in
 * the full build only — routes files and operations through the meshio++ WASM
 * worker.
 *
 * Three payload sources, tried in order:
 *   1. an embedded base64 VTP in `#vtp-payload`, written by
 *      `meshioplusplus.view(backend="browser")`;
 *   2. a `?url=` query parameter, for shareable links;
 *   3. the file picker / drag-and-drop, which needs the WASM worker.
 *
 * Under `__VIEWER_EMBEDDED__` only (1) is compiled in, so Rollup drops the
 * worker and the WASM binary from the wheel-bundled build.
 */
import { Legend } from './render/legend';
import { initPrefs } from './features/prefs';
import { type AnalysisController, initAnalysis } from './features/analysis';
import { type ViewController, initView } from './features/view';
import { setupOrientationCube } from './render/orientation';
import { Renderer, SOLID_COLOR, colorKey } from './render/renderer';
import type { RepresentationMode } from './render/renderer';
import { Picker } from './render/pick';
import { downloadPng } from './render/screenshot';
import { OpsPanel } from './ops/panel';
import type { RefineTarget } from './ops/panel';
import { ICONS } from './ui/icons';
import { glyph } from './ui/glyphs';
import { $, maybe, setOptions, show } from './ui/dom';
import { initTheme, isThemeMode } from './ui/theme';
import { initSections, setCollapsed } from './ui/section';
import { initMenus, setMenuChecked } from './ui/menu';
import { findByShortcut, registerAction } from './app/actions';
import { initSidebar, showShortcuts } from './features/shell';
import { setStatusBar } from './features/statusbar';
import { SelectionController } from './features/selection';
import { GroupsController } from './features/groups';
import { Visibility } from './app/visibility';
import type { RegionsResult } from './worker/protocol';
import { SCENE_THEMES } from './ui/theme';
import { opsPreserveCellIds } from './worker/opsIds';
import type { ArrayEntry, PickInfo, ScalarRange, ViewerStatus } from './types';
import type { SurfaceBuffers } from '@meshioplusplus/wasm';

import type { MeshMeta, OpSpec, RenderResult } from './worker/protocol';

const renderer = new Renderer($('render'));
const theme = initTheme();
theme.onChange((_resolved, scene) => renderer.setSceneTheme(scene.background, scene.body));
const opsPanel = new OpsPanel({
    onChange: (ops) => applyOps(ops),
    refineTargets: () => refineTargets(),
});
let arrays: ArrayEntry[] = [];
/** True when the open file has no volume cells and no vertex or line blocks. */
let surfaceOnlyPolygons = false;
let currentFile: File | null = null;

/**
 * The Playwright hook. Deliberately a documented part of the app: asserting on
 * this instead of on pixels keeps the end-to-end tests GPU-independent.
 */
window.__viewerState = {
    status: 'booting',
    numPoints: 0,
    numCells: 0,
    arrays: [],
    name: '',
    picked: null,
    selection: { tool: 'off', active: null, sets: [], hidden: 0 },
    groups: { status: 'none', groups: [] },
    analysis: {
        table: { open: false, kind: 'cells', rows: 0, sortedBy: null, direction: 'asc' },
        histogram: { array: null, counts: [] },
    },
    view: {
        parallel: false,
        direction: [0, 0, 1],
        logScale: false,
        bands: 0,
        clip: false,
        deformed: 0,
        measuring: false,
        measured: null,
        bookmarks: [],
    },
};

function setStatus(status: ViewerStatus, message = ''): void {
    window.__viewerState.status = status;
    const el = $('status');
    el.dataset.status = status;
    el.textContent = message;
    el.hidden = !message;
}

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

// --- display ------------------------------------------------------------- //

/** A worker's surface, or the VTP the offline embedded page carries. */
function show3d(
    source: SurfaceBuffers | ArrayBuffer,
    meta: MeshMeta | null,
    name: string,
    keepSelection = false,
    idsPreserved = true
): void {
    const info = source instanceof ArrayBuffer ? renderer.load(source) : renderer.loadSurface(source);
    arrays = info.arrays;

    Object.assign(window.__viewerState, {
        numPoints: info.numPoints,
        numCells: info.numCells,
        arrays: arrays.map((a) => a.label),
        name,
    });

    fillColorBy(SOLID_COLOR);
    show($('color-section'));
    legend.hide();
    picker.clear();
    renderPick(null);
    if (keepSelection) visibility.pruneTo(info.numCells);
    else visibility.reset();
    selection.attachSurface(keepSelection);
    groups.attachSurface(info.numCells, keepSelection, idsPreserved);
    syncConvertGroups();
    view.attachSurface();
    analysis.attachSurface();
    surfaceOnlyPolygons =
        !!meta &&
        meta.cellBlocks.length > 0 &&
        meta.cellBlocks.every((b) => /^(triangle|quad|polygon)/.test(b.type));
    renderInfo(meta, info.numCells, name);
    setStatusBar({ file: name, points: info.numPoints, cells: info.numCells, selected: 0 });
    setStatus('ready');
}

function renderInfo(meta: MeshMeta | null, renderedCells: number, name: string): void {
    const rows: [string, string][] = [];
    const n = (v: number) => v.toLocaleString('en-US');
    if (name) rows.push(['File', name]);
    if (meta) {
        rows.push(['Format', meta.format]);
        rows.push(['Points', n(meta.numPoints)]);
        rows.push(['Cells', n(meta.numCells)]);
        const blocks = (meta.cellBlocks ?? [])
            .map((b) => `${b.type} × ${n(b.numCells)}`)
            .join(', ');
        if (blocks) rows.push(['Cell types', blocks]);
        const names = [...(meta.pointDataNames ?? []), ...(meta.cellDataNames ?? [])];
        if (names.length) rows.push(['Data arrays', names.join(', ')]);
    }
    rows.push(['Rendered', `${n(renderedCells)} surface cells`]);

    $('info').replaceChildren(
        ...rows.flatMap(([key, value]) => {
            const dt = document.createElement('dt');
            dt.textContent = key;
            const dd = document.createElement('dd');
            dd.textContent = value;
            return [dt, dd];
        })
    );
    show($('info-section'));
}

/**
 * File cells (global, block-major indices) behind surface cells, or null when
 * that cannot be told. A solid's facets carry their owner cell; a surface mesh
 * of polygon blocks only is numbered like the file.
 */
function fileCellsOf(surfaceCells: ArrayLike<number>): number[] | null {
    const data = renderer.data;
    if (!data) return null;
    const parent = data.getCellData().getArrayByName('surface:parent_cell');
    if (parent) {
        const values = parent.getData();
        return [...new Set(Array.from(surfaceCells, (i) => Number(values[i])))].sort((a, b) => a - b);
    }
    if (surfaceOnlyPolygons) return Array.from(surfaceCells);
    return null;
}

/** What the Refine form can be limited to: the file's cell regions, saved groups and the selection. */
function refineTargets(): RefineTarget[] {
    if (opsPanel.pipeline.some((op) => !opsPreserveCellIds([op]))) return [];
    const out: RefineTarget[] = [];
    const sel = selection.idsReport();
    if (sel && sel.surfaceCells.length > 0) {
        const cells = fileCellsOf(sel.surfaceCells);
        if (cells && cells.length) {
            out.push({
                value: 'selection',
                label: `Selection (${cells.length.toLocaleString('en-US')} cells)`,
                spec: { cells },
            });
        }
    }
    for (const g of groups.listForRefine()) {
        if (g.user) {
            const cells = fileCellsOf(g.surfaceCells);
            if (cells && cells.length) out.push({ value: `group:${g.name}`, label: `Group ${g.name}`, spec: { cells } });
        } else {
            out.push({ value: `region:${g.name}`, label: `Region ${g.name}`, spec: { region: g.name } });
        }
    }
    return out;
}

/** The saved-groups checkbox of the Convert section: shown only when it can apply. */
function syncConvertGroups(): void {
    const usable = !__VIEWER_EMBEDDED__ && groups.userGroups().length > 0 && opsPanel.pipeline.length === 0;
    show($('convert-groups-row'), usable);
    show($('convert-groups-hint'), usable);
}

/** (Re)build the colour-by menu: the arrays, then any groups, keeping the choice if it survives. */
function fillColorBy(selected?: string): void {
    setOptions(
        $<HTMLSelectElement>('color-by'),
        [
            { value: SOLID_COLOR, label: 'Solid colour' },
            ...arrays.map((a) => ({ value: colorKey(a), label: a.label })),
            ...groups.colorOptions(),
        ],
        selected
    );
}

function applyColorBy(key: string): void {
    if (groups.owns(key)) {
        groups.applyColoring(key);
        return;
    }
    groups.clearColoring();
    const range: ScalarRange | null = renderer.setColorBy(key);
    const entry = arrays.find((a) => colorKey(a) === key);
    if (!range || !entry) {
        legend.hide();
        return;
    }
    legend.show(entry.label, renderer.colormap, range);
}

// --- static wiring -------------------------------------------------------- //

$('color-by').addEventListener('change', (e) =>
    applyColorBy((e.target as HTMLSelectElement).value)
);
$('colormap').addEventListener('change', (e) => {
    const name = (e.target as HTMLSelectElement).value;
    renderer.setColormap(name);
    legend.setColormap(name);
});
$('representation').addEventListener('change', (e) =>
    renderer.setRepresentation(
        Number((e.target as HTMLSelectElement).value) as RepresentationMode
    )
);
$('opacity').addEventListener('input', (e) =>
    renderer.setOpacity(Number((e.target as HTMLInputElement).value) / 100)
);
$('edges').addEventListener('change', (e) =>
    renderer.setEdgeVisibility((e.target as HTMLInputElement).checked)
);
decorate('h-mesh', 'info');
decorate('h-color', 'colormap');
prependGlyph('tool-reset', 'maximize');
prependGlyph('tool-screenshot', 'camera');
prependGlyph('tool-inspect', 'crosshair');
prependGlyph('tool-select', 'boxSelect');
decorate('ops-undo', 'undo');
decorate('ops-revert', 'revert');
decorate('op-quality-apply', 'quality');
decorate('op-clean-apply', 'orphan');
decorate('op-smooth-apply', 'smooth');
decorate('op-refine-apply', 'refine');
decorate('op-partition-apply', 'partition');
decorate('op-section-apply', 'cut');
decorate('convert-download', 'export');

$('tool-reset').addEventListener('click', () => renderer.resetCamera());
$('tool-screenshot').addEventListener('click', () => {
    const base = (window.__viewerState.name || 'mesh').replace(/\.[^.]*$/, '');
    void downloadPng(renderer, `${base}.png`);
});

const canvas = renderer.view.getCanvas();
if (canvas) setupOrientationCube(renderer, canvas, () => view.refresh());

// --- click to inspect ----------------------------------------------------- //

const picker = new Picker(renderer);
const visibility = new Visibility((mask) => renderer.setVisibleMask(mask));
/** Set once the worker module is loaded; stays null on the offline page. */
let regionsFetcher: ((epoch: number, cells: number) => Promise<RegionsResult>) | null = null;

const selection: SelectionController = new SelectionController({
    renderer,
    picker,
    visibility,
    palette: () => SCENE_THEMES[theme.resolved].palette,
    onToolOn: () => {
        setInspecting(false);
        view.setMeasuring(false);
    },
    onSaveGroup: (name, mask, color) => groups.addUserGroup(name, mask, color),
    onShowAll: () => groups.showAll(),
    onState: (state) => {
        window.__viewerState.selection = state;
    },
});
const groups: GroupsController = new GroupsController({
    renderer,
    visibility,
    legend,
    palette: () => SCENE_THEMES[theme.resolved].palette,
    selectMask: (mask, how) => selection.selectMask(mask, how),
    fetch: () => regionsFetcher,
    onListChanged: () => {
        fillColorBy();
        syncConvertGroups();
    },
    onState: (state) => {
        window.__viewerState.groups = state;
    },
});
visibility.onChange(() => selection.refresh());
let inspecting = false;

const analysis: AnalysisController = initAnalysis({
    renderer,
    select: (mask) => selection.selectMask(mask),
    download: (name, text, mime) => {
        const url = URL.createObjectURL(new Blob([text], { type: mime }));
        const a = document.createElement('a');
        a.href = url;
        a.download = name;
        a.click();
        URL.revokeObjectURL(url);
    },
    onState: (state) => {
        window.__viewerState.analysis = state;
    },
    fileBase: () => (window.__viewerState.name || 'mesh').replace(/\.[^.]*$/, ''),
});

const view: ViewController = initView({
    renderer,
    legend,
    picker,
    onMeasureOn: () => {
        setInspecting(false);
        selection.setTool('off');
    },
    notify: (message) => setStatus('ready', message),
    onState: (state) => {
        window.__viewerState.view = state;
    },
});

function formatComponents(values: number[]): string {
    return values
        .map((v) => (Number.isFinite(v) ? String(Number(v.toPrecision(6))) : String(v)))
        .join(', ');
}

function renderPick(info: PickInfo | null): void {
    window.__viewerState.picked = info;
    const popover = $('pick-popover');
    if (!info) {
        popover.hidden = true;
        return;
    }

    const body = $('pick-body');
    const rows: HTMLElement[] = [];
    const add = (key: string, value: string) => {
        const dt = document.createElement('dt');
        dt.textContent = key;
        const dd = document.createElement('dd');
        dd.textContent = value;
        rows.push(dt, dd);
    };
    const heading = (text: string) => {
        const el = document.createElement('div');
        el.className = 'pick-heading';
        el.textContent = text;
        rows.push(el);
    };
    const note = (text: string) => {
        const el = document.createElement('div');
        el.className = 'pick-note';
        el.textContent = text;
        rows.push(el);
    };

    heading('Cell');
    add('Surface cell', String(info.cellId));
    add('Type', info.cellType);
    if (info.parentCell !== null) {
        // After an operation the id names a cell of the operated mesh, not of
        // the file -- say which rather than let it be read as the original.
        const applied = opsPanel.pipeline.length > 0;
        add(applied ? 'Cell (after ops)' : 'Original cell', String(info.parentCell));
    }
    for (const v of info.cellValues) add(v.name, formatComponents(v.components));
    if (info.cellValues.length === 0) note('no cell data');

    heading('Nearest point');
    add('Point', String(info.pointId));
    add('Position', formatComponents(info.position));
    for (const v of info.pointValues) add(v.name, formatComponents(v.components));
    if (info.pointValues.length === 0) note('no point data');

    body.replaceChildren(...rows);
    popover.hidden = false;
}

function setInspecting(on: boolean): void {
    inspecting = on;
    const button = $('tool-inspect');
    button.setAttribute('aria-pressed', String(on));
    if (!on) {
        picker.clear();
        renderPick(null);
    }
}

$('tool-inspect').addEventListener('click', () => {
    if (!inspecting && !picker.canPick()) {
        setStatus('ready', 'this mesh is too large to inspect');
        return;
    }
    setInspecting(!inspecting);
    // The inspector, measuring and the selection tool all own a plain click.
    if (inspecting) {
        selection.setTool('off');
        view.setMeasuring(false);
    }
});
$('pick-close').addEventListener('click', () => setInspecting(false));

if (canvas) {
    // The listener goes on the *container*, not the canvas: vtk.js's interactor
    // takes pointer capture on the container, so `click` is dispatched there
    // and a canvas-bound listener never fires at all. Coordinates are still
    // measured against the canvas, which is what the picker works in.
    //
    // Click, never move: a CPU ray/cell intersection per mouse-move would make
    // a large mesh feel broken, and hover inspection is not worth that.
    $('render').addEventListener('click', (event: MouseEvent) => {
        if (!inspecting) return;
        const rect = canvas.getBoundingClientRect();
        renderPick(
            picker.pickAt(event.clientX - rect.left, event.clientY - rect.top, rect.height)
        );
    });
}

// --- chrome: menus, shortcuts, sections ------------------------------------ //

const sidebar = initSidebar(() => renderer.resize());
const themeSelect = $<HTMLSelectElement>('theme-select');
themeSelect.value = theme.mode;
themeSelect.addEventListener('change', () => {
    if (isThemeMode(themeSelect.value)) theme.setMode(themeSelect.value);
});

/** Change a control the way a user would, so every listener on it runs. */
function driveControl(el: HTMLInputElement | HTMLSelectElement, value: string | boolean): void {
    if (el instanceof HTMLInputElement && typeof value === 'boolean') el.checked = value;
    else el.value = String(value);
    el.dispatchEvent(new Event('change', { bubbles: true }));
}

function syncViewMenu(): void {
    setMenuChecked('view.edges', $<HTMLInputElement>('edges').checked);
    setMenuChecked('view.wireframe', $<HTMLSelectElement>('representation').value === '1');
    setMenuChecked('view.inspect', inspecting);
}
$('edges').addEventListener('change', syncViewMenu);
$('representation').addEventListener('change', syncViewMenu);
$('tool-inspect').addEventListener('click', () => queueMicrotask(syncViewMenu));
$('pick-close').addEventListener('click', () => queueMicrotask(syncViewMenu));

const WEB_ONLY = () => !__VIEWER_EMBEDDED__;

registerAction({
    id: 'file.open',
    label: 'Open a mesh',
    shortcut: 'Ctrl+O',
    group: 'File',
    enabled: WEB_ONLY,
    run: () => $<HTMLInputElement>('file').click(),
});
registerAction({
    id: 'file.screenshot',
    label: 'Save a PNG screenshot',
    shortcut: 'P',
    group: 'File',
    run: () => $('tool-screenshot').click(),
});
registerAction({
    id: 'file.convert',
    label: 'Convert to another format',
    group: 'File',
    enabled: WEB_ONLY,
    run: () => {
        const section = $('convert-section');
        if (section.hidden) return;
        setCollapsed(section, false);
        $('convert-format').focus();
    },
});
registerAction({
    id: 'view.fit',
    label: 'Fit the mesh to the window',
    shortcut: 'F',
    group: 'View',
    run: () => renderer.resetCamera(),
});
registerAction({
    id: 'view.edges',
    label: 'Show or hide edges',
    shortcut: 'E',
    group: 'View',
    run: () => driveControl($<HTMLInputElement>('edges'), !$<HTMLInputElement>('edges').checked),
});
registerAction({
    id: 'view.wireframe',
    label: 'Toggle wireframe',
    shortcut: 'W',
    group: 'View',
    run: () =>
        driveControl(
            $<HTMLSelectElement>('representation'),
            $<HTMLSelectElement>('representation').value === '1' ? '2' : '1'
        ),
});
registerAction({
    id: 'view.inspect',
    label: 'Inspect cells by clicking',
    shortcut: 'Shift+I',
    group: 'View',
    run: () => $('tool-inspect').click(),
});
registerAction({
    id: 'view.sidebar',
    label: 'Show or hide the sidebar',
    shortcut: 'Ctrl+B',
    group: 'View',
    run: () => sidebar.toggle(),
});
registerAction({
    id: 'help.shortcuts',
    label: 'Show the keyboard shortcuts',
    shortcut: '?',
    group: 'Help',
    run: showShortcuts,
});
registerAction({
    id: 'help.docs',
    label: 'Open the documentation',
    group: 'Help',
    run: () => window.open('https://loumalouomega.github.io/meshioplusplus/viewer', '_blank', 'noopener'),
});

initPrefs();

document.addEventListener('keydown', (e) => {
    const target = e.target as HTMLElement | null;
    // Typing in a field, a select or a dialog is never a shortcut.
    if (target?.closest('input, select, textarea, dialog, [contenteditable="true"]')) return;
    const action = findByShortcut(e);
    if (!action) return;
    e.preventDefault();
    action.run();
});

initMenus();
initSections();
syncViewMenu();

// --- payload resolution --------------------------------------------------- //

/** Prefix a control with one of the line glyphs, if the control exists. */
function prependGlyph(id: string, name: Parameters<typeof glyph>[0]): void {
    const el = maybe(id);
    if (!el) return;
    el.insertAdjacentHTML('afterbegin', glyph(name));
}

/** Prefix a control with its inline SVG icon, if the control exists. */
function decorate(id: string, icon: keyof typeof ICONS): void {
    const el = maybe(id);
    if (!el) return;
    const span = document.createElement('span');
    span.className = 'toolbar-icon';
    span.innerHTML = ICONS[icon];
    el.prepend(span);
}

function decodeBase64(text: string): ArrayBuffer {
    const binary = atob(text);
    const bytes = new Uint8Array(binary.length);
    for (let i = 0; i < binary.length; i++) bytes[i] = binary.charCodeAt(i);
    return bytes.buffer;
}

/** Stats Python precomputed, which an offline page cannot derive itself. */
function embeddedStats(): Record<string, unknown> | null {
    const slot = maybe('stats-payload');
    const text = slot?.textContent?.trim();
    if (!text) return null;
    try {
        return JSON.parse(text) as Record<string, unknown>;
    } catch {
        return null;
    }
}

function embeddedConfig(): { colorBy?: string } | null {
    const slot = maybe('view-config');
    const text = slot?.textContent?.trim();
    if (!text) return null;
    try {
        return JSON.parse(text) as { colorBy?: string };
    } catch {
        return null;
    }
}

/** Whatever the page was given, however it arrived. */
function present(
    result: RenderResult,
    name: string,
    keepSelection = false,
    idsPreserved = true
): void {
    show3d(result.surface, result.meta, name, keepSelection, idsPreserved);
    opsPanel.showReport(result.report);
}

let applyOpsImpl: (ops: OpSpec[]) => void = () => {};
function applyOps(ops: OpSpec[]): void {
    applyOpsImpl(ops);
}

async function boot(): Promise<void> {
    const embedded = $('vtp-payload').textContent?.trim();
    if (embedded) {
        setStatus('loading');
        show3d(decodeBase64(embedded), null, document.title);

        const stats = embeddedStats();
        if (stats) {
            const dl = $('info');
            for (const [key, value] of Object.entries(stats)) {
                const dt = document.createElement('dt');
                dt.textContent = key;
                const dd = document.createElement('dd');
                dd.textContent =
                    typeof value === 'number' ? String(Number(value.toPrecision(6))) : String(value);
                dl.append(dt, dd);
            }
        }
        // Python precomputed whatever the panel can offer offline.
        opsPanel.setAvailable(true);

        const config = embeddedConfig();
        if (config?.colorBy) {
            $<HTMLSelectElement>('color-by').value = config.colorBy;
            applyColorBy(config.colorBy);
        }
        return;
    }

    if (__VIEWER_EMBEDDED__) {
        // The embedded build has no other way to get a mesh, and an empty slot
        // means the page was written out wrong -- say so rather than sit blank.
        setStatus('error', 'this page carries no mesh payload');
        return;
    }

    const { applyOps: apply, convertTo, getRegions, init, openFile } = await import('./worker/client');
    regionsFetcher = (epoch, cells) => getRegions(epoch, cells);

    applyOpsImpl = (ops) => {
        opsPanel.setBusy(true);
        setStatus('loading', ops.length ? 'applying operations…' : 'reverting…');
        apply(ops)
            .then((result) => {
                const preserved = opsPreserveCellIds(ops);
                present(result, currentFile?.name ?? window.__viewerState.name, preserved, preserved);
                opsPanel.setBusy(false);
            })
            .catch((e: Error) => {
                opsPanel.setBusy(false);
                setStatus('error', e.message);
            });
    };

    async function open(file: File): Promise<void> {
        currentFile = file;
        setStatus('loading', `reading ${file.name}…`);
        try {
            const result = await openFile(file, (stage) =>
                setStatus(
                    'loading',
                    stage === 'read' ? 'reading…' : 'extracting surface…'
                )
            );
            // Prefer the file's own bounding box, but fall back to the
            // rendered extent: VTU's fast metadata path reports no bbox, and
            // without one the section control would silently do nothing.
            present(result, file.name);
            opsPanel.reset(
                result.bounds ?? renderer.bounds(),
                result.meta.numCells,
                // Only point arrays: a cell field is piecewise constant and has
                // no level set, so isosurface would reject it by name.
                result.meta.pointDataNames ?? [],
                // Cell arrays drive refine's predicate, which wants exactly the
                // piecewise-constant fields isosurface cannot use.
                result.meta.cellDataNames ?? []
            );
            opsPanel.setAvailable(true);
            show($('convert-section'));
        } catch (e) {
            setStatus('error', (e as Error).message);
        }
    }

    wireFileInput(open);
    wireDropZone(open);

    $('convert-download').addEventListener('click', () => {
        const outFormat = $<HTMLSelectElement>('convert-format').value;
        setStatus('loading', `converting to ${outFormat}…`);
        const wanted =
            !$('convert-groups-row').hidden && $<HTMLInputElement>('convert-groups').checked
                ? groups.userGroups()
                : [];
        convertTo(outFormat, opsPanel.pipeline, wanted.length ? wanted : undefined)
            .then(({ bytes, groups: written }) => {
                const base = (currentFile?.name ?? 'mesh').replace(/\.[^.]*$/, '');
                download(bytes, `${base}.${outFormat}`);
                if (written && written.requested > 0) {
                    const kept = written.kept.length;
                    setStatus(
                        'ready',
                        kept === written.requested
                            ? `wrote ${kept} group${kept === 1 ? '' : 's'}: ${written.kept.join(', ')}`
                            : `${outFormat} kept ${kept} of ${written.requested} groups` +
                                  (kept ? `: ${written.kept.join(', ')}` : '')
                    );
                } else {
                    setStatus('ready');
                }
            })
            .catch((e: Error) => setStatus('error', e.message));
    });

    setStatus('loading', 'starting meshio++…');
    try {
        const { formats } = await init();
        const writers = formats.writers ?? [];
        if (writers.length) {
            setOptions(
                $<HTMLSelectElement>('convert-format'),
                writers.map((f) => ({ value: f, label: f })),
                writers.includes('vtu') ? 'vtu' : writers[0]
            );
        }
        setStatus('idle');
    } catch (e) {
        setStatus('error', `could not start meshio++: ${(e as Error).message}`);
        return;
    }

    const url = new URLSearchParams(location.search).get('url');
    if (url) {
        try {
            const response = await fetch(url);
            if (!response.ok) throw new Error(`${response.status} ${response.statusText}`);
            const blob = await response.blob();
            await open(new File([blob], url.split('/').pop() || 'mesh'));
        } catch (e) {
            setStatus('error', `could not fetch ${url}: ${(e as Error).message}`);
        }
    }

    wireSamples(open);
}

function wireFileInput(open: (file: File) => void): void {
    const input = $<HTMLInputElement>('file');
    input.addEventListener('change', () => {
        const file = input.files?.[0];
        // Reset so re-picking the same file fires `change` again.
        input.value = '';
        if (file) void open(file);
    });
}

function wireDropZone(open: (file: File) => void): void {
    const overlay = $('drop-overlay');
    let depth = 0;
    window.addEventListener('dragenter', (e) => {
        e.preventDefault();
        // dragenter/dragleave fire per child element, so a plain toggle
        // flickers; count depth instead.
        if (++depth === 1) overlay.hidden = false;
    });
    window.addEventListener('dragover', (e) => e.preventDefault());
    window.addEventListener('dragleave', (e) => {
        e.preventDefault();
        if (--depth <= 0) {
            depth = 0;
            overlay.hidden = true;
        }
    });
    window.addEventListener('drop', (e) => {
        e.preventDefault();
        depth = 0;
        overlay.hidden = true;
        const file = e.dataTransfer?.files?.[0];
        if (file) void open(file);
    });
}

/** Bundled example meshes: tiny, generated, and exempted from Git LFS. */
const SAMPLES = [
    { file: 'cube.vtp', label: 'Cube (surface)' },
    { file: 'block.vtu', label: 'Block (volume)' },
    { file: 'wave.vtp', label: 'Wave (point data)' },
    { file: 'zones.msh', label: 'Zones (regions)' },
];

function wireSamples(open: (file: File) => Promise<void>): void {
    $('samples').replaceChildren(
        ...SAMPLES.map((sample) => {
            const button = document.createElement('button');
            button.type = 'button';
            button.className = 'sample';
            button.textContent = sample.label;
            button.addEventListener('click', async () => {
                setStatus('loading', `loading ${sample.label}…`);
                const response = await fetch(
                    `${import.meta.env.BASE_URL}samples/${sample.file}`
                );
                const blob = await response.blob();
                await open(new File([blob], sample.file));
            });
            return button;
        })
    );
}

function download(bytes: ArrayBuffer, filename: string): void {
    const url = URL.createObjectURL(new Blob([bytes]));
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    a.click();
    URL.revokeObjectURL(url);
}

void boot();
