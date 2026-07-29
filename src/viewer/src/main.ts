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
import { COLORMAPS, DEFAULT_COLORMAP } from './render/colormaps';
import { Legend } from './render/legend';
import { setupOrientationCube } from './render/orientation';
import { Renderer, SOLID_COLOR, colorKey } from './render/renderer';
import type { RepresentationMode } from './render/renderer';
import { Picker } from './render/pick';
import { downloadPng } from './render/screenshot';
import { OpsPanel } from './ops/panel';
import { ICONS } from './ui/icons';
import { $, maybe, setOptions, show } from './ui/dom';
import type { ArrayEntry, PickInfo, ScalarRange, ViewerStatus } from './types';
import type { MeshMeta, OpSpec, RenderResult } from './worker/protocol';

const renderer = new Renderer($('render'));
const opsPanel = new OpsPanel({ onChange: (ops) => applyOps(ops) });
let arrays: ArrayEntry[] = [];
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

function show3d(vtp: ArrayBuffer, meta: MeshMeta | null, name: string): void {
    const info = renderer.load(vtp);
    arrays = info.arrays;

    Object.assign(window.__viewerState, {
        numPoints: info.numPoints,
        numCells: info.numCells,
        arrays: arrays.map((a) => a.label),
        name,
    });

    setOptions(
        $<HTMLSelectElement>('color-by'),
        [
            { value: SOLID_COLOR, label: 'Solid colour' },
            ...arrays.map((a) => ({ value: colorKey(a), label: a.label })),
        ],
        SOLID_COLOR
    );
    show($('color-section'));
    legend.hide();
    picker.clear();
    renderPick(null);
    renderInfo(meta, info.numCells, name);
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

function applyColorBy(key: string): void {
    const range: ScalarRange | null = renderer.setColorBy(key);
    const entry = arrays.find((a) => colorKey(a) === key);
    if (!range || !entry) {
        legend.hide();
        return;
    }
    legend.show(entry.label, renderer.colormap, range);
}

// --- static wiring -------------------------------------------------------- //

setOptions(
    $<HTMLSelectElement>('colormap'),
    COLORMAPS.map((name) => ({ value: name, label: name })),
    DEFAULT_COLORMAP
);

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
decorate('tool-reset', 'reset');
decorate('tool-screenshot', 'screenshot');
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
if (canvas) setupOrientationCube(renderer, canvas);

// --- click to inspect ----------------------------------------------------- //

const picker = new Picker(renderer);
let inspecting = false;

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

// --- payload resolution --------------------------------------------------- //

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
function present(result: RenderResult, name: string): void {
    show3d(result.vtp, result.meta, name);
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

    const { applyOps: apply, convertTo, init, openFile } = await import('./worker/client');

    applyOpsImpl = (ops) => {
        opsPanel.setBusy(true);
        setStatus('loading', ops.length ? 'applying operations…' : 'reverting…');
        apply(ops)
            .then((result) => {
                present(result, currentFile?.name ?? window.__viewerState.name);
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
        convertTo(outFormat, opsPanel.pipeline)
            .then(({ bytes }) => {
                const base = (currentFile?.name ?? 'mesh').replace(/\.[^.]*$/, '');
                download(bytes, `${base}.${outFormat}`);
                setStatus('ready');
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
