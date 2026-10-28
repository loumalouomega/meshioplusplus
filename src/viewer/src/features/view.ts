/**
 * Viewing tools: colour scale, standard views, projection, camera bookmarks,
 * lighting, deformed shape, clip plane and the measure tool.
 *
 * Ideas and defaults come from the Kratos MDPA Preview extension's nav dock,
 * lighting and bookmarks panels (same author). Everything here is pure client
 * code, so it ships in the offline page too.
 */
import { registerAction } from '../app/actions';
import { VIEW_KEYS, VIEW_NAMES, type ViewName, parseBookmarks, parseCameraJson, withBookmark } from '../render/camera';
import type { Bookmark } from '../render/camera';
import { colormapGroups } from '../render/colormaps';
import type { Legend } from '../render/legend';
import { MeasureOverlay, measureBetween } from '../render/measure';
import type { Picker } from '../render/pick';
import type { Renderer } from '../render/renderer';
import { bandCount } from '../render/scale';
import type { Vector3, ViewState } from '../types';
import { $, show } from '../ui/dom';
import { glyph } from '../ui/glyphs';

const BOOKMARK_KEY = 'meshioplusplus.viewer.bookmarks';

export const DEFAULT_LIGHTING = { ambient: 10, diffuse: 100, specular: 0 } as const;

export interface ViewOptions {
    renderer: Renderer;
    legend: Legend;
    picker: Picker;
    /** The measure tool took over the plain click: turn the other click tools off. */
    onMeasureOn: () => void;
    /** Tell the user something short. */
    notify: (message: string) => void;
    onState: (state: ViewState) => void;
}

export interface ViewController {
    /** A new surface was loaded: refresh the vector list, drop the deformation, keep the rest. */
    attachSurface(): void;
    /** The camera changed from outside (the orientation cube): republish the test-hook state. */
    refresh(): void;
    setMeasuring(on: boolean): void;
    readonly measuring: boolean;
}

function loadBookmarks(): Bookmark[] {
    try {
        const raw = localStorage.getItem(BOOKMARK_KEY);
        return raw ? parseBookmarks(JSON.parse(raw) as unknown) : [];
    } catch {
        return [];
    }
}

function saveBookmarks(list: Bookmark[]): void {
    try {
        localStorage.setItem(BOOKMARK_KEY, JSON.stringify(list));
    } catch {
        // best effort
    }
}

const fmt = (v: number) => String(Number(v.toPrecision(6)));

export function initView(opts: ViewOptions): ViewController {
    const { renderer, legend, picker } = opts;

    let deformed = 0;
    let measuring = false;
    let marks: Vector3[] = [];
    let bookmarks = loadBookmarks();
    let scaleBands = 0;

    function publish(): void {
        const cam = renderer.getCameraState();
        const d = cam.position.map((v, i) => v - cam.focalPoint[i]);
        const len = Math.hypot(d[0], d[1], d[2]) || 1;
        const round = (v: number) => Math.round((v / len) * 1000) / 1000 + 0;
        opts.onState({
            parallel: cam.parallel,
            direction: [round(d[0]), round(d[1]), round(d[2])],
            logScale: renderer.logActive,
            bands: scaleBands,
            clip: $<HTMLInputElement>('clip-on').checked,
            deformed,
            measuring,
            measured: marks.length === 2 ? measureBetween(marks[0], marks[1]).distance : null,
            bookmarks: bookmarks.map((b) => b.name),
        });
    }

    // --- colour map list, scale and bands ---------------------------------- //

    const colormap = $<HTMLSelectElement>('colormap');
    function fillColormaps(): void {
        const safe = $<HTMLInputElement>('colormap-safe').checked;
        const keep = renderer.colormap;
        const groups = colormapGroups(safe);
        const present = groups.some((g) => g.names.includes(keep));
        const nodes: HTMLElement[] = groups.map((g) => {
            const group = document.createElement('optgroup');
            group.label = g.category[0].toUpperCase() + g.category.slice(1);
            group.append(...g.names.map((n) => new Option(n, n, false, n === keep)));
            return group;
        });
        if (!present) nodes.unshift(new Option(keep, keep, false, true));
        colormap.replaceChildren(...nodes);
        colormap.value = keep;
    }
    $('colormap-safe').addEventListener('change', fillColormaps);
    fillColormaps();

    function applyScale(): void {
        const wantLog = $<HTMLSelectElement>('color-scale').value === 'log';
        const bands = bandCount($<HTMLSelectElement>('color-bands').value);
        scaleBands = bands;
        renderer.setScalarOptions({ log: wantLog, bands });
        legend.setScale(renderer.logActive, bands);
        if (wantLog && !renderer.logActive) opts.notify('Log scale needs a strictly positive range.');
        publish();
    }
    $('color-scale').addEventListener('change', applyScale);
    $('color-bands').addEventListener('change', applyScale);
    // A new range may make the log scale usable (or not): re-evaluate.
    $('color-by').addEventListener('change', () => queueMicrotask(applyScale));

    // --- standard views and projection ------------------------------------ //

    const projection = $<HTMLSelectElement>('view-projection');
    function setView(name: ViewName): void {
        renderer.resetCamera(name);
        publish();
    }
    for (const button of document.querySelectorAll<HTMLButtonElement>('#view-buttons [data-view]')) {
        const name = button.dataset.view as ViewName;
        if (VIEW_NAMES.includes(name)) button.addEventListener('click', () => setView(name));
    }
    projection.addEventListener('change', () => {
        renderer.setParallel(projection.value === 'ortho');
        publish();
    });
    for (const [key, name] of Object.entries(VIEW_KEYS)) {
        registerAction({
            id: `view.camera.${name}`,
            label: `Look from the ${name === 'iso' ? 'isometric corner' : name}`,
            shortcut: key,
            group: 'Camera',
            run: () => setView(name),
        });
    }
    registerAction({
        id: 'view.projection',
        label: 'Switch perspective / orthographic',
        shortcut: 'O',
        group: 'Camera',
        run: () => {
            projection.value = projection.value === 'ortho' ? 'persp' : 'ortho';
            projection.dispatchEvent(new Event('change', { bubbles: true }));
        },
    });

    // --- bookmarks --------------------------------------------------------- //

    const bmNote = $('bm-note');
    function renderBookmarks(): void {
        $('bm-list').replaceChildren(
            ...bookmarks.map((b) => {
                const li = document.createElement('li');
                const go = document.createElement('button');
                go.type = 'button';
                go.className = 'bm-go';
                go.textContent = b.name;
                go.title = 'Go to this view';
                go.addEventListener('click', () => {
                    renderer.setCameraState(b.camera);
                    projection.value = b.camera.parallel ? 'ortho' : 'persp';
                    publish();
                });
                const del = document.createElement('button');
                del.type = 'button';
                del.className = 'icon-btn';
                del.title = `Delete ${b.name}`;
                del.setAttribute('aria-label', `Delete ${b.name}`);
                del.innerHTML = glyph('x');
                del.addEventListener('click', () => {
                    bookmarks = bookmarks.filter((x) => x.name !== b.name);
                    saveBookmarks(bookmarks);
                    renderBookmarks();
                    publish();
                });
                li.append(go, del);
                return li;
            })
        );
    }
    renderBookmarks();
    $('bm-save').addEventListener('click', () => {
        const input = $<HTMLInputElement>('bm-name');
        const name = input.value.trim() || `View ${bookmarks.length + 1}`;
        bookmarks = withBookmark(bookmarks, { name, camera: renderer.getCameraState() });
        saveBookmarks(bookmarks);
        input.value = '';
        renderBookmarks();
        publish();
    });
    $('bm-copy').addEventListener('click', () => {
        const text = JSON.stringify(renderer.getCameraState());
        bmNote.textContent = 'Camera copied.';
        navigator.clipboard?.writeText(text).catch(() => {
            bmNote.textContent = 'Copy blocked by the browser: select it from the prompt instead.';
            window.prompt('Camera JSON', text);
        });
    });
    $('bm-paste').addEventListener('click', () => {
        const apply = (text: string) => {
            const state = parseCameraJson(text);
            if (!state) {
                bmNote.textContent = 'That is not a camera saved by this viewer.';
                return;
            }
            renderer.setCameraState(state);
            projection.value = state.parallel ? 'ortho' : 'persp';
            bmNote.textContent = 'Camera applied.';
            publish();
        };
        if (navigator.clipboard?.readText) {
            navigator.clipboard.readText().then(apply, () => apply(window.prompt('Paste the camera JSON') ?? ''));
        } else {
            apply(window.prompt('Paste the camera JSON') ?? '');
        }
    });

    // --- lighting ---------------------------------------------------------- //

    function applyLighting(): void {
        const v = (id: string) => Number($<HTMLInputElement>(id).value) / 100;
        const specular = v('light-specular');
        renderer.setLighting({
            ambient: v('light-ambient'),
            diffuse: v('light-diffuse'),
            specular,
            specularPower: 20,
        });
    }
    for (const id of ['light-ambient', 'light-diffuse', 'light-specular']) {
        $(id).addEventListener('input', applyLighting);
    }
    $('light-cull').addEventListener('change', (e) =>
        renderer.setBackfaceCulling((e.target as HTMLInputElement).checked)
    );
    $('light-reset').addEventListener('click', () => {
        $<HTMLInputElement>('light-ambient').value = String(DEFAULT_LIGHTING.ambient);
        $<HTMLInputElement>('light-diffuse').value = String(DEFAULT_LIGHTING.diffuse);
        $<HTMLInputElement>('light-specular').value = String(DEFAULT_LIGHTING.specular);
        $<HTMLInputElement>('light-cull').checked = false;
        applyLighting();
        renderer.setBackfaceCulling(false);
    });

    // --- deformed shape ---------------------------------------------------- //

    const warpArray = $<HTMLSelectElement>('warp-array');
    const warpNote = $('warp-note');
    function applyWarp(): void {
        const name = warpArray.value || null;
        const drawn = renderer.setDeform(name, Number($<HTMLInputElement>('warp-scale').value));
        deformed = name ? drawn : 0;
        warpNote.textContent = name
            ? `Largest displacement drawn: ${fmt(drawn)}. 100% moves the largest by 10% of the mesh size.`
            : 'Pick a vector point array to move the mesh by it.';
        publish();
    }
    warpArray.addEventListener('change', applyWarp);
    $('warp-scale').addEventListener('input', applyWarp);

    // --- clip plane -------------------------------------------------------- //

    function applyClip(): void {
        if (!$<HTMLInputElement>('clip-on').checked) {
            renderer.setClip(null);
            publish();
            return;
        }
        const axis = Number(document.querySelector<HTMLInputElement>('input[name="clip-axis"]:checked')?.value ?? 2);
        renderer.setClip({
            axis: axis === 0 || axis === 1 ? axis : 2,
            t: Number($<HTMLInputElement>('clip-pos').value) / 100,
            flip: $<HTMLInputElement>('clip-flip').checked,
        });
        publish();
    }
    $('clip-on').addEventListener('change', applyClip);
    $('clip-pos').addEventListener('input', applyClip);
    $('clip-flip').addEventListener('change', applyClip);
    for (const radio of document.querySelectorAll('input[name="clip-axis"]')) {
        radio.addEventListener('change', applyClip);
    }

    // --- measure ------------------------------------------------------------ //

    const overlay = new MeasureOverlay(renderer);
    const result = $('measure-result');

    function renderMeasure(): void {
        publish();
        if (marks.length === 0) {
            result.textContent = '';
            return;
        }
        const rows = marks.map((p, i) => `<div><span>Point ${i + 1}</span><b>${p.map(fmt).join(', ')}</b></div>`);
        if (marks.length === 2) {
            const m = measureBetween(marks[0], marks[1]);
            rows.push(
                `<div><span>Distance</span><b>${fmt(m.distance)}</b></div>`,
                `<div><span>Δ x, y, z</span><b>${m.delta.map(fmt).join(', ')}</b></div>`
            );
        }
        result.innerHTML = rows.join('');
    }

    function setMeasuring(on: boolean): void {
        measuring = on;
        $('tool-measure').setAttribute('aria-pressed', String(on));
        if (on) {
            opts.onMeasureOn();
            const details = $<HTMLDetailsElement>('view-measure');
            details.open = true;
            marks = [];
            overlay.clear();
            renderMeasure();
        } else {
            marks = [];
            overlay.clear();
            renderMeasure();
        }
    }
    $('tool-measure').addEventListener('click', () => setMeasuring(!measuring));
    registerAction({
        id: 'view.measure',
        label: 'Measure the distance between two points',
        shortcut: 'M',
        group: 'Tools',
        run: () => $('tool-measure').click(),
    });
    registerAction({
        id: 'view.measure.stop',
        label: 'Stop measuring',
        shortcut: 'Esc',
        group: 'Tools',
        enabled: () => measuring,
        run: () => setMeasuring(false),
    });

    const canvas = renderer.view.getCanvas();
    $('render').addEventListener('click', (event: MouseEvent) => {
        if (!measuring || !canvas) return;
        const rect = canvas.getBoundingClientRect();
        const info = picker.pickAt(event.clientX - rect.left, event.clientY - rect.top, rect.height);
        if (!info) return;
        picker.clear();
        if (marks.length >= 2) marks = [];
        marks.push($<HTMLInputElement>('measure-snap').checked ? info.position : info.hit);
        overlay.set(marks);
        renderMeasure();
    });

    // --- glyphs -------------------------------------------------------------- //

    const measureButton = $('tool-measure');
    measureButton.insertAdjacentHTML('afterbegin', glyph('ruler'));

    return {
        attachSurface(): void {
            const vectors = renderer.vectorArrays();
            const keep = warpArray.value;
            warpArray.replaceChildren(new Option('none', ''), ...vectors.map((n) => new Option(n, n)));
            if (vectors.includes(keep)) warpArray.value = keep;
            show($('view-section'));
            applyWarp();
            applyScale();
            applyClip();
            if (measuring) setMeasuring(true);
        },
        refresh: publish,
        setMeasuring,
        get measuring() {
            return measuring;
        },
    };
}
