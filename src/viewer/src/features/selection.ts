/**
 * The selection tool: sets of surface cells you build by clicking, dragging a
 * box or a lasso, thresholding a field, or growing/inverting what you have.
 *
 * Pure logic lives in `selection/*` (masks, projection, thresholds, tested in
 * Node); this file is the glue: pointer gestures on the canvas, the sidebar
 * panel, the highlight overlay and cell hiding.
 *
 * A selection is a mask over the *surface* cells. Ids reported to the user
 * are surface cell ids, with the originating file cell (`surface:parent_cell`)
 * alongside when the surface came from a solid.
 */
import { registerAction } from '../app/actions';
import type { Visibility } from '../app/visibility';
import { arrayOf } from '../render/scalars';
import type { Picker } from '../render/pick';
import type { Renderer } from '../render/renderer';
import {
    type Combine,
    type Mask,
    combine,
    connectedMask,
    countMask,
    emptyMask,
    growMask,
    idsOf,
    invertMask,
    maskStats,
    shrinkMask,
    thresholdCellValues,
    thresholdPointValues,
} from '../selection/mask';
import { SelectionOverlay } from '../selection/overlay';
import { type Shape, cellsInShape, facingMask, projectPoints } from '../selection/project';
import type { ArrayEntry, SelectionState, Vector3 } from '../types';
import { $, show } from '../ui/dom';
import { cssColor } from '../ui/theme';
import { setStatusBar } from './statusbar';

export type SelectionTool = 'off' | 'click' | 'box' | 'lasso';

interface SelectionSet {
    id: number;
    name: string;
    color: Vector3;
    mask: Mask;
}

export interface SelectionOptions {
    renderer: Renderer;
    picker: Picker;
    /** Cell visibility shared with the groups panel; this tool owns the `selection` layer. */
    visibility: Visibility;
    /** Promote a set to a named group. */
    onSaveGroup: (name: string, mask: Mask, color: Vector3) => void;
    /** "Show all" also reveals what the groups panel hid. */
    onShowAll: () => void;
    /** Colours for new sets (the scene palette of the current scheme). */
    palette: () => Vector3[];
    /** The selection tool took over the pointer: turn the inspector off. */
    onToolOn: () => void;
    /** Publish state for the test hook. */
    onState: (state: SelectionState) => void;
}

/** A drag shorter than this many pixels is a click. */
const DRAG_PIXELS = 4;

export class SelectionController {
    private readonly o: SelectionOptions;
    private readonly overlay: SelectionOverlay;
    private sets: SelectionSet[] = [];
    private activeId: number | null = null;
    private nextId = 1;
    private tool: SelectionTool = 'off';
    private lastTool: Exclude<SelectionTool, 'off'> = 'box';
    private arrays: ArrayEntry[] = [];

    constructor(options: SelectionOptions) {
        this.o = options;
        this.overlay = new SelectionOverlay(options.renderer);
        this.wirePanel();
        this.wireGestures();
        this.registerActions();
    }

    // --- surface lifecycle ----------------------------------------------------

    /**
     * A surface was (re)loaded. `keep` keeps the sets when the cells are the
     * same ones in the same order (an operation that preserves cell ids);
     * otherwise everything starts over, because a mask over different cells
     * would highlight the wrong ones.
     */
    attachSurface(keep: boolean): void {
        const n = this.cellCount();
        const same = keep && this.sets.length > 0 && this.sets.every((s) => s.mask.length === n);
        this.overlay.clear();
        if (!same) {
            this.sets = [];
            this.activeId = null;
            this.o.visibility.set('selection', null);
            if (n > 0) this.newSet();
        }
        this.arrays = this.o.renderer.arrays().filter((a) => this.valuesOf(a) !== null);
        this.fillArrayList();
        show($('selection-section'), n > 0);
        this.redraw();
    }

    /** Redraw highlights after the visible cells changed from outside. */
    refresh(): void {
        this.redraw();
    }

    /** Fold a mask (for example a group's cells) into the active set. */
    selectMask(mask: Mask, how: Combine = 'replace'): void {
        this.apply(mask, how);
    }

    private cellCount(): number {
        return this.o.renderer.surfaceGeometry()?.topo.numCells ?? 0;
    }

    // --- tool state ---------------------------------------------------------------

    get currentTool(): SelectionTool {
        return this.tool;
    }

    setTool(tool: SelectionTool): void {
        if (tool === this.tool) return;
        this.tool = tool;
        if (tool !== 'off') {
            this.lastTool = tool;
            this.o.onToolOn();
        }
        for (const b of document.querySelectorAll<HTMLElement>('#sel-tools [data-tool]')) {
            b.setAttribute('aria-pressed', String(b.dataset.tool === tool));
        }
        $('tool-select').setAttribute('aria-pressed', String(tool !== 'off'));
        $('stage').classList.toggle('selecting', tool !== 'off');
        this.cancelGesture();
        this.publish();
    }

    toggleTool(): void {
        this.setTool(this.tool === 'off' ? this.lastTool : 'off');
    }

    // --- sets -------------------------------------------------------------------------

    private active(): SelectionSet | null {
        return this.sets.find((s) => s.id === this.activeId) ?? null;
    }

    private newSet(): SelectionSet {
        const palette = this.o.palette();
        const set: SelectionSet = {
            id: this.nextId,
            name: `Selection ${this.nextId}`,
            color: palette[(this.nextId - 1) % palette.length] ?? [1, 0.8, 0],
            mask: emptyMask(this.cellCount()),
        };
        this.nextId++;
        this.sets.push(set);
        this.activeId = set.id;
        return set;
    }

    private setActiveMask(mask: Mask): void {
        const set = this.active() ?? this.newSet();
        set.mask = mask;
        this.redraw();
    }

    /** Fold a gesture's cells into the active set the way the modifiers ask. */
    private apply(incoming: Mask, how: Combine): void {
        const set = this.active() ?? this.newSet();
        this.setActiveMask(combine(set.mask, incoming, how));
    }

    private redraw(): void {
        for (const set of this.sets) this.overlay.set(set.id, set.color, set.mask, this.o.visibility.combined());
        this.renderSets();
        this.renderStats();
        const total = this.active() ? countMask(this.active()!.mask) : 0;
        setStatusBar({ selected: total });
        this.publish();
    }

    private publish(): void {
        this.o.onState({
            tool: this.tool,
            active: this.activeId,
            sets: this.sets.map((s) => ({ id: s.id, name: s.name, count: countMask(s.mask) })),
            hidden: this.o.visibility.hiddenCount(),
        });
    }

    // --- gestures -----------------------------------------------------------------------

    private gesture: { kind: 'box' | 'lasso'; points: number[]; id: number } | null = null;

    private wireGestures(): void {
        const root = $('render');
        let down: { x: number; y: number } | null = null;

        // Capture phase, before vtk.js's interactor sees the press: box and
        // lasso own a plain left drag so it does not also rotate the camera.
        // Holding Alt hands the drag back to the camera.
        root.addEventListener(
            'pointerdown',
            (e) => {
                if (this.tool === 'off' || e.button !== 0) return;
                down = { x: e.clientX, y: e.clientY };
                if (this.tool === 'click' || e.altKey) return;
                e.stopImmediatePropagation();
                e.preventDefault();
                const p = this.local(e);
                this.gesture = { kind: this.tool, points: [p.x, p.y], id: e.pointerId };
                root.setPointerCapture(e.pointerId);
                this.drawGesture(p);
            },
            true
        );
        root.addEventListener('pointermove', (e) => {
            if (!this.gesture || e.pointerId !== this.gesture.id) return;
            this.drawGesture(this.local(e));
        });
        root.addEventListener('pointerup', (e) => {
            const g = this.gesture;
            if (!g || e.pointerId !== g.id) return;
            const p = this.local(e);
            this.gesture = null;
            this.hideGesture();
            const dragged = down
                ? Math.hypot(e.clientX - down.x, e.clientY - down.y) >= DRAG_PIXELS
                : true;
            if (!dragged) return;
            const how: Combine = e.shiftKey ? 'add' : e.ctrlKey || e.metaKey ? 'subtract' : 'replace';
            this.commitGesture(g.kind, g.points, p, how);
        });
        root.addEventListener('pointercancel', () => this.cancelGesture());

        // A click (not a drag that ended on the canvas) picks one cell.
        root.addEventListener('click', (e: MouseEvent) => {
            if (this.tool !== 'click') return;
            if (down && Math.hypot(e.clientX - down.x, e.clientY - down.y) >= DRAG_PIXELS) return;
            const canvas = this.o.renderer.view.getCanvas();
            if (!canvas) return;
            const rect = canvas.getBoundingClientRect();
            const id = this.o.picker.cellAt(e.clientX - rect.left, e.clientY - rect.top, rect.height);
            const n = this.cellCount();
            const picked = emptyMask(n);
            if (id !== null) picked[id] = 1;
            // Plain click replaces, Shift adds, Ctrl toggles; clicking empty
            // space with no modifier clears.
            const how: Combine = e.shiftKey ? 'add' : e.ctrlKey || e.metaKey ? 'toggle' : 'replace';
            this.apply(picked, how);
        });

        document.addEventListener('keydown', (e) => {
            if (e.key === 'Escape' && this.gesture) {
                e.preventDefault();
                this.cancelGesture();
            }
        });
    }

    private local(e: PointerEvent): { x: number; y: number } {
        const canvas = this.o.renderer.view.getCanvas();
        const rect = (canvas ?? $('render')).getBoundingClientRect();
        return { x: e.clientX - rect.left, y: e.clientY - rect.top };
    }

    private drawGesture(p: { x: number; y: number }): void {
        const g = this.gesture;
        if (!g) return;
        const overlay = $('sel-overlay');
        // `hidden` is not an IDL property of SVG elements, so use the attribute.
        overlay.removeAttribute('hidden');
        const rect = $('sel-rect');
        const poly = $('sel-poly');
        if (g.kind === 'box') {
            const [x0, y0] = g.points;
            rect.setAttribute('x', String(Math.min(x0, p.x)));
            rect.setAttribute('y', String(Math.min(y0, p.y)));
            rect.setAttribute('width', String(Math.abs(p.x - x0)));
            rect.setAttribute('height', String(Math.abs(p.y - y0)));
            rect.style.display = '';
            poly.style.display = 'none';
        } else {
            const last = g.points.length;
            if (Math.hypot(p.x - g.points[last - 2], p.y - g.points[last - 1]) >= 3) {
                g.points.push(p.x, p.y);
            }
            poly.setAttribute('points', g.points.join(' '));
            poly.style.display = '';
            rect.style.display = 'none';
        }
    }

    private hideGesture(): void {
        $('sel-overlay').setAttribute('hidden', '');
    }

    private cancelGesture(): void {
        if (this.gesture) {
            try {
                $('render').releasePointerCapture(this.gesture.id);
            } catch {
                // already released
            }
        }
        this.gesture = null;
        this.hideGesture();
    }

    /** The cells a finished box or lasso covers. */
    private commitGesture(
        kind: 'box' | 'lasso',
        points: number[],
        end: { x: number; y: number },
        how: Combine
    ): void {
        const cam = this.o.renderer.cameraProjection();
        const geometry = this.o.renderer.surfaceGeometry();
        if (!cam || !geometry) return;
        const shape: Shape =
            kind === 'box'
                ? { kind: 'rect', rect: { x0: points[0], y0: points[1], x1: end.x, y1: end.y } }
                : { kind: 'polygon', points };
        const proj = projectPoints(geometry.centroids, cam.matrix, cam.width, cam.height);

        let eligible: Mask | null = this.o.visibility.combined();
        if ($<HTMLInputElement>('sel-front').checked) {
            const facing = facingMask(
                geometry.centroids,
                geometry.normals,
                cam.parallel ? null : cam.eye,
                cam.parallel ? cam.direction : null
            );
            if (eligible) {
                const both = new Uint8Array(facing.length);
                for (let i = 0; i < both.length; i++) both[i] = facing[i] & eligible[i];
                eligible = both;
            } else {
                eligible = facing;
            }
        }
        const incoming = cellsInShape(proj, shape, eligible ?? undefined);
        // A plain drag that selects nothing still replaces, i.e. clears.
        this.apply(incoming, how);
    }

    // --- operations on the active set --------------------------------------------------

    private modify(fn: (mask: Mask) => Mask): void {
        const set = this.active();
        if (!set) return;
        this.setActiveMask(fn(set.mask));
    }

    invert(): void {
        this.modify((m) => {
            const out = invertMask(m);
            const visible = this.o.visibility.combined();
            if (visible) for (let i = 0; i < out.length; i++) out[i] &= visible[i];
            return out;
        });
    }

    grow(): void {
        const g = this.o.renderer.surfaceGeometry();
        if (g) this.modify((m) => growMask(m, g.topo, g.pointCells));
    }

    shrink(): void {
        const g = this.o.renderer.surfaceGeometry();
        if (g) this.modify((m) => shrinkMask(m, g.topo, g.pointCells));
    }

    connected(): void {
        const g = this.o.renderer.surfaceGeometry();
        if (g) this.modify((m) => connectedMask(m, g.topo, g.pointCells));
    }

    clear(): void {
        this.modify((m) => emptyMask(m.length));
    }

    hideSelected(): void {
        const set = this.active();
        if (!set || countMask(set.mask) === 0) return;
        const n = this.cellCount();
        const current = this.o.visibility.get('selection');
        const next = current ? current.slice() : new Uint8Array(n).fill(1);
        for (let i = 0; i < n; i++) if (set.mask[i]) next[i] = 0;
        this.o.visibility.set('selection', next);
    }

    isolate(): void {
        const set = this.active();
        if (!set || countMask(set.mask) === 0) return;
        this.o.visibility.set('selection', set.mask.slice());
    }

    showAll(): void {
        this.o.visibility.set('selection', null);
        this.o.onShowAll();
    }

    /** Bring the active set into view. */
    frame(): void {
        const set = this.active();
        if (set && countMask(set.mask) > 0) this.o.renderer.frameCells(set.mask);
    }

    saveAsGroup(): void {
        const set = this.active();
        if (set && countMask(set.mask) > 0) this.o.onSaveGroup(set.name, set.mask.slice(), set.color);
    }

    // --- selection by value ---------------------------------------------------------------

    private valuesOf(entry: ArrayEntry): ReturnType<typeof arrayOf> {
        const a = arrayOf(this.o.renderer.data, entry);
        return a && a.getNumberOfComponents() >= 1 ? a : null;
    }

    private fillArrayList(): void {
        const select = $<HTMLSelectElement>('sel-array');
        select.replaceChildren(
            ...this.arrays.map(
                (a, i) =>
                    new Option(`${a.label} (${a.location === 'cell' ? 'cells' : 'points'})`, String(i))
            )
        );
        $('sel-value').hidden = this.arrays.length === 0;
        this.fillRange();
    }

    private fillRange(): void {
        const entry = this.arrays[Number($<HTMLSelectElement>('sel-array').value)];
        if (!entry) return;
        const array = this.valuesOf(entry);
        if (!array) return;
        const range = this.o.renderer.rangeOf(entry);
        $<HTMLInputElement>('sel-min').value = range ? String(Number(range.min.toPrecision(6))) : '';
        $<HTMLInputElement>('sel-max').value = range ? String(Number(range.max.toPrecision(6))) : '';
        $<HTMLSelectElement>('sel-rule').disabled = entry.location === 'cell';
    }

    private selectByValue(): void {
        const entry = this.arrays[Number($<HTMLSelectElement>('sel-array').value)];
        const geometry = this.o.renderer.surfaceGeometry();
        const array = entry ? this.valuesOf(entry) : null;
        if (!entry || !array || !geometry) return;
        const lo = Number($<HTMLInputElement>('sel-min').value);
        const hi = Number($<HTMLInputElement>('sel-max').value);
        if (!Number.isFinite(lo) || !Number.isFinite(hi)) return;
        const values = array.getData();
        const comps = array.getNumberOfComponents();
        let incoming: Mask;
        if (entry.location === 'cell') {
            incoming = thresholdCellValues(values, comps, entry.component, lo, hi, geometry.topo.numCells);
        } else {
            const rule = $<HTMLSelectElement>('sel-rule').value === 'any' ? 'any' : 'all';
            incoming = thresholdPointValues(values, comps, entry.component, lo, hi, geometry.topo, rule);
        }
        const visible = this.o.visibility.combined();
        if (visible) for (let i = 0; i < incoming.length; i++) incoming[i] &= visible[i];
        this.apply(incoming, $<HTMLSelectElement>('sel-how').value as Combine);
    }

    // --- ids -------------------------------------------------------------------------------------

    /** Surface ids of the active set and, when known, their file cells. */
    idsReport(): { name: string; surfaceCells: number[]; parentCells: number[] | null } | null {
        const set = this.active();
        if (!set) return null;
        const ids = idsOf(set.mask);
        const parent = this.o.renderer.data?.getCellData().getArrayByName('surface:parent_cell');
        let parents: number[] | null = null;
        if (parent) {
            const data = parent.getData();
            parents = [...new Set(Array.from(ids, (i) => Number(data[i])))].sort((a, b) => a - b);
        }
        return { name: set.name, surfaceCells: Array.from(ids), parentCells: parents };
    }

    private download(text: string, name: string, type: string): void {
        const url = URL.createObjectURL(new Blob([text], { type }));
        const a = document.createElement('a');
        a.href = url;
        a.download = name;
        a.click();
        URL.revokeObjectURL(url);
    }

    private exportIds(format: 'csv' | 'json'): void {
        const report = this.idsReport();
        if (!report || report.surfaceCells.length === 0) return;
        const base = report.name.replace(/[^\w.-]+/g, '_');
        if (format === 'json') {
            this.download(JSON.stringify(report, null, 2), `${base}.json`, 'application/json');
            return;
        }
        const data = this.o.renderer.data?.getCellData().getArrayByName('surface:parent_cell');
        const rows = report.surfaceCells.map((id) => (data ? `${id},${Number(data.getData()[id])}` : String(id)));
        const header = data ? 'surface_cell,parent_cell' : 'surface_cell';
        this.download([header, ...rows].join('\n') + '\n', `${base}.csv`, 'text/csv');
    }

    // --- panel ------------------------------------------------------------------------------------

    private wirePanel(): void {
        for (const b of document.querySelectorAll<HTMLButtonElement>('#sel-tools [data-tool]')) {
            b.addEventListener('click', () => {
                const t = b.dataset.tool as SelectionTool;
                this.setTool(this.tool === t ? 'off' : t);
            });
        }
        $('tool-select').addEventListener('click', () => this.toggleTool());
        $('sel-new').addEventListener('click', () => {
            this.newSet();
            this.redraw();
        });
        $('sel-clear').addEventListener('click', () => this.clear());
        $('sel-invert').addEventListener('click', () => this.invert());
        $('sel-grow').addEventListener('click', () => this.grow());
        $('sel-shrink').addEventListener('click', () => this.shrink());
        $('sel-connected').addEventListener('click', () => this.connected());
        $('sel-hide').addEventListener('click', () => this.hideSelected());
        $('sel-isolate').addEventListener('click', () => this.isolate());
        $('sel-show').addEventListener('click', () => this.showAll());
        $('sel-frame').addEventListener('click', () => this.frame());
        $('sel-group').addEventListener('click', () => this.saveAsGroup());
        $('sel-array').addEventListener('change', () => this.fillRange());
        $('sel-by-value').addEventListener('click', () => this.selectByValue());
        $('sel-csv').addEventListener('click', () => this.exportIds('csv'));
        $('sel-json').addEventListener('click', () => this.exportIds('json'));
        $('sel-copy').addEventListener('click', () => {
            const report = this.idsReport();
            if (report) void navigator.clipboard?.writeText(report.surfaceCells.join(',')).catch(() => {});
        });
    }

    private renderSets(): void {
        const list = $('sel-sets');
        list.replaceChildren(
            ...this.sets.map((set) => {
                const li = document.createElement('li');
                li.dataset.id = String(set.id);
                if (set.id === this.activeId) li.setAttribute('aria-current', 'true');
                const swatch = document.createElement('span');
                swatch.className = 'sel-swatch';
                swatch.style.background = cssColor(set.color);
                const name = document.createElement('span');
                name.className = 'sel-name';
                name.textContent = set.name;
                name.title = 'Double-click to rename';
                const count = document.createElement('span');
                count.className = 'sel-count';
                count.textContent = countMask(set.mask).toLocaleString('en-US');
                const remove = document.createElement('button');
                remove.type = 'button';
                remove.className = 'sel-remove';
                remove.title = 'Delete this set';
                remove.textContent = '×';
                li.append(swatch, name, count, remove);

                li.addEventListener('click', () => {
                    if (this.activeId === set.id) return;
                    this.activeId = set.id;
                    this.redraw();
                });
                name.addEventListener('dblclick', (ev) => {
                    ev.stopPropagation();
                    const input = document.createElement('input');
                    input.type = 'text';
                    input.className = 'sel-rename';
                    input.value = set.name;
                    const done = (keep: boolean) => {
                        if (keep && input.value.trim()) set.name = input.value.trim();
                        this.renderSets();
                        this.publish();
                    };
                    input.addEventListener('keydown', (k) => {
                        if (k.key === 'Enter') done(true);
                        else if (k.key === 'Escape') done(false);
                        k.stopPropagation();
                    });
                    input.addEventListener('blur', () => done(true));
                    name.replaceWith(input);
                    input.focus();
                    input.select();
                });
                remove.addEventListener('click', (ev) => {
                    ev.stopPropagation();
                    this.overlay.remove(set.id);
                    this.sets = this.sets.filter((s) => s.id !== set.id);
                    if (this.activeId === set.id) this.activeId = this.sets[this.sets.length - 1]?.id ?? null;
                    this.redraw();
                });
                return li;
            })
        );
    }

    private renderStats(): void {
        const box = $('sel-stats');
        const set = this.active();
        const data = this.o.renderer.data;
        if (!set || !data || countMask(set.mask) === 0) {
            box.replaceChildren();
            return;
        }
        const table = document.createElement('table');
        const head = document.createElement('tr');
        for (const t of ['', 'min', 'mean', 'max']) {
            const th = document.createElement('th');
            th.textContent = t;
            head.append(th);
        }
        const rows: HTMLElement[] = [];
        for (const entry of this.arrays.filter((a) => a.location === 'cell')) {
            const array = this.valuesOf(entry);
            if (!array) continue;
            const s = maskStats(array.getData(), array.getNumberOfComponents(), entry.component, set.mask);
            if (s.count === 0) continue;
            const row = document.createElement('tr');
            const label = document.createElement('th');
            label.textContent = entry.label;
            row.append(label);
            for (const v of [s.min, s.mean, s.max]) {
                const td = document.createElement('td');
                td.textContent = String(Number(v.toPrecision(5)));
                row.append(td);
            }
            rows.push(row);
        }
        if (rows.length === 0) {
            box.replaceChildren();
            return;
        }
        table.append(head, ...rows);
        box.replaceChildren(table);
    }

    // --- actions ---------------------------------------------------------------------------------------

    private registerActions(): void {
        const needs = () => this.cellCount() > 0;
        const add = (id: string, label: string, shortcut: string | undefined, run: () => void) =>
            registerAction({ id, label, shortcut, group: 'Selection', enabled: needs, run });
        add('select.toggle', 'Turn the selection tool on or off', 'S', () => this.toggleTool());
        add('select.box', 'Box selection', 'B', () => this.setTool(this.tool === 'box' ? 'off' : 'box'));
        add('select.lasso', 'Lasso selection', 'L', () => this.setTool(this.tool === 'lasso' ? 'off' : 'lasso'));
        add('select.invert', 'Invert the selection', 'I', () => this.invert());
        add('select.grow', 'Grow the selection by one ring', 'G', () => this.grow());
        add('select.shrink', 'Shrink the selection by one ring', 'Shift+G', () => this.shrink());
        add('select.hide', 'Hide the selected cells', 'H', () => this.hideSelected());
        add('select.isolate', 'Hide everything but the selection', 'X', () => this.isolate());
        add('select.showAll', 'Show every cell again', 'Shift+H', () => this.showAll());
        add('select.clear', 'Clear the active selection', 'Shift+C', () => this.clear());
    }
}
