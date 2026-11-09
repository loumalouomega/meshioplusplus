/**
 * The measuring mode (roadmap 7.9.1): click a feature of a scan to grow and
 * fit a primitive on it, read its values and sigma, the angle and distance
 * between two fitted features, the scan's deviation from one of them with the
 * share within a tolerance, and its wall thickness.
 *
 * A click goes to the worker's `fitFeature` request (`fitPrimitives` in seed
 * mode on the surface as drawn); the deviation view (a seeded `measure`
 * element, then `deviation` from it) and the thickness view are ordinary
 * pipeline operations appended to the Ops panel, so they show as
 * chips, undo like any other, and colour the surface by their field.
 */
import { registerAction } from '../app/actions';
import { type FittedElement, elementFromRow, relate } from '../analysis/elements';
import { FeatureOverlay } from '../render/feature';
import type { Picker } from '../render/pick';
import type { Renderer } from '../render/renderer';
import type { MeasuringState, Vector3 } from '../types';
import { $, show } from '../ui/dom';
import { glyph } from '../ui/glyphs';
import type { FeatureResult, OpReport, OpSpec } from '../worker/protocol';

export interface MeasuringOptions {
    renderer: Renderer;
    picker: Picker;
    /** The worker's fit, or null on a page without one (the embedded viewer). */
    fit: () => ((seed: Vector3, types: string[], tolerance: number) => Promise<FeatureResult>) | null;
    /** Append operations to the Ops panel's pipeline (which re-renders). */
    append: (ops: OpSpec[]) => void;
    /** Colour the surface by an array, as the colour-by menu would. */
    colorBy: (location: 'point' | 'cell', name: string) => void;
    /** The mode took the plain click: turn the other click tools off. */
    onModeOn: () => void;
    notify: (message: string) => void;
    onState: (state: MeasuringState) => void;
}

export interface MeasuringController {
    /** A new surface was drawn (a file or a pipeline change). */
    attachSurface(): void;
    /** A new file: forget the features. */
    reset(): void;
    /** The step report of the last render, for the deviation and thickness rows. */
    showReport(report: OpReport): void;
    setActive(on: boolean): void;
    readonly active: boolean;
}

interface Fitted {
    element: FittedElement;
    triangles: Float32Array;
    seed: Vector3;
    types: string[];
    tolerance: number;
}

function fmt(v: number | null | undefined, digits = 5): string {
    if (v === null || v === undefined || !Number.isFinite(v)) return '–';
    return String(Number(v.toPrecision(digits)));
}

export function initMeasuring(opts: MeasuringOptions): MeasuringController {
    const { renderer, picker } = opts;
    const overlay = new FeatureOverlay(renderer);
    const list = $('feat-list');
    const relationBox = $('feat-relation');
    const resultBox = $('feat-view-result');
    const button = $('tool-fit');
    button.insertAdjacentHTML('afterbegin', glyph('shapes'));

    let active = false;
    let busy = false;
    let error: string | null = null;
    let fitted: Fitted[] = [];
    /** Indices into `fitted`, oldest first, at most two. */
    let selected: number[] = [];
    let counter = 0;
    let view: MeasuringState['view'] = null;
    /** What the next render's report is expected to carry. */
    let pendingView: 'deviation' | 'thickness' | null = null;

    function relation(): MeasuringState['relation'] {
        if (selected.length !== 2) return null;
        const [a, b] = selected.map((i) => fitted[i].element);
        const r = relate(a, b);
        return { a: a.name, b: b.name, angleDeg: r.angleDeg, distance: r.distance, what: r.what };
    }

    function publish(): void {
        opts.onState({
            active,
            busy,
            error,
            elements: fitted.map(({ element: e }) => ({
                name: e.name,
                kind: e.kind,
                radius: e.radius,
                angleDeg: e.angleDeg,
                sigma: e.sigma,
                formError: e.formError,
                numCells: e.numCells,
            })),
            selected: selected.map((i) => fitted[i].element.name),
            relation: relation(),
            view,
        });
    }

    function describe(e: FittedElement): string {
        const rows = [`<div><span>${e.name}</span><b>${e.kind}</b></div>`];
        if (e.kind === 'sphere' || e.kind === 'cylinder')
            rows.push(`<div><span>Radius</span><b>${fmt(e.radius)}</b></div>`,
                `<div><span>Diameter</span><b>${fmt(2 * e.radius)}</b></div>`);
        if (e.kind === 'torus')
            rows.push(`<div><span>Radii</span><b>${fmt(e.radius)}, ${fmt(e.radius2)}</b></div>`);
        if (e.kind === 'cone') rows.push(`<div><span>Half-angle</span><b>${fmt(e.angleDeg)}°</b></div>`);
        const dirLabel = e.kind === 'plane' ? 'Normal' : 'Axis';
        if (e.kind !== 'sphere') rows.push(`<div><span>${dirLabel}</span><b>${e.direction.map((v) => fmt(v, 4)).join(', ')}</b></div>`);
        rows.push(
            `<div><span>${e.kind === 'plane' ? 'Point' : 'Centre'}</span><b>${e.point.map((v) => fmt(v, 5)).join(', ')}</b></div>`,
            `<div><span>σ</span><b>${fmt(e.sigma, 3)}</b></div>`,
            `<div><span>Form error</span><b>${fmt(e.formError, 3)}</b></div>`,
            `<div><span>Cells / used / rejected</span><b>${e.numCells} / ${e.numUsed} / ${e.numRejected}</b></div>`
        );
        return rows.join('');
    }

    function render(): void {
        list.replaceChildren(
            ...fitted.map((f, i) => {
                const item = document.createElement('button');
                item.type = 'button';
                item.className = 'feat-item sel-stats';
                item.setAttribute('aria-pressed', String(selected.includes(i)));
                item.dataset.name = f.element.name;
                item.innerHTML = describe(f.element);
                item.addEventListener('click', () => toggle(i));
                return item;
            })
        );
        const r = relation();
        if (!r) {
            relationBox.innerHTML =
                fitted.length > 1 ? '<p class="hint">Select two features to relate them.</p>' : '';
        } else {
            const rows = [`<div><span>${r.a} – ${r.b}</span><b>${r.what}</b></div>`];
            rows.push(`<div><span>Angle</span><b>${r.angleDeg === null ? '–' : `${fmt(r.angleDeg, 4)}°`}</b></div>`);
            rows.push(`<div><span>Distance</span><b>${r.distance === null ? 'they meet' : fmt(r.distance)}</b></div>`);
            relationBox.innerHTML = rows.join('');
        }
        const status = error ? `<p class="hint feat-error">${error}</p>` : busy ? '<p class="hint">fitting…</p>' : '';
        $('feat-status').innerHTML = status;
        overlay.set(
            fitted.map((f) => f.triangles),
            new Set(selected)
        );
        if (view?.kind === 'deviation') {
            resultBox.innerHTML = [
                `<div><span>Deviation from</span><b>${view.feature}</b></div>`,
                `<div><span>Within ± ${fmt(view.tolerance)}</span><b>${fmt(100 * (view.shareWithin ?? NaN), 4)} %</b></div>`,
                `<div><span>RMS</span><b>${fmt(view.rms, 3)}</b></div>`,
                `<div><span>Min / max</span><b>${fmt(view.min, 3)} / ${fmt(view.max, 3)}</b></div>`,
            ].join('');
        } else if (view?.kind === 'thickness') {
            resultBox.innerHTML = [
                `<div><span>Thickness</span><b>${view.method}</b></div>`,
                `<div><span>Min / mean</span><b>${fmt(view.min, 4)} / ${fmt(view.mean, 4)}</b></div>`,
                `<div><span>Missed</span><b>${view.missed ?? '–'}</b></div>`,
            ].join('');
        } else {
            resultBox.innerHTML = '';
        }
        $<HTMLButtonElement>('feat-deviation').disabled = selected.length === 0 || busy;
        publish();
    }

    function toggle(i: number): void {
        selected = selected.includes(i) ? selected.filter((j) => j !== i) : [...selected, i].slice(-2);
        render();
    }

    function setActive(on: boolean): void {
        active = on;
        button.setAttribute('aria-pressed', String(on));
        if (on) {
            opts.onModeOn();
            show($('features-section'));
        }
        publish();
    }

    async function fitAt(seed: Vector3): Promise<void> {
        const fit = opts.fit();
        if (!fit || busy) return;
        const kind = $<HTMLSelectElement>('feat-type').value;
        const types = kind === 'auto' ? ['plane', 'sphere', 'cylinder', 'cone'] : [kind];
        const tolerance = Number($<HTMLInputElement>('feat-tol').value) || 0.002;
        busy = true;
        error = null;
        render();
        try {
            const r = await fit(seed, types, tolerance);
            const element = elementFromRow(`F${counter + 1}`, r.patchType, r.params, {
                rms: r.rms,
                sigma: r.sigma,
                formError: r.formError,
                numCells: r.numCells,
                numUsed: r.numUsed,
                numRejected: r.numRejected,
            });
            counter += 1;
            fitted.push({ element, triangles: r.triangles, seed, types, tolerance });
            selected = [...selected, fitted.length - 1].slice(-2);
        } catch (e) {
            error = (e as Error).message;
        } finally {
            busy = false;
            render();
        }
    }

    button.addEventListener('click', () => setActive(!active));
    registerAction({
        id: 'measuring.fit',
        label: 'Fit a feature where you click',
        shortcut: 'G',
        group: 'Tools',
        run: () => button.click(),
    });

    const canvas = renderer.view.getCanvas();
    $('render').addEventListener('click', (event: MouseEvent) => {
        if (!active || !canvas) return;
        const rect = canvas.getBoundingClientRect();
        const info = picker.pickAt(event.clientX - rect.left, event.clientY - rect.top, rect.height);
        if (!info) return;
        picker.clear();
        void fitAt(info.hit);
    });

    $('feat-clear').addEventListener('click', () => {
        fitted = [];
        selected = [];
        counter = 0;
        error = null;
        render();
    });

    $('feat-deviation').addEventListener('click', () => {
        const last = selected[selected.length - 1];
        if (last === undefined) return;
        const f = fitted[last];
        const tolerance = Number($<HTMLInputElement>('feat-dev-tol').value) || 0;
        pendingView = 'deviation';
        view = { kind: 'deviation', feature: f.element.name, tolerance, method: null };
        // The mesh's own element table names the feature: a seeded element
        // fitted by `measure`, then the deviation from it.
        const kind = f.types.length === 1 ? f.types[0] : 'auto';
        const spec = { elements: [{ name: f.element.name, fit: kind, seeds: [f.seed], tolerance: f.tolerance }] };
        opts.append([
            { op: 'measure', spec: JSON.stringify(spec) },
            { op: 'deviation', element: f.element.name, tolerance, bins: 20 },
        ]);
    });

    $('feat-thickness').addEventListener('click', () => {
        const method = $<HTMLSelectElement>('feat-thk-method').value as 'ray' | 'sphere' | 'cone';
        pendingView = 'thickness';
        view = { kind: 'thickness', feature: null, tolerance: null, method };
        opts.append([{ op: 'thickness', method }]);
    });

    return {
        attachSurface(): void {
            show($('features-section'));
            render();
        },
        reset(): void {
            fitted = [];
            selected = [];
            counter = 0;
            error = null;
            view = null;
            pendingView = null;
            render();
        },
        showReport(report: OpReport): void {
            if (!pendingView || !view) return;
            const step = [...report.steps].reverse().find((s) => s.op === pendingView);
            if (!step) return;
            const num = (key: string): number | null => {
                const v = (step as unknown as Record<string, unknown>)[key];
                return typeof v === 'number' ? v : null;
            };
            if (pendingView === 'deviation') {
                view = { ...view, shareWithin: num('shareWithin'), rms: num('rms'), min: num('min'), max: num('max') };
                opts.colorBy('point', 'deviation:distance');
            } else {
                view = { ...view, min: num('min'), mean: num('mean'), missed: num('numMissed') };
                opts.colorBy('cell', 'thickness:distance');
            }
            pendingView = null;
            render();
        },
        setActive,
        get active() {
            return active;
        },
    };
}
