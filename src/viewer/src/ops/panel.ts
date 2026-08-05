/**
 * The operations panel.
 *
 * Runs meshio++'s own mesh operations on the open file and re-renders the
 * result — quality, clean, smooth, refine, partition, a sectioning cut and
 * the isosurface of a point field.
 * This is what distinguishes the viewer from a generic vtk.js app: the whole
 * point of the demo is that a real mesh library is running in the browser.
 *
 * The panel owns the pipeline as a plain list. Apply pushes, Undo pops, Revert
 * clears, and every change re-sends the **whole list** to the worker, which
 * replays it against the original file bytes. So undo is exact rather than
 * approximate, and no operation needs an inverse.
 *
 * Under `__VIEWER_EMBEDDED__` (the wheel-bundled page) there is no WASM to run
 * anything, so the panel renders in a read-only "baked" mode instead: whatever
 * Python precomputed is offered as a colour-by, and the rest is disabled with
 * a reason.
 */
import { $, maybe, setOptions, show } from '../ui/dom';
import type { OpReport, OpSpec } from '../worker/protocol';
import { OP_DEFAULTS, describeOp } from '../worker/protocol';
import type { Vector3 } from '../types';

export interface OpsPanelCallbacks {
    /** The pipeline changed; re-render with this list. */
    onChange(ops: OpSpec[]): void;
}

/** Refining is exponential; beyond this, ask first. */
const CONFIRM_CELLS = 2_000_000;

export class OpsPanel {
    private readonly root: HTMLElement;
    private readonly chips: HTMLElement;
    private readonly warnings: HTMLElement;
    private readonly callbacks: OpsPanelCallbacks;

    private ops: OpSpec[] = [];
    private bounds: { min: Vector3; max: Vector3 } | null = null;
    private cellCount = 0;
    private busy = false;

    constructor(callbacks: OpsPanelCallbacks) {
        this.callbacks = callbacks;
        this.root = $('ops-section');
        this.chips = $('ops-chips');
        this.warnings = $('ops-warnings');

        if (__VIEWER_EMBEDDED__) {
            this.renderBakedNotice();
            return;
        }

        $('ops-undo').addEventListener('click', () => this.undo());
        $('ops-revert').addEventListener('click', () => this.revert());

        this.wireQuality();
        this.wireClean();
        this.wireSmooth();
        this.wireRefine();
        this.wirePartition();
        this.wireSection();
        this.wireIsosurface();
        this.wireGradient();
        this.wireVoxelize();
        this.wireComputeSdf();
        this.renderChips();
    }

    /** The current pipeline, for the convert-and-download path. */
    get pipeline(): OpSpec[] {
        return [...this.ops];
    }

    setAvailable(available: boolean): void {
        show(this.root, available);
    }

    /** A new mesh was opened; the pipeline no longer applies to it. */
    reset(
        bounds: { min: Vector3; max: Vector3 } | null,
        cellCount: number,
        pointArrays: string[] = [],
        cellArrays: string[] = []
    ): void {
        this.ops = [];
        this.bounds = bounds;
        this.cellCount = cellCount;
        this.warnings.replaceChildren();
        this.resetSectionSlider();
        this.setIsosurfaceArrays(pointArrays);
        this.setGradientArrays(pointArrays);
        this.setRefineArrays(cellArrays);
        this.renderChips();
    }

    setBusy(busy: boolean): void {
        this.busy = busy;
        this.root.classList.toggle('busy', busy);
        for (const el of this.root.querySelectorAll('button, input, select')) {
            (el as HTMLButtonElement).disabled = busy;
        }
    }

    /** Show the counters and caveats the C++ side reported. */
    showReport(report: OpReport): void {
        this.warnings.replaceChildren(
            ...report.warnings.map((text) => {
                const li = document.createElement('li');
                li.textContent = text;
                return li;
            })
        );
        this.renderChips(report);
    }

    // --- pipeline ---------------------------------------------------------

    private push(spec: OpSpec): void {
        if (this.busy) return;
        this.ops.push(spec);
        this.renderChips();
        this.callbacks.onChange(this.pipeline);
    }

    /** Replace the trailing op of the same kind, or append. Used by section. */
    private replaceOrPush(spec: OpSpec): void {
        if (this.busy) return;
        const last = this.ops[this.ops.length - 1];
        if (last && last.op === spec.op) this.ops[this.ops.length - 1] = spec;
        else this.ops.push(spec);
        this.renderChips();
        this.callbacks.onChange(this.pipeline);
    }

    private undo(): void {
        if (this.busy || this.ops.length === 0) return;
        this.ops.pop();
        this.renderChips();
        this.callbacks.onChange(this.pipeline);
    }

    private revert(): void {
        if (this.busy || this.ops.length === 0) return;
        this.ops = [];
        this.renderChips();
        this.callbacks.onChange(this.pipeline);
    }

    private removeAt(index: number): void {
        if (this.busy) return;
        this.ops.splice(index, 1);
        this.renderChips();
        this.callbacks.onChange(this.pipeline);
    }

    private renderChips(report?: OpReport): void {
        const undo = maybe<HTMLButtonElement>('ops-undo');
        const revert = maybe<HTMLButtonElement>('ops-revert');
        if (undo) undo.disabled = this.ops.length === 0;
        if (revert) revert.disabled = this.ops.length === 0;

        if (this.ops.length === 0) {
            const empty = document.createElement('p');
            empty.className = 'hint';
            empty.textContent = 'No operations applied.';
            this.chips.replaceChildren(empty);
            return;
        }

        this.chips.replaceChildren(
            ...this.ops.map((spec, i) => {
                const chip = document.createElement('span');
                chip.className = 'op-chip';
                chip.append(describeOp(spec));

                const counters = report?.steps[i];
                if (counters) {
                    const detail = Object.entries(counters)
                        .filter(([k, v]) => k !== 'op' && typeof v === 'number' && v !== 0)
                        .map(([k, v]) => `${k} ${v}`)
                        .join(', ');
                    if (detail) chip.title = detail;
                }

                const remove = document.createElement('button');
                remove.type = 'button';
                remove.className = 'op-chip-remove';
                remove.title = 'Remove this operation';
                remove.textContent = '×';
                remove.addEventListener('click', () => this.removeAt(i));
                chip.append(remove);
                return chip;
            })
        );
    }

    // --- per-operation wiring --------------------------------------------

    private wireQuality(): void {
        $('op-quality-apply').addEventListener('click', () =>
            this.push({ ...OP_DEFAULTS.quality })
        );
    }

    private wireClean(): void {
        $('op-clean-apply').addEventListener('click', () =>
            this.push({
                ...OP_DEFAULTS.clean,
                weld: $<HTMLInputElement>('op-clean-weld').checked,
                atol: Number($<HTMLInputElement>('op-clean-atol').value) || 1e-8,
                removeOrphans: $<HTMLInputElement>('op-clean-orphans').checked,
                dropDegenerate: $<HTMLInputElement>('op-clean-degenerate').checked,
                dropDuplicateCells: $<HTMLInputElement>('op-clean-duplicates').checked,
            })
        );
    }

    private wireSmooth(): void {
        $('op-smooth-apply').addEventListener('click', () =>
            this.push({
                ...OP_DEFAULTS.smooth,
                method: $<HTMLSelectElement>('op-smooth-method').value as
                    | 'laplacian'
                    | 'taubin',
                iterations: Number($<HTMLInputElement>('op-smooth-iterations').value) || 10,
                fixBoundary: $<HTMLInputElement>('op-smooth-boundary').checked,
            })
        );
    }

    /**
     * Offer the mesh's own cell arrays as refine predicates, plus the metrics a
     * preceding `quality` chip would produce: composing the two is the point of
     * the feature, and those names do not exist on the file itself.
     */
    private setRefineArrays(names: string[]): void {
        const select = maybe<HTMLSelectElement>('op-refine-array');
        if (!select) return;
        const offered = [...names];
        for (const metric of ['scaled_jacobian', 'aspect_ratio', 'skewness']) {
            const full = `quality:${metric}`;
            if (!offered.includes(full)) offered.push(full);
        }
        setOptions(
            select,
            [
                { value: '', label: 'every cell' },
                ...offered.map((name) => ({ value: name, label: name })),
            ],
            ''
        );
    }

    private wireRefine(): void {
        const levels = $<HTMLInputElement>('op-refine-levels');
        const note = $('op-refine-note');
        const array = maybe<HTMLSelectElement>('op-refine-array');
        const projected = () => {
            const n = Number(levels.value) || 1;
            // A predicate refines only the cells that match, and the closure
            // adds a bounded amount around them, so the uniform bound is only
            // honest when no predicate is set.
            if (array && array.value) return 0;
            // Every supported type splits 8-for-1 in 3D, 4-for-1 in 2D; 8 is
            // the honest upper bound to warn against.
            return this.cellCount * Math.pow(8, n);
        };
        const updateNote = () => {
            const cells = projected();
            if (!this.cellCount) {
                note.textContent = '';
            } else if (cells === 0) {
                // A predicate's growth depends on how many cells match and how
                // far the closure has to reach, neither of which is knowable
                // before running it. Say that rather than print a wrong number.
                note.textContent = 'only the matching cells, plus a conforming closure';
            } else {
                note.textContent = `≈ ${cells.toLocaleString('en-US')} cells`;
            }
        };
        levels.addEventListener('input', updateNote);
        array?.addEventListener('change', updateNote);
        updateNote();

        $('op-refine-apply').addEventListener('click', () => {
            const cells = projected();
            if (
                cells > CONFIRM_CELLS &&
                !window.confirm(
                    `Refining produces about ${cells.toLocaleString('en-US')} cells, ` +
                        'which may take a while and use a lot of memory. Continue?'
                )
            ) {
                return;
            }
            this.push({
                ...OP_DEFAULTS.refine,
                levels: Number(levels.value) || 1,
                array: array?.value ?? '',
                compare: ($<HTMLSelectElement>('op-refine-compare').value as '<' | '>') ?? '<',
                value: Number($<HTMLInputElement>('op-refine-value').value) || 0,
                closure: $<HTMLSelectElement>('op-refine-closure').value as
                    | 'redgreen'
                    | 'propagate'
                    | 'balanced',
            });
        });
    }

    private wirePartition(): void {
        $('op-partition-apply').addEventListener('click', () =>
            this.push({
                ...OP_DEFAULTS.partition,
                nparts: Number($<HTMLInputElement>('op-partition-nparts').value) || 2,
                method: $<HTMLSelectElement>('op-partition-method').value as
                    | 'auto'
                    | 'sfc'
                    | 'kahip',
            })
        );
    }

    private wireSection(): void {
        const slider = $<HTMLInputElement>('op-section-position');
        const apply = () => {
            if (!this.bounds) {
                // Silently doing nothing here looks exactly like a broken
                // button, which is how this was first found.
                this.warnings.replaceChildren(
                    Object.assign(document.createElement('li'), {
                        textContent: 'the mesh extent is unknown, so it cannot be sectioned',
                    })
                );
                return;
            }
            const axis = Number(
                (
                    this.root.querySelector<HTMLInputElement>(
                        'input[name="op-section-axis"]:checked'
                    ) ?? { value: '2' }
                ).value
            );
            const flip = $<HTMLInputElement>('op-section-flip').checked;
            const t = Number(slider.value) / 100;

            // World coordinates, not the slider fraction: the section must
            // stay put even if an earlier operation changes the mesh extent.
            const point: Vector3 = [0, 0, 0];
            point[axis] =
                this.bounds.min[axis] + t * (this.bounds.max[axis] - this.bounds.min[axis]);
            const normal: Vector3 = [0, 0, 0];
            normal[axis] = flip ? -1 : 1;

            this.replaceOrPush({ ...OP_DEFAULTS.section, point, normal, mode: 'all' });
        };
        $('op-section-apply').addEventListener('click', apply);
        $('op-section-flip').addEventListener('change', () => {
            // Only re-apply if a section is already active, so ticking the box
            // before choosing a position does not create one.
            if (this.ops.some((o) => o.op === 'section')) apply();
        });
    }

    private wireIsosurface(): void {
        $('op-iso-apply').addEventListener('click', () => {
            const array = $<HTMLSelectElement>('op-iso-array').value;
            if (!array) {
                // Only cell fields (or none) — say so rather than looking broken,
                // the way the section control does with an unknown extent.
                this.warnings.replaceChildren(
                    Object.assign(document.createElement('li'), {
                        textContent:
                            'this mesh carries no point_data, and a cell field is ' +
                            'piecewise constant — it has no level set',
                    })
                );
                return;
            }
            const component = $<HTMLInputElement>('op-iso-component').value;
            this.push({
                ...OP_DEFAULTS.isosurface,
                array,
                isovalue: Number($<HTMLInputElement>('op-iso-value').value) || 0,
                component: component === '' ? -1 : Number(component),
            });
        });
    }

    /** Offer the mesh's own point arrays; a cell field has no level set. */
    private setIsosurfaceArrays(names: string[]): void {
        const select = maybe<HTMLSelectElement>('op-iso-array');
        if (!select) return;
        setOptions(
            select,
            names.map((name) => ({ value: name, label: name })),
            names[0] ?? ''
        );
        const apply = maybe<HTMLButtonElement>('op-iso-apply');
        if (apply) apply.disabled = names.length === 0;
    }

    private wireGradient(): void {
        $('op-grad-apply').addEventListener('click', () => {
            const array = $<HTMLSelectElement>('op-grad-array').value;
            if (!array) {
                // The same reason isosurface refuses: a cell field is piecewise
                // constant, so it has no derivative either.
                this.warnings.replaceChildren(
                    Object.assign(document.createElement('li'), {
                        textContent:
                            'this mesh carries no point_data, and a cell field is ' +
                            'piecewise constant — it has no derivative',
                    })
                );
                return;
            }
            this.push({
                ...OP_DEFAULTS.gradient,
                array,
                operator: $<HTMLSelectElement>('op-grad-operator')
                    .value as 'gradient' | 'divergence' | 'curl',
                method: $<HTMLSelectElement>('op-grad-method').value as
                    | 'green-gauss'
                    | 'least-squares',
                location: $<HTMLSelectElement>('op-grad-location').value as 'point' | 'cell',
            });
        });
    }

    private wireVoxelize(): void {
        $('op-vox-apply').addEventListener('click', () => {
            const n = Math.max(1, Number($<HTMLInputElement>('op-vox-resolution').value) || 32);
            this.push({
                ...OP_DEFAULTS.voxelize,
                resolution: [n, n, n],
                fill: $<HTMLSelectElement>('op-vox-fill').value as
                    | 'all'
                    | 'surface'
                    | 'inside',
            });
        });
    }

    private wireComputeSdf(): void {
        const structure = maybe<HTMLSelectElement>('op-sdf-structure');
        const sync = () => {
            const octree = structure?.value === 'octree';
            // resolution/cellSize size a VOXEL grid; an octree's finest cell is
            // root/2^depth and is therefore already determined, so the two sets
            // of inputs are mutually exclusive rather than merely unused.
            const show = (id: string, on: boolean) => {
                const el = maybe<HTMLElement>(id);
                if (el?.parentElement) el.parentElement.hidden = !on;
            };
            show('op-sdf-resolution', !octree);
            show('op-sdf-root', octree);
            show('op-sdf-depth', octree);
        };
        structure?.addEventListener('change', sync);
        sync();

        $('op-sdf-apply').addEventListener('click', () => {
            const octree = structure?.value === 'octree';
            const n = Math.max(
                1,
                Number($<HTMLInputElement>('op-sdf-resolution').value) || 32,
            );
            this.push({
                ...OP_DEFAULTS.computeSdf,
                structure: octree ? 'octree' : 'voxel',
                resolution: octree ? [] : [n, n, n],
                rootResolution: Math.max(
                    1,
                    Number($<HTMLInputElement>('op-sdf-root').value) || 8,
                ),
                maxDepth: Math.max(0, Number($<HTMLInputElement>('op-sdf-depth').value) || 3),
            });
        });
    }

    /** Offer the mesh's own point arrays; a cell field has no derivative. */
    private setGradientArrays(names: string[]): void {
        const select = maybe<HTMLSelectElement>('op-grad-array');
        if (!select) return;
        setOptions(
            select,
            names.map((name) => ({ value: name, label: name })),
            names[0] ?? ''
        );
        const apply = maybe<HTMLButtonElement>('op-grad-apply');
        if (apply) apply.disabled = names.length === 0;
    }

    private resetSectionSlider(): void {
        const slider = maybe<HTMLInputElement>('op-section-position');
        if (slider) slider.value = '50';
        const flip = maybe<HTMLInputElement>('op-section-flip');
        if (flip) flip.checked = false;
    }

    private renderBakedNotice(): void {
        const note = document.createElement('p');
        note.className = 'hint';
        note.textContent =
            'This is an offline page. Operations were computed when it was written.';
        this.chips.replaceChildren(note);
        for (const el of this.root.querySelectorAll('button, input, select')) {
            const control = el as HTMLButtonElement;
            control.disabled = true;
            control.title = 'Available in the online viewer';
        }
    }
}

/** The colour-by options a baked (offline) page can offer. */
export function bakedArrays(arrayLabels: string[]): string[] {
    return arrayLabels.filter((label) => label.startsWith('quality:'));
}

export { setOptions };
