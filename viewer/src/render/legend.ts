/**
 * The colour legend, drawn in the DOM rather than with `vtkScalarBarActor`.
 *
 * The gradient comes from the same stop table that builds the colour transfer
 * function (see `colormaps.ts`), so the legend cannot drift from the render —
 * and that agreement is unit-tested rather than eyeballed. A DOM legend also
 * gets the page's font and theme for free, where the vtk.js actor renders its
 * own text into the scene at whatever size the camera happens to give it.
 *
 * The range inputs are the reason this is interactive at all: an outlier cell
 * otherwise compresses every interesting value into one colour, and the only
 * fix is to narrow the range by hand.
 */
import type { ScalarRange } from '../types';
import { cssGradient, formatValue } from './colormaps';

export interface LegendCallbacks {
    /** A new explicit range was entered or a rescale button was pressed. */
    onRange(min: number, max: number): void;
    /** Reset to the full range of the array. */
    onRescaleToData(): void;
}

const TICK_COUNT = 5;

export class Legend {
    private readonly root: HTMLElement;
    private readonly bar: HTMLElement;
    private readonly ticks: HTMLElement;
    private readonly title: HTMLElement;
    private readonly note: HTMLElement;
    private readonly minInput: HTMLInputElement;
    private readonly maxInput: HTMLInputElement;

    constructor(root: HTMLElement, callbacks: LegendCallbacks) {
        this.root = root;
        root.classList.add('legend');
        root.hidden = true;
        root.innerHTML = `
            <div class="legend-title"></div>
            <div class="legend-body">
              <div class="legend-bar"></div>
              <div class="legend-ticks"></div>
            </div>
            <div class="legend-range">
              <input type="number" class="legend-min" step="any" aria-label="Range minimum" />
              <input type="number" class="legend-max" step="any" aria-label="Range maximum" />
            </div>
            <button type="button" class="legend-rescale">Rescale to data</button>
            <div class="legend-note"></div>`;

        this.title = must(root, '.legend-title');
        this.bar = must(root, '.legend-bar');
        this.ticks = must(root, '.legend-ticks');
        this.note = must(root, '.legend-note');
        this.minInput = must<HTMLInputElement>(root, '.legend-min');
        this.maxInput = must<HTMLInputElement>(root, '.legend-max');

        const commit = () => {
            const min = Number(this.minInput.value);
            const max = Number(this.maxInput.value);
            // Reject rather than silently swap: a user who typed them the
            // wrong way round wants to know, not to get a different picture.
            if (!Number.isFinite(min) || !Number.isFinite(max) || min >= max) {
                this.note.textContent = 'Minimum must be below maximum.';
                return;
            }
            callbacks.onRange(min, max);
        };
        this.minInput.addEventListener('change', commit);
        this.maxInput.addEventListener('change', commit);
        must(root, '.legend-rescale').addEventListener('click', () =>
            callbacks.onRescaleToData()
        );
    }

    hide(): void {
        this.root.hidden = true;
    }

    /** Show the legend for `label`, coloured by `colormap`, spanning `range`. */
    show(label: string, colormap: string, range: ScalarRange): void {
        this.root.hidden = false;
        this.title.textContent = label;
        this.bar.style.background = cssGradient(colormap, 'to top');
        this.update(range);
    }

    /** Refresh the numbers without rebuilding the gradient. */
    update(range: ScalarRange): void {
        this.minInput.value = String(range.min);
        this.maxInput.value = String(range.max);

        const labels: string[] = [];
        for (let i = TICK_COUNT - 1; i >= 0; i--) {
            const t = i / (TICK_COUNT - 1);
            labels.push(formatValue(range.min + t * (range.max - range.min)));
        }
        this.ticks.replaceChildren(
            ...labels.map((text) => {
                const span = document.createElement('span');
                span.textContent = text;
                return span;
            })
        );

        this.note.textContent =
            range.nanCount > 0
                ? `${range.nanCount.toLocaleString('en-US')} non-finite (grey)`
                : '';
    }

    /** Repaint the gradient after a colormap change. */
    setColormap(colormap: string): void {
        this.bar.style.background = cssGradient(colormap, 'to top');
    }
}

function must<T extends HTMLElement = HTMLElement>(root: HTMLElement, sel: string): T {
    const el = root.querySelector<T>(sel);
    if (!el) throw new Error(`legend: missing ${sel}`);
    return el;
}
