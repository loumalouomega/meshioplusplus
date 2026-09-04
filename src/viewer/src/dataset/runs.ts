/**
 * Run history and comparison: the table of every training run the companion
 * process knows, and an overlay of a few runs' validation curves with their
 * hyperparameters side by side.
 *
 * **Three runs is a validated cap, not a round number.** Overlaid curves
 * cross, so their colours must be distinguishable in *every* pair rather
 * than only between neighbours — and under that rule the categorical
 * palette's first three slots pass on this page's surface while the fourth
 * (yellow beside orange) fails both the colour-vision and the
 * normal-vision floors. Selecting a fourth run is refused with that reason
 * rather than silently cycling a colour.
 */

import { isToolError, type ApiClient } from './api';
import { SERIES_COLORS, renderLineChart, type ChartSeries } from './chart';
import { stableStringify } from './diff';
import type { JobSummary, MetricRow } from './jobs';
import { $, svg$ } from '../ui/dom';

/** How many runs may be overlaid at once (see the module comment). */
export const MAX_COMPARED = 3;

/** The hyperparameters a comparison lines up, in the order shown. */
export const COMPARED_FIELDS: { key: keyof JobSummary; label: string }[] = [
    { key: 'status', label: 'status' },
    { key: 'fields', label: 'input fields' },
    { key: 'target_fields', label: 'targets' },
    { key: 'train_split', label: 'train split' },
    { key: 'valid_split', label: 'valid split' },
    { key: 'epochs', label: 'epochs' },
    { key: 'batch_size', label: 'batch' },
    { key: 'learning_rate', label: 'lr' },
    { key: 'seed', label: 'seed' },
    { key: 'processor_size', label: 'layers' },
    { key: 'hidden_dim', label: 'hidden' },
    { key: 'best_valid_loss', label: 'best valid' },
    { key: 'final_valid_loss', label: 'final valid' },
    { key: 'duration_seconds', label: 'duration' },
];

export type RunSort = 'job_id' | 'best_valid_loss' | 'epochs' | 'duration_seconds';

/** Newest first by default; a numeric column sorts ascending with unknown
 * values last (a run that never reported one is not "the best"). */
export function sortRuns(runs: JobSummary[], key: RunSort, descending: boolean): JobSummary[] {
    const out = [...runs];
    out.sort((a, b) => {
        if (key === 'job_id') return a.job_id.localeCompare(b.job_id);
        const x = a[key] as number | null;
        const y = b[key] as number | null;
        if (x === null || x === undefined) return 1;
        if (y === null || y === undefined) return -1;
        return x - y;
    });
    return descending ? out.reverse() : out;
}

export function filterRuns(runs: JobSummary[], needle: string): JobSummary[] {
    const wanted = needle.trim().toLowerCase();
    if (!wanted) return runs;
    return runs.filter((run) =>
        [run.job_id, run.status, run.manifest ?? '', ...(run.tags ?? []), ...(run.fields ?? [])]
            .join(' ')
            .toLowerCase()
            .includes(wanted),
    );
}

function format(value: unknown): string {
    if (value === null || value === undefined || value === '') return '—';
    if (typeof value === 'number') {
        if (!Number.isFinite(value)) return '—';
        return Number.isInteger(value) ? String(value) : String(Number(value.toPrecision(4)));
    }
    if (Array.isArray(value)) return value.length ? value.join(', ') : '—';
    return String(value);
}

/** Rows of `[label, value per run]`, flagged when the runs disagree. */
export function comparisonRows(
    runs: JobSummary[],
): { label: string; values: string[]; differs: boolean }[] {
    return COMPARED_FIELDS.map(({ key, label }) => {
        const raw = runs.map((run) => run[key]);
        const values = raw.map(format);
        const differs = new Set(raw.map((v) => stableStringify(v))).size > 1;
        return { label, values, differs };
    });
}

function el(tag: string, cls = '', text?: string): HTMLElement {
    const node = document.createElement(tag);
    if (cls) node.className = cls;
    if (text !== undefined) node.textContent = text;
    return node;
}

export interface RunsPanelHooks {
    getApi(): ApiClient | null;
    /** Follow one run in the run panel (the training panel's own view). */
    onOpen(jobId: string): void;
    setState(patch: { runs?: JobSummary[]; compare?: string[] }): void;
    setStatus(message: string): void;
    fail(error: unknown): void;
}

/** The history table plus the comparison overlay. */
export class RunsPanel {
    private runs: JobSummary[] = [];
    private compare: string[] = [];
    private curves = new Map<string, MetricRow[]>();
    private sort: RunSort = 'job_id';
    private descending = true;

    constructor(private readonly hooks: RunsPanelHooks) {
        $('runs-close').addEventListener('click', () => this.close());
        $('runs-filter').addEventListener('input', () => this.render());
        $<HTMLSelectElement>('runs-sort').addEventListener('change', (e) => {
            this.sort = (e.target as HTMLSelectElement).value as RunSort;
            this.render();
        });
        $('runs-order').addEventListener('click', () => {
            this.descending = !this.descending;
            this.render();
        });
    }

    /** Load the run list and show the panel. */
    async open(): Promise<void> {
        const api = this.hooks.getApi();
        if (!api) return;
        const listed = await api.tool<{ jobs: JobSummary[] }>('train_list', {});
        if (isToolError(listed)) throw new Error(listed.error);
        this.runs = listed.jobs;
        this.hooks.setState({ runs: this.runs });
        // Drop selections whose run is gone, and any curve that went with it.
        this.compare = this.compare.filter((id) => this.runs.some((r) => r.job_id === id));
        $('runs-wrap').hidden = false;
        this.render();
        await this.refreshCurves();
    }

    close(): void {
        $('runs-wrap').hidden = true;
    }

    private async refreshCurves(): Promise<void> {
        const api = this.hooks.getApi();
        if (!api) return;
        for (const jobId of this.compare) {
            if (this.curves.has(jobId)) continue;
            const metrics = await api.tool<{ rows: MetricRow[] }>('train_metrics', {
                job_id: jobId,
            });
            if (isToolError(metrics)) continue;
            this.curves.set(jobId, metrics.rows);
        }
        this.render();
    }

    private toggle(jobId: string): void {
        if (this.compare.includes(jobId)) {
            this.compare = this.compare.filter((id) => id !== jobId);
        } else if (this.compare.length >= MAX_COMPARED) {
            this.hooks.setStatus(
                `comparing at most ${MAX_COMPARED} runs: beyond that the curve colours ` +
                    'stop being reliably distinguishable',
            );
            return;
        } else {
            this.compare = [...this.compare, jobId];
        }
        this.hooks.setState({ compare: this.compare });
        this.render();
        void this.refreshCurves().catch(this.hooks.fail);
    }

    private render(): void {
        const needle = $<HTMLInputElement>('runs-filter').value;
        const shown = sortRuns(filterRuns(this.runs, needle), this.sort, this.descending);
        $('runs-count').textContent = `(${shown.length}/${this.runs.length})`;
        $('runs-order').textContent = this.descending ? '↓' : '↑';

        const header = el('tr');
        for (const label of ['', 'run', 'status', 'epochs', 'best valid', 'duration', 'tags']) {
            header.append(el('th', '', label));
        }
        const rows = shown.map((run) => {
            const tr = el('tr');
            if (this.compare.includes(run.job_id)) tr.classList.add('selected');
            const pick = document.createElement('input');
            pick.type = 'checkbox';
            pick.checked = this.compare.includes(run.job_id);
            pick.title = 'Compare this run';
            pick.addEventListener('change', () => this.toggle(run.job_id));
            const cell = el('td');
            cell.append(pick);
            tr.append(cell);
            const open = document.createElement('button');
            open.type = 'button';
            open.className = 'run-open';
            open.textContent = run.job_id;
            open.addEventListener('click', () => this.hooks.onOpen(run.job_id));
            const idCell = el('td');
            idCell.append(open);
            tr.append(idCell);
            const status = el('td');
            status.append(
                el(
                    'span',
                    `badge ${run.status === 'finished' ? 'ok' : run.status === 'running' ? 'none' : 'bad'}`,
                    run.status,
                ),
            );
            tr.append(status);
            tr.append(el('td', '', format(run.epochs)));
            tr.append(el('td', '', format(run.best_valid_loss)));
            tr.append(el('td', '', run.duration_seconds ? `${Math.round(run.duration_seconds)} s` : '—'));
            const tags = el('td');
            for (const tag of run.tags ?? []) tags.append(el('span', 'chip tag', tag));
            tr.append(tags);
            return tr;
        });
        $('runs-table').replaceChildren(header, ...rows);
        this.renderComparison();
    }

    private renderComparison(): void {
        const selected = this.compare
            .map((id) => this.runs.find((r) => r.job_id === id))
            .filter((run): run is JobSummary => !!run);
        const wrap = $('runs-compare');
        if (selected.length < 1) {
            wrap.hidden = true;
            return;
        }
        wrap.hidden = false;
        const series: ChartSeries[] = selected.map((run, i) => ({
            name: run.job_id.slice(-4),
            color: SERIES_COLORS[i % SERIES_COLORS.length]!,
            points: (this.curves.get(run.job_id) ?? [])
                .filter((row) => row.valid_loss !== null)
                .map((row) => [row.epoch, row.valid_loss as number] as [number, number]),
        }));
        renderLineChart(svg$('runs-chart'), series, {
            logY: true,
            xLabel: 'epoch',
            yLabel: 'validation loss',
            tooltip: $('runs-tooltip'),
            formatX: (x) => `epoch ${Math.round(x)}`,
        });
        $('runs-legend').replaceChildren(
            ...selected.map((run, i) => {
                const item = el('span', 'split-item');
                const swatch = el('i', 'split-swatch');
                swatch.style.background = SERIES_COLORS[i % SERIES_COLORS.length]!;
                item.append(swatch, run.job_id);
                return item;
            }),
        );
        const header = el('tr');
        header.append(el('th', '', ''));
        for (const run of selected) header.append(el('th', '', run.job_id.slice(-4)));
        const rows = comparisonRows(selected).map((row) => {
            const tr = el('tr');
            if (row.differs) tr.classList.add('differs');
            tr.append(el('th', '', row.label));
            for (const value of row.values) tr.append(el('td', '', value));
            return tr;
        });
        $('runs-diff').replaceChildren(header, ...rows);
    }
}
