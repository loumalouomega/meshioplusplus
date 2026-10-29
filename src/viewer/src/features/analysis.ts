/**
 * The data table and the statistics (histogram) panel.
 *
 * Both read the full surface that is on screen, so they work on the offline
 * page too. The table is windowed (only the visible rows are in the DOM) and
 * sorts and exports through the pure model in `analysis/table.ts`; a clicked
 * row, or a clicked histogram bar, becomes a selection.
 *
 * The table's shape -- cells/points, sortable columns, a column picker, CSV --
 * is the Kratos MDPA Preview extension's Data table.
 */
import { registerAction } from '../app/actions';
import { binRange, histogram, type Histogram } from '../analysis/histogram';
import {
    type Column,
    type Table,
    csvChunks,
    formatCell,
    makeTable,
    sortedOrder,
} from '../analysis/table';
import { colorKey, INTERNAL_PREFIXES } from '../render/scalars';
import type { Renderer } from '../render/renderer';
import type { Mask } from '../selection/mask';
import { maskFromIds, thresholdCellValues, thresholdPointValues } from '../selection/mask';
import type { AnalysisState, ArrayEntry } from '../types';
import { $, setOptions, show } from '../ui/dom';

const ROW_H = 22;
const COL_W = 104;
const OVERSCAN = 6;
/** Sorting reads every row of a column; above this the sort is refused rather than made slow. */
const MAX_SORT_ROWS = 3_000_000;

type Kind = 'cells' | 'points';

export interface AnalysisOptions {
    renderer: Renderer;
    /** Replace the active selection set with these cells. */
    select: (mask: Mask) => void;
    /** Save text as a download. */
    download: (name: string, text: string, mime: string) => void;
    onState: (state: AnalysisState) => void;
    fileBase: () => string;
}

export interface AnalysisController {
    attachSurface(): void;
}

export function initAnalysis(opts: AnalysisOptions): AnalysisController {
    const { renderer } = opts;

    // ------------------------------------------------------------------ table
    let kind: Kind = 'cells';
    let table: Table | null = null;
    let order: Uint32Array | null = null;
    let sortCol = -1;
    let direction: 'asc' | 'desc' = 'asc';
    let hidden = new Set<string>();
    let open = false;

    let hist: Histogram | null = null;
    let histArray: string | null = null;
    let histEntry: ArrayEntry | null = null;
    let entries: ArrayEntry[] = [];

    function publish(): void {
        opts.onState({
            table: {
                open,
                kind,
                rows: table?.rowCount ?? 0,
                sortedBy: sortCol >= 0 ? (table?.columns[sortCol]?.name ?? null) : null,
                direction,
            },
            histogram: { array: histArray, counts: hist ? Array.from(hist.counts) : [] },
        });
    }

    function buildTable(): Table | null {
        const data = renderer.data;
        const geometry = renderer.surfaceGeometry();
        if (!data || !geometry) return null;
        const attrs = kind === 'cells' ? data.getCellData() : data.getPointData();
        const arrays = [];
        for (let i = 0; i < attrs.getNumberOfArrays(); i++) {
            const a = attrs.getArrayByIndex(i);
            if (!a) continue;
            const name = a.getName();
            // The provenance array is useful in a table; the other plumbing is not.
            if (name !== 'surface:parent_cell' && INTERNAL_PREFIXES.some((p) => name.startsWith(p))) continue;
            arrays.push({ name: name.replace(/^surface:/, ''), values: a.getData(), components: a.getNumberOfComponents() });
        }
        if (kind === 'cells') {
            const { offsets } = geometry.topo;
            const lead: Column[] = [
                { key: 'id', name: 'cell', get: (i) => i },
                { key: 'corners', name: 'corners', get: (i) => offsets[i + 1] - offsets[i] },
            ];
            return makeTable(geometry.topo.numCells, lead, arrays);
        }
        const coords = () => data.getPoints().getData();
        const lead: Column[] = [
            { key: 'id', name: 'point', get: (i) => i },
            { key: 'x', name: 'x', get: (i) => coords()[i * 3] },
            { key: 'y', name: 'y', get: (i) => coords()[i * 3 + 1] },
            { key: 'z', name: 'z', get: (i) => coords()[i * 3 + 2] },
        ];
        return makeTable(data.getNumberOfPoints(), lead, arrays);
    }

    function visibleColumns(): number[] {
        return (table?.columns ?? []).map((c, i) => (hidden.has(c.key) ? -1 : i)).filter((i) => i >= 0);
    }

    function rowAt(slot: number): number {
        return order ? order[slot] : slot;
    }

    const scroll = $('tbl-scroll');
    const head = $('tbl-head');
    const rowsEl = $('tbl-rows');
    const body = $('tbl-body');

    function renderHead(): void {
        const cols = visibleColumns();
        head.style.width = `${cols.length * COL_W}px`;
        head.replaceChildren(
            ...cols.map((ci) => {
                const col = table!.columns[ci];
                const b = document.createElement('button');
                b.type = 'button';
                b.className = 'tbl-th';
                b.dataset.col = col.key;
                b.title = `Sort by ${col.name}`;
                const arrow = ci === sortCol ? (direction === 'asc' ? ' ▲' : ' ▼') : '';
                b.textContent = col.name + arrow;
                b.setAttribute('aria-sort', ci === sortCol ? (direction === 'asc' ? 'ascending' : 'descending') : 'none');
                b.addEventListener('click', () => sortBy(ci));
                return b;
            })
        );
    }

    function renderRows(): void {
        if (!table) return;
        const cols = visibleColumns();
        const total = table.rowCount;
        body.style.height = `${total * ROW_H}px`;
        const top = Math.max(0, scroll.scrollTop - head.offsetHeight);
        const first = Math.max(0, Math.floor(top / ROW_H) - OVERSCAN);
        const last = Math.min(total, Math.ceil((top + scroll.clientHeight) / ROW_H) + OVERSCAN);
        const frag = document.createDocumentFragment();
        for (let slot = first; slot < last; slot++) {
            const row = rowAt(slot);
            const div = document.createElement('div');
            div.className = 'tbl-row';
            div.dataset.row = String(row);
            div.style.top = `${slot * ROW_H}px`;
            div.style.width = `${cols.length * COL_W}px`;
            for (const ci of cols) {
                const cell = document.createElement('span');
                cell.textContent = formatCell(table.columns[ci].get(row));
                div.append(cell);
            }
            frag.append(div);
        }
        rowsEl.replaceChildren(frag);
    }

    function rebuild(): void {
        table = open ? buildTable() : null;
        order = null;
        sortCol = -1;
        $('tbl-count').textContent = table ? `${table.rowCount.toLocaleString('en-US')} ${kind}` : '';
        if (table) {
            renderColumnPicker();
            renderHead();
            renderRows();
        }
        publish();
    }

    function renderColumnPicker(): void {
        $('tbl-column-list').replaceChildren(
            ...table!.columns.map((c) => {
                const label = document.createElement('label');
                label.className = 'check';
                const box = document.createElement('input');
                box.type = 'checkbox';
                box.checked = !hidden.has(c.key);
                box.addEventListener('change', () => {
                    if (box.checked) hidden.delete(c.key);
                    else hidden.add(c.key);
                    renderHead();
                    renderRows();
                });
                label.append(box, document.createTextNode(` ${c.name}`));
                return label;
            })
        );
    }

    function sortBy(ci: number): void {
        if (!table) return;
        if (table.rowCount > MAX_SORT_ROWS) return;
        direction = ci === sortCol && direction === 'asc' ? 'desc' : 'asc';
        sortCol = ci;
        order = sortedOrder(table, ci, direction);
        scroll.scrollTop = 0;
        renderHead();
        renderRows();
        publish();
    }

    function selectRow(row: number): void {
        const geometry = renderer.surfaceGeometry();
        if (!geometry) return;
        const n = geometry.topo.numCells;
        const ids = kind === 'cells' ? [row] : geometry.pointCells.cells.subarray(
            geometry.pointCells.start[row],
            geometry.pointCells.start[row + 1]
        );
        const mask = maskFromIds(n, Array.from(ids));
        opts.select(mask);
        renderer.frameCells(mask);
    }

    function setOpen(next: boolean): void {
        open = next;
        show($('table-panel'), next);
        $('tool-table').setAttribute('aria-pressed', String(next));
        if (next) rebuild();
        else {
            table = null;
            publish();
        }
        renderer.resize();
    }

    scroll.addEventListener('scroll', renderRows);
    rowsEl.addEventListener('click', (e) => {
        const div = (e.target as HTMLElement).closest<HTMLElement>('.tbl-row');
        if (div) selectRow(Number(div.dataset.row));
    });
    for (const b of document.querySelectorAll<HTMLButtonElement>('#tbl-kind [data-kind]')) {
        b.addEventListener('click', () => {
            kind = b.dataset.kind === 'points' ? 'points' : 'cells';
            for (const o of document.querySelectorAll<HTMLElement>('#tbl-kind [data-kind]')) {
                o.setAttribute('aria-pressed', String(o === b));
            }
            hidden = new Set();
            rebuild();
        });
    }
    $('tbl-close').addEventListener('click', () => setOpen(false));
    $('tool-table').addEventListener('click', () => setOpen(!open));
    $('tbl-csv').addEventListener('click', () => {
        if (!table) return;
        const rows = order ?? Uint32Array.from({ length: table.rowCount }, (_, i) => i);
        const text = [...csvChunks(table, rows, visibleColumns())].join('');
        opts.download(`${opts.fileBase()}-${kind}.csv`, text, 'text/csv');
    });
    registerAction({
        id: 'view.table',
        label: 'Show or hide the data table',
        shortcut: 'T',
        group: 'Tools',
        run: () => $('tool-table').click(),
    });

    // ------------------------------------------------------------- statistics
    const histEl = $('hist');

    function entryValues(entry: ArrayEntry): { values: ArrayLike<number>; components: number } | null {
        const data = renderer.data;
        if (!data) return null;
        const attrs = entry.location === 'cell' ? data.getCellData() : data.getPointData();
        const a = attrs.getArrayByName(entry.name);
        return a ? { values: a.getData(), components: a.getNumberOfComponents() } : null;
    }

    function drawHistogram(): void {
        const entry = histEntry;
        const v = entry ? entryValues(entry) : null;
        if (!entry || !v) {
            hist = null;
            histEl.replaceChildren();
            $('hist-summary').textContent = '';
            $('hist-min').textContent = '';
            $('hist-max').textContent = '';
            publish();
            return;
        }
        const bins = Number($<HTMLSelectElement>('stats-bins').value) || 20;
        hist = histogram(v.values, v.components, { component: entry.component, bins });
        const top = Math.max(1, ...hist.counts);
        histEl.replaceChildren(
            ...Array.from(hist.counts, (count, i) => {
                const bar = document.createElement('button');
                bar.type = 'button';
                bar.className = 'hist-bar';
                const [a, b] = binRange(hist!, i);
                bar.title = `${a.toPrecision(4)} … ${b.toPrecision(4)}: ${count.toLocaleString('en-US')} ${entry.location === 'cell' ? 'cells' : 'points'}`;
                bar.style.setProperty('--h', `${Math.round((count / top) * 100)}%`);
                bar.addEventListener('click', () => selectBin(i));
                return bar;
            })
        );
        const f = (x: number) => String(Number(x.toPrecision(5)));
        $('hist-min').textContent = f(hist.min);
        $('hist-max').textContent = f(hist.max);
        $('hist-summary').innerHTML =
            `<div><span>Count</span><b>${hist.count.toLocaleString('en-US')}</b></div>` +
            `<div><span>Mean</span><b>${f(hist.mean)}</b></div>` +
            `<div><span>Std</span><b>${f(hist.std)}</b></div>` +
            (hist.nonFinite ? `<div><span>Non-finite</span><b>${hist.nonFinite.toLocaleString('en-US')}</b></div>` : '');
        publish();
    }

    function selectBin(i: number): void {
        const entry = histEntry;
        const geometry = renderer.surfaceGeometry();
        const v = entry ? entryValues(entry) : null;
        if (!hist || !entry || !geometry || !v) return;
        const [lo, hi] = binRange(hist, i);
        const mask =
            entry.location === 'cell'
                ? thresholdCellValues(v.values, v.components, entry.component, lo, hi, geometry.topo.numCells)
                : thresholdPointValues(v.values, v.components, entry.component, lo, hi, geometry.topo, 'any');
        opts.select(mask);
    }

    $('stats-array').addEventListener('change', () => {
        const key = $<HTMLSelectElement>('stats-array').value;
        histEntry = entries.find((e) => colorKey(e) === key) ?? null;
        histArray = histEntry?.label ?? null;
        drawHistogram();
    });
    $('stats-bins').addEventListener('change', drawHistogram);
    // Follow the colour-by choice, since that is usually the array being looked at.
    $('color-by').addEventListener('change', () => {
        const key = $<HTMLSelectElement>('color-by').value;
        if (entries.some((e) => colorKey(e) === key)) {
            $<HTMLSelectElement>('stats-array').value = key;
            $('stats-array').dispatchEvent(new Event('change'));
        }
    });

    return {
        attachSurface(): void {
            entries = renderer.arrays();
            show($('stats-section'), entries.length > 0);
            const keep = histEntry ? colorKey(histEntry) : undefined;
            setOptions(
                $<HTMLSelectElement>('stats-array'),
                entries.map((e) => ({ value: colorKey(e), label: e.label })),
                keep
            );
            const key = $<HTMLSelectElement>('stats-array').value;
            histEntry = entries.find((e) => colorKey(e) === key) ?? null;
            histArray = histEntry?.label ?? null;
            drawHistogram();
            if (open) rebuild();
        },
    };
}
