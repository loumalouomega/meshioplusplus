/**
 * The data table's model: columns over rows, sorting and CSV.
 *
 * Pure -- the panel (features/analysis.ts) only draws what this describes, so
 * ordering, naming and export are unit-tested without a DOM.
 *
 * The idea (windowed rows, sort, column picker, CSV) is the Kratos MDPA Preview
 * extension's Data table; this is a typed-array version of it.
 */

export interface Column {
    /** Stable key, unique within the table. */
    key: string;
    name: string;
    /** Value in row `row`; non-finite values print blank. */
    get(row: number): number;
}

export interface Table {
    rowCount: number;
    columns: Column[];
}

export interface ArraySpec {
    name: string;
    values: ArrayLike<number>;
    components: number;
}

const AXES = ['x', 'y', 'z'];

/** Column name for component `k` of a `n`-component array. */
export function componentName(name: string, n: number, k: number): string {
    return n === 1 ? name : n <= 3 ? `${name}_${AXES[k]}` : `${name}_${k}`;
}

/** The columns of one array: a column per component, plus the magnitude of a vector. */
export function arrayColumns(spec: ArraySpec): Column[] {
    const { name, values, components: n } = spec;
    const out: Column[] = [];
    for (let k = 0; k < n; k++) {
        out.push({ key: `${name}#${k}`, name: componentName(name, n, k), get: (row) => Number(values[row * n + k]) });
    }
    if (n > 1) {
        out.push({
            key: `${name}#mag`,
            name: `|${name}|`,
            get: (row) => {
                let s = 0;
                for (let k = 0; k < n; k++) s += Number(values[row * n + k]) ** 2;
                return Math.sqrt(s);
            },
        });
    }
    return out;
}

export function makeTable(rowCount: number, leading: Column[], arrays: ArraySpec[]): Table {
    return { rowCount, columns: [...leading, ...arrays.flatMap(arrayColumns)] };
}

/** Row indices ordered by a column; non-finite values go last either way. Stable. */
export function sortedOrder(table: Table, columnIndex: number, direction: 'asc' | 'desc'): Uint32Array {
    const col = table.columns[columnIndex];
    const n = table.rowCount;
    const order = new Uint32Array(n);
    for (let i = 0; i < n; i++) order[i] = i;
    if (!col) return order;
    const keys = new Float64Array(n);
    for (let i = 0; i < n; i++) keys[i] = col.get(i);
    const sign = direction === 'asc' ? 1 : -1;
    return order.sort((a, b) => {
        const x = keys[a];
        const y = keys[b];
        const fx = Number.isFinite(x);
        const fy = Number.isFinite(y);
        if (fx !== fy) return fx ? -1 : 1;
        if (!fx) return a - b;
        return x === y ? a - b : x < y ? -sign : sign;
    });
}

/** Rows (original indices) whose `column` lies in `[lo, hi]`. */
export function filterRows(table: Table, columnIndex: number, lo: number, hi: number): Uint32Array {
    const col = table.columns[columnIndex];
    const out: number[] = [];
    if (col) {
        for (let i = 0; i < table.rowCount; i++) {
            const v = col.get(i);
            if (Number.isFinite(v) && v >= lo && v <= hi) out.push(i);
        }
    }
    return Uint32Array.from(out);
}

export function csvEscape(text: string): string {
    return /[",\n\r]/.test(text) ? `"${text.replace(/"/g, '""')}"` : text;
}

/** Shortest decimal that round-trips a double; blank when not finite. */
export function formatCell(v: number): string {
    return Number.isFinite(v) ? String(v) : '';
}

/**
 * The table as CSV, in chunks (a multi-million-row table is past a string's
 * length limit). `rows` selects and orders the rows; `columns` the columns.
 */
export function* csvChunks(
    table: Table,
    rows: ArrayLike<number>,
    columns: number[],
    chunkRows = 5000
): Generator<string> {
    const cols = columns.map((i) => table.columns[i]).filter((c): c is Column => !!c);
    yield `${cols.map((c) => csvEscape(c.name)).join(',')}\n`;
    for (let start = 0; start < rows.length; start += chunkRows) {
        const lines: string[] = [];
        const end = Math.min(rows.length, start + chunkRows);
        for (let i = start; i < end; i++) lines.push(cols.map((c) => formatCell(c.get(rows[i]))).join(','));
        yield `${lines.join('\n')}\n`;
    }
}
