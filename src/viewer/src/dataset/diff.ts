/**
 * Manifest diffing: a structural, entry-keyed comparison of two manifest
 * documents (what changed, per entry, per field) plus a plain line diff of
 * their serialized text. Manifests are hand-editable JSON with a stable
 * serialization (doc/datasets.md), so both views are meaningful: the
 * structural one answers "which cases moved split", the line one shows the
 * exact edit. Pure data apart from `renderDiff`; unit-tested in
 * `tests/viewer/unit/diff.test.mjs`.
 */

import type { Manifest, ManifestEntry } from './manifest';
import type { DiffSummary } from './types';

export interface FieldChange {
    field: string;
    before: unknown;
    after: unknown;
}

export interface EntryChange {
    id: string;
    fields: FieldChange[];
}

export interface ManifestDiff {
    /** Name / Description / Metadata changes at the document level. */
    header: FieldChange[];
    /** Entry ids present only in B / only in A, in B's / A's order. */
    added: string[];
    removed: string[];
    /** Entries in both whose fields differ. */
    changed: EntryChange[];
}

/** JSON with keys sorted at every level, so key order never reads as a change. */
export function stableStringify(value: unknown): string {
    if (Array.isArray(value)) return `[${value.map(stableStringify).join(',')}]`;
    if (value && typeof value === 'object') {
        const obj = value as Record<string, unknown>;
        const keys = Object.keys(obj).sort();
        return `{${keys.map((k) => `${JSON.stringify(k)}:${stableStringify(obj[k])}`).join(',')}}`;
    }
    return JSON.stringify(value) ?? 'null';
}

const ENTRY_FIELDS: (keyof ManifestEntry)[] = [
    'source',
    'split',
    'tags',
    'group',
    'notes',
    'metadata',
];

function fieldChanges<T extends object>(a: T, b: T, fields: (keyof T)[]): FieldChange[] {
    const out: FieldChange[] = [];
    for (const field of fields) {
        if (stableStringify(a[field]) !== stableStringify(b[field])) {
            out.push({ field: String(field), before: a[field], after: b[field] });
        }
    }
    return out;
}

export function diffManifests(a: Manifest, b: Manifest): ManifestDiff {
    const byIdA = new Map(a.entries.map((e) => [e.id, e]));
    const byIdB = new Map(b.entries.map((e) => [e.id, e]));
    const changed: EntryChange[] = [];
    for (const entry of a.entries) {
        const other = byIdB.get(entry.id);
        if (!other) continue;
        const fields = fieldChanges(entry, other, ENTRY_FIELDS);
        if (fields.length) changed.push({ id: entry.id, fields });
    }
    return {
        header: fieldChanges(a, b, ['name', 'description', 'metadata']),
        added: b.entries.filter((e) => !byIdA.has(e.id)).map((e) => e.id),
        removed: a.entries.filter((e) => !byIdB.has(e.id)).map((e) => e.id),
        changed,
    };
}

export function summarizeDiff(labelA: string, labelB: string, diff: ManifestDiff): DiffSummary {
    return {
        a: labelA,
        b: labelB,
        headerChanged: diff.header.length,
        added: diff.added.length,
        removed: diff.removed.length,
        changed: diff.changed.length,
    };
}

export interface DiffLine {
    kind: 'same' | 'add' | 'del';
    text: string;
}

/** Above this many cell comparisons the line diff is skipped (null): the
 * O(n*m) table would not be worth it, and the structural diff still tells
 * the story. Manifests of a few thousand lines stay well under it. */
export const MAX_LINE_DIFF_CELLS = 25_000_000;

/** A longest-common-subsequence line diff of two texts. */
export function diffLines(aText: string, bText: string): DiffLine[] | null {
    const a = aText.split('\n');
    const b = bText.split('\n');
    if (a.length * b.length > MAX_LINE_DIFF_CELLS) return null;
    // Trim the common prefix/suffix first: the usual manifest edit touches a
    // few lines in the middle, so this keeps the table small.
    let start = 0;
    while (start < a.length && start < b.length && a[start] === b[start]) start += 1;
    let endA = a.length;
    let endB = b.length;
    while (endA > start && endB > start && a[endA - 1] === b[endB - 1]) {
        endA -= 1;
        endB -= 1;
    }
    const n = endA - start;
    const m = endB - start;
    // lcs[i][j] = LCS length of a[start+i..endA) and b[start+j..endB).
    const width = m + 1;
    const lcs = new Uint32Array((n + 1) * width);
    for (let i = n - 1; i >= 0; i -= 1) {
        for (let j = m - 1; j >= 0; j -= 1) {
            lcs[i * width + j] =
                a[start + i] === b[start + j]
                    ? lcs[(i + 1) * width + j + 1] + 1
                    : Math.max(lcs[(i + 1) * width + j], lcs[i * width + j + 1]);
        }
    }
    const out: DiffLine[] = [];
    for (let k = 0; k < start; k += 1) out.push({ kind: 'same', text: a[k] });
    let i = 0;
    let j = 0;
    while (i < n && j < m) {
        if (a[start + i] === b[start + j]) {
            out.push({ kind: 'same', text: a[start + i] });
            i += 1;
            j += 1;
        } else if (lcs[(i + 1) * width + j] >= lcs[i * width + j + 1]) {
            out.push({ kind: 'del', text: a[start + i] });
            i += 1;
        } else {
            out.push({ kind: 'add', text: b[start + j] });
            j += 1;
        }
    }
    for (; i < n; i += 1) out.push({ kind: 'del', text: a[start + i] });
    for (; j < m; j += 1) out.push({ kind: 'add', text: b[start + j] });
    for (let k = endA; k < a.length; k += 1) out.push({ kind: 'same', text: a[k] });
    return out;
}

function el(tag: string, cls: string, text?: string): HTMLElement {
    const node = document.createElement(tag);
    node.className = cls;
    if (text !== undefined) node.textContent = text;
    return node;
}

function describe(value: unknown): string {
    if (value === null || value === undefined) return '—';
    if (typeof value === 'string') return value;
    return stableStringify(value);
}

/** Render the structural diff and (when available) the line diff, with
 * unchanged context collapsed to the three lines around each change. */
export function renderDiff(
    container: HTMLElement,
    diff: ManifestDiff,
    lines: DiffLine[] | null,
): void {
    const parts: HTMLElement[] = [];
    if (diff.header.length) {
        const box = el('div', 'diff-section');
        box.append(el('h3', '', 'Document'));
        for (const c of diff.header) {
            box.append(
                el('p', 'diff-chg', `${c.field}: ${describe(c.before)} → ${describe(c.after)}`),
            );
        }
        parts.push(box);
    }
    if (diff.added.length || diff.removed.length || diff.changed.length) {
        const box = el('div', 'diff-section');
        box.append(el('h3', '', 'Entries'));
        for (const id of diff.added) box.append(el('p', 'diff-add', `+ ${id}`));
        for (const id of diff.removed) box.append(el('p', 'diff-del', `− ${id}`));
        for (const change of diff.changed) {
            const item = el('div', 'diff-entry');
            item.append(el('p', 'diff-chg', `~ ${change.id}`));
            for (const c of change.fields) {
                item.append(
                    el(
                        'p',
                        'diff-field',
                        `${c.field}: ${describe(c.before)} → ${describe(c.after)}`,
                    ),
                );
            }
            box.append(item);
        }
        parts.push(box);
    }
    if (!parts.length) parts.push(el('p', 'hint', 'The two manifests are identical.'));
    if (lines && lines.some((l) => l.kind !== 'same')) {
        const pre = el('pre', 'diff-lines');
        const context = 3;
        const keep = new Set<number>();
        lines.forEach((l, i) => {
            if (l.kind === 'same') return;
            for (let k = i - context; k <= i + context; k += 1) keep.add(k);
        });
        let gap = false;
        lines.forEach((l, i) => {
            if (!keep.has(i)) {
                if (!gap) pre.append(el('span', 'diff-gap', '⋯\n'));
                gap = true;
                return;
            }
            gap = false;
            const marker = l.kind === 'add' ? '+' : l.kind === 'del' ? '−' : ' ';
            pre.append(el('span', `diff-${l.kind}`, `${marker} ${l.text}\n`));
        });
        parts.push(pre);
    }
    container.replaceChildren(...parts);
}
