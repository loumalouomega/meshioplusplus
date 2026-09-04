/**
 * The overview depth of the dataset page: one card per manifest found at
 * the workspace root, with entry counts, split balance, tags/groups,
 * last-modified, health badges and a thumbnail.
 *
 * Card construction and sorting are pure (`tests/viewer/unit/overview.test.mjs`);
 * only `renderCards`/`splitBar` touch the DOM.
 *
 * Colours follow the dataviz reference palette's dark-surface steps: split
 * hues are categorical slots assigned in FIXED order by sorted split name
 * (never cycled — a ninth split folds into "other"), so the same split keeps
 * the same colour on every card; badge colours are the fixed status palette.
 */

import type { ServerManifest } from './api';
import { healthBadges } from './health';
import { parseManifest, splitCounts, type Manifest } from './manifest';
import type { ManifestCard, SplitBalance } from './types';

/** Categorical slots 1-8 (dark surface). */
export const SPLIT_COLORS = [
    '#3987e5',
    '#d95926',
    '#199e70',
    '#c98500',
    '#d55181',
    '#008300',
    '#9085e9',
    '#e66767',
];
/** Splits past the eighth slot, and the "unassigned" share. */
export const OTHER_COLOR = '#898781';
export const UNASSIGNED_COLOR = '#383835';

/** The fixed split -> colour assignment over every split name in play. */
export function splitPaletteOf(names: Iterable<string>): Map<string, string> {
    const sorted = [...new Set(names)].filter((n) => n !== '').sort();
    const palette = new Map<string, string>();
    sorted.forEach((name, i) => palette.set(name, SPLIT_COLORS[i] ?? OTHER_COLOR));
    palette.set('', UNASSIGNED_COLOR);
    return palette;
}

export function cardFromManifest(
    path: string,
    manifest: Manifest,
    lastModified: number | null,
): ManifestCard {
    const tags = new Set<string>();
    const groups = new Set<string>();
    for (const entry of manifest.entries) {
        for (const t of entry.tags) tags.add(t);
        if (entry.group) groups.add(entry.group);
    }
    return {
        path,
        name: manifest.name,
        description: manifest.description,
        numEntries: manifest.entries.length,
        splits: Object.fromEntries(splitCounts(manifest)),
        tags: [...tags].sort(),
        groups: [...groups].sort(),
        lastModified,
        parseError: null,
        health: null,
        thumbnail: null,
        sha256: null,
        dirty: false,
        serverPath: null,
        serverOnly: false,
    };
}

/** A card for a manifest only the companion process knows (found by
 * `dataset_find`, not in the picked directory). Its `path` is namespaced so
 * it can never collide with a workspace card. */
export function cardFromServerManifest(m: ServerManifest): ManifestCard {
    return {
        path: `server:${m.relpath}`,
        name: m.name,
        description: null,
        numEntries: m.num_entries,
        splits: { ...m.splits },
        tags: [],
        groups: [],
        lastModified: m.mtime,
        parseError: null,
        health: null,
        thumbnail: null,
        sha256: m.sha256,
        dirty: false,
        serverPath: m.path,
        serverOnly: true,
    };
}

/** A card from a manifest's text; a strict-parse failure is a flagged card,
 * never a missing one — a broken manifest is exactly what the overview
 * should surface. */
export function buildCard(path: string, text: string, lastModified: number | null): ManifestCard {
    try {
        return cardFromManifest(path, parseManifest(text), lastModified);
    } catch (e) {
        return {
            path,
            name: null,
            description: null,
            numEntries: 0,
            splits: {},
            tags: [],
            groups: [],
            lastModified,
            parseError: e instanceof Error ? e.message : String(e),
            health: null,
            thumbnail: null,
            sha256: null,
            dirty: false,
            serverPath: null,
            serverOnly: false,
        };
    }
}

export type CardSort = 'name' | 'entries' | 'modified' | 'health';

/** Worst first for the health sort: broken, bad, unscanned, healthy. */
function healthRank(card: ManifestCard): number {
    if (card.parseError) return 0;
    if (!card.health) return 2;
    return card.health.badEntries.length ? 1 : 3;
}

export function sortCards(cards: ManifestCard[], key: CardSort): ManifestCard[] {
    const out = [...cards];
    const byName = (a: ManifestCard, b: ManifestCard) =>
        (a.name ?? a.path).localeCompare(b.name ?? b.path);
    switch (key) {
        case 'entries':
            out.sort((a, b) => b.numEntries - a.numEntries || byName(a, b));
            break;
        case 'modified':
            out.sort((a, b) => (b.lastModified ?? -1) - (a.lastModified ?? -1) || byName(a, b));
            break;
        case 'health':
            out.sort((a, b) => healthRank(a) - healthRank(b) || byName(a, b));
            break;
        default:
            out.sort(byName);
    }
    return out;
}

/** `{split: count}` -> balance rows (unassigned last), for a card. */
export function balanceOf(splits: Record<string, number>): SplitBalance[] {
    const total = Object.values(splits).reduce((n, c) => n + c, 0);
    const rows = Object.entries(splits).map(([split, count]) => ({
        split,
        count,
        fraction: total ? count / total : 0,
    }));
    rows.sort((x, y) => {
        if (x.split === '') return 1;
        if (y.split === '') return -1;
        return x.split.localeCompare(y.split);
    });
    return rows;
}

function el(tag: string, cls: string, text?: string): HTMLElement {
    const node = document.createElement(tag);
    if (cls) node.className = cls;
    if (text !== undefined) node.textContent = text;
    return node;
}

/** A proportional bar of the splits plus its legend (colour AND label —
 * identity is never colour alone). */
export function splitBar(rows: SplitBalance[], palette: Map<string, string>): HTMLElement {
    const wrap = el('div', 'split-wrap');
    const bar = el('div', 'split-bar');
    const legend = el('div', 'split-legend');
    for (const row of rows) {
        const color = palette.get(row.split) ?? OTHER_COLOR;
        const seg = el('span', 'split-seg');
        seg.style.flex = String(Math.max(row.fraction, 0.02));
        seg.style.background = color;
        seg.title = `${row.split || 'unassigned'}: ${row.count}`;
        bar.append(seg);
        const item = el('span', 'split-item');
        const swatch = el('i', 'split-swatch');
        swatch.style.background = color;
        item.append(swatch, `${row.split || 'unassigned'} ${row.count}`);
        legend.append(item);
    }
    wrap.append(bar, legend);
    return wrap;
}

export function formatModified(ms: number | null): string {
    if (ms === null) return 'modified: unknown';
    const d = new Date(ms);
    return `modified ${d.toLocaleDateString()} ${d.toLocaleTimeString([], {
        hour: '2-digit',
        minute: '2-digit',
    })}`;
}

export interface CardHandlers {
    onOpen(path: string): void;
    onScan(path: string): void;
    /** Present only while a companion process is connected. */
    onScanServer?(path: string): void;
}

export function renderCards(
    container: HTMLElement,
    cards: ManifestCard[],
    palette: Map<string, string>,
    handlers: CardHandlers,
): void {
    const nodes = cards.map((card) => {
        const node = el('article', 'card');
        node.dataset.path = card.path;
        if (card.parseError) node.classList.add('broken');
        if (card.dirty) node.classList.add('dirty');
        if (card.serverOnly) node.classList.add('server-only');

        const thumb = el('div', 'card-thumb');
        if (card.thumbnail) {
            const img = document.createElement('img');
            img.src = card.thumbnail;
            img.alt = '';
            thumb.append(img);
        } else {
            thumb.append(
                el(
                    'span',
                    '',
                    card.parseError
                        ? 'unreadable'
                        : card.serverOnly
                          ? 'on the companion process only'
                          : 'no preview yet',
                ),
            );
        }
        node.append(thumb);

        const body = el('div', 'card-body');
        body.append(el('h3', 'card-title', card.name ?? card.path));
        const sub = [card.serverOnly ? (card.serverPath ?? card.path) : card.path];
        sub.push(`${card.numEntries} entr${card.numEntries === 1 ? 'y' : 'ies'}`);
        sub.push(formatModified(card.lastModified));
        body.append(el('p', 'card-sub', sub.join(' · ')));
        if (card.description) body.append(el('p', 'card-desc', card.description));
        if (card.parseError) {
            body.append(el('p', 'error', card.parseError));
        } else {
            body.append(splitBar(balanceOf(card.splits), palette));
            if (card.tags.length || card.groups.length || card.serverPath) {
                const chips = el('div', 'card-chips');
                for (const t of card.tags) chips.append(el('span', 'chip tag', t));
                for (const g of card.groups) chips.append(el('span', 'chip group', g));
                if (card.serverPath) {
                    chips.append(
                        el('span', 'chip server', card.serverOnly ? 'server only' : 'on server'),
                    );
                }
                body.append(chips);
            }
            const badges = el('div', 'card-badges');
            if (card.health) {
                for (const b of healthBadges(card.health)) {
                    badges.append(el('span', `badge ${b.level}`, b.label));
                }
            } else {
                badges.append(el('span', 'badge none', 'not scanned'));
            }
            if (card.dirty) badges.append(el('span', 'badge warn', 'unsaved edits'));
            body.append(badges);
        }
        const actions = el('div', 'card-actions');
        const open = el('button', 'card-open', 'Open') as HTMLButtonElement;
        open.type = 'button';
        open.disabled = !!card.parseError || card.serverOnly;
        open.title = card.serverOnly ? 'Pick the directory holding this manifest to curate it' : '';
        open.addEventListener('click', () => handlers.onOpen(card.path));
        const scan = el('button', 'card-scan', 'Scan') as HTMLButtonElement;
        scan.type = 'button';
        scan.disabled = !!card.parseError || card.serverOnly;
        scan.addEventListener('click', () => handlers.onScan(card.path));
        actions.append(open, scan);
        if (card.serverPath && handlers.onScanServer && !card.parseError) {
            const remote = el('button', 'card-scan-server', 'Scan (server)') as HTMLButtonElement;
            remote.type = 'button';
            remote.title = 'Scan on the companion process (full quality metrics, nothing staged here)';
            remote.addEventListener('click', () => handlers.onScanServer?.(card.path));
            actions.append(remote);
        }
        body.append(actions);
        node.append(body);
        return node;
    });
    container.replaceChildren(...nodes);
}
