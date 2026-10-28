/**
 * The Groups panel: the file's named regions and cell blocks, plus the groups
 * you make from a selection.
 *
 * Groups come from the worker (`getRegions`), already mapped onto the drawn
 * surface (see `groups/regionMap.ts`). Each is a set of surface cells that can
 * be shown or hidden, isolated, selected, framed, recoloured, and used to
 * colour the whole surface. User groups (saved selections) work everywhere,
 * including the offline page, where there is no worker to read regions with.
 */
import type { Visibility } from '../app/visibility';
import {
    type RegionInfo,
    groupLabels,
    groupMask,
    smallestFirst,
} from '../groups/regionMap';
import type { Legend } from '../render/legend';
import type { Renderer } from '../render/renderer';
import type { Mask } from '../selection/mask';
import { idsOf, type Combine } from '../selection/mask';
import type { GroupsState, Vector3 } from '../types';
import type { RegionsResult } from '../worker/protocol';
import { glyph } from '../ui/glyphs';
import { $, show } from '../ui/dom';
import { cssColor } from '../ui/theme';

export type GroupList = 'regions' | 'blocks';

interface GroupEntry {
    id: number;
    list: GroupList;
    info: RegionInfo;
    user: boolean;
    color: Vector3;
    visible: boolean;
}

export interface GroupsOptions {
    renderer: Renderer;
    visibility: Visibility;
    legend: Legend;
    palette: () => Vector3[];
    selectMask: (mask: Mask, how: Combine) => void;
    /** The regions reader, or null where there is no worker (the offline page). */
    fetch: () => ((epoch: number, numSurfaceCells: number) => Promise<RegionsResult>) | null;
    /** The colour-by menu needs rebuilding (groups appeared or went away). */
    onListChanged: () => void;
    onState: (state: GroupsState) => void;
}

const COLOR_KEYS: Record<GroupList, string> = { regions: 'group:regions', blocks: 'group:blocks' };

/** A palette colour for group `i`, darkened for each full cycle so colours stay distinct. */
export function groupColor(palette: Vector3[], i: number): Vector3 {
    const base = palette[i % palette.length] ?? [0.5, 0.5, 0.5];
    const cycle = Math.floor(i / palette.length);
    const k = Math.max(0.35, 1 - 0.2 * cycle);
    return [base[0] * k, base[1] * k, base[2] * k];
}

function toHex(color: Vector3): string {
    const h = (v: number) =>
        Math.round(Math.min(1, Math.max(0, v)) * 255)
            .toString(16)
            .padStart(2, '0');
    return `#${h(color[0])}${h(color[1])}${h(color[2])}`;
}

function fromHex(hex: string): Vector3 {
    const n = parseInt(hex.slice(1), 16);
    return [((n >> 16) & 255) / 255, ((n >> 8) & 255) / 255, (n & 255) / 255];
}

export class GroupsController {
    private readonly o: GroupsOptions;
    private groups: GroupEntry[] = [];
    private nextId = 1;
    private epoch = 0;
    private numCells = 0;
    private status: GroupsState['status'] = 'none';
    private note = '';
    /** Per list: the group currently isolated, if any. */
    private isolated: Record<GroupList, number | null> = { regions: null, blocks: null };
    private coloring: GroupList | null = null;

    constructor(options: GroupsOptions) {
        this.o = options;
        $('groups-show-all').addEventListener('click', () => this.showAll());
    }

    // --- surface lifecycle -------------------------------------------------------

    /**
     * A surface was (re)loaded. `idsPreserved` says whether its cells are still
     * the file's, in order (no operation that changes them): only then do the
     * file's groups apply. `keep` keeps what was already loaded.
     */
    attachSurface(numCells: number, keep: boolean, idsPreserved: boolean): void {
        const sameSurface = keep && numCells === this.numCells;
        this.numCells = numCells;
        this.coloring = null;

        // User groups are masks over these exact cells.
        if (!sameSurface) this.groups = this.groups.filter(() => false);
        else this.groups = this.groups.filter((g) => !g.user || (g.info.surfaceCells?.length ?? 0) > 0);

        if (!idsPreserved) {
            this.groups = this.groups.filter((g) => g.user && sameSurface);
            this.epoch++;
            this.status = 'unavailable';
            this.note =
                'Groups describe the file as opened. They are hidden while an operation that changes the cells is applied.';
            this.resetVisibility();
            this.render();
            return;
        }

        const haveFileGroups = this.groups.some((g) => !g.user);
        if (sameSurface && haveFileGroups) {
            this.resetVisibility();
            this.applyVisibility();
            this.render();
            return;
        }

        this.groups = this.groups.filter((g) => g.user);
        this.resetVisibility();
        const fetcher = this.o.fetch();
        if (!fetcher || numCells === 0) {
            this.status = this.groups.length ? 'ready' : 'none';
            this.note = '';
            this.render();
            return;
        }

        const epoch = ++this.epoch;
        this.status = 'loading';
        this.note = 'Reading groups…';
        this.render();
        fetcher(epoch, numCells)
            .then((result) => {
                if (epoch !== this.epoch || result.numSurfaceCells !== this.numCells) return;
                this.receive(result);
            })
            .catch((e: Error) => {
                if (epoch !== this.epoch) return;
                this.status = 'failed';
                this.note = `Could not read groups: ${e.message}`;
                this.render();
            });
    }

    private receive(result: RegionsResult): void {
        if (!result.ok) {
            this.status = result.skipped ? 'unavailable' : 'failed';
            this.note = result.skipped ?? `Groups unavailable: ${result.reason ?? 'unknown reason'}`;
            this.render();
            return;
        }
        const palette = this.o.palette();
        const users = this.groups.filter((g) => g.user);
        const fresh: GroupEntry[] = [];
        let regionRank = users.length;
        let blockRank = 0;
        for (const info of result.groups) {
            const list: GroupList = info.source === 'block' ? 'blocks' : 'regions';
            const rank = list === 'regions' ? regionRank++ : blockRank++;
            fresh.push({
                id: this.nextId++,
                list,
                info,
                user: false,
                color: groupColor(palette, rank),
                visible: true,
            });
        }
        this.groups = [...fresh, ...users];
        const regions = fresh.filter((g) => g.list === 'regions').length;
        this.status = 'ready';
        this.note = regions === 0 ? 'This file defines no named regions; its cell blocks are listed below.' : '';
        this.render();
        this.o.onListChanged();
    }

    // --- user groups ------------------------------------------------------------------

    addUserGroup(name: string, mask: Mask, color: Vector3): void {
        const cells = idsOf(mask);
        if (cells.length === 0) return;
        const unique = this.groups.filter((g) => g.list === 'regions' && g.info.name.startsWith(name)).length;
        const label = unique ? `${name} (${unique + 1})` : name;
        this.groups.push({
            id: this.nextId++,
            list: 'regions',
            user: true,
            color,
            visible: true,
            info: {
                name: label,
                kind: 'cell',
                source: 'region',
                dim: 2,
                tag: -1,
                numEntries: cells.length,
                mapped: cells.length,
                surfaceCells: cells,
                surfacePoints: null,
            },
        });
        if (this.status === 'none') this.status = 'ready';
        this.render();
        this.o.onListChanged();
    }

    /** Groups Refine can be limited to: file regions of cells (by name) and saved groups (by cells). */
    listForRefine(): { name: string; user: boolean; surfaceCells: Uint32Array }[] {
        return this.groups
            .filter((g) => g.list === 'regions' && g.info.kind === 'cell' && g.info.surfaceCells)
            .map((g) => ({ name: g.info.name, user: g.user, surfaceCells: g.info.surfaceCells as Uint32Array }));
    }

    /** The user groups, for exporting them with the mesh. */
    userGroups(): { name: string; surfaceCells: Uint32Array }[] {
        return this.groups
            .filter((g) => g.user && g.info.surfaceCells)
            .map((g) => ({ name: g.info.name, surfaceCells: g.info.surfaceCells as Uint32Array }));
    }

    // --- visibility -----------------------------------------------------------------------

    private inList(list: GroupList): GroupEntry[] {
        return this.groups.filter((g) => g.list === list);
    }

    private resetVisibility(): void {
        this.isolated = { regions: null, blocks: null };
        for (const g of this.groups) g.visible = true;
        this.o.visibility.set('groups:regions', null);
        this.o.visibility.set('groups:blocks', null);
    }

    private applyVisibility(): void {
        for (const list of ['regions', 'blocks'] as const) {
            const entries = this.inList(list);
            const n = this.numCells;
            const iso = this.isolated[list];
            let layer: Mask | null = null;
            if (iso !== null) {
                const g = entries.find((e) => e.id === iso);
                layer = g ? groupMask(g.info, n) : null;
            } else if (entries.some((g) => !g.visible)) {
                layer = new Uint8Array(n).fill(1);
                for (const g of entries) {
                    if (g.visible || !g.info.surfaceCells) continue;
                    for (const c of g.info.surfaceCells) layer[c] = 0;
                }
            }
            this.o.visibility.set(`groups:${list}`, layer);
        }
    }

    showAll(): void {
        this.resetVisibility();
        this.render();
    }

    private toggle(g: GroupEntry): void {
        if (this.isolated[g.list] !== null) {
            // Leaving isolate: everything comes back, and the click is spent.
            this.isolated[g.list] = null;
            for (const e of this.inList(g.list)) e.visible = true;
        } else {
            g.visible = !g.visible;
        }
        this.applyVisibility();
        this.render();
    }

    private isolate(g: GroupEntry): void {
        this.isolated[g.list] = this.isolated[g.list] === g.id ? null : g.id;
        this.applyVisibility();
        this.render();
    }

    // --- colouring ---------------------------------------------------------------------------

    /** Colour-by entries for the groups that exist. */
    colorOptions(): { value: string; label: string }[] {
        const out: { value: string; label: string }[] = [];
        if (this.inList('regions').some((g) => g.info.surfaceCells))
            out.push({ value: COLOR_KEYS.regions, label: 'Groups: regions' });
        if (this.inList('blocks').some((g) => g.info.surfaceCells))
            out.push({ value: COLOR_KEYS.blocks, label: 'Groups: cell blocks' });
        return out;
    }

    /** Whether a colour-by value belongs to this panel. */
    owns(value: string): boolean {
        return value === COLOR_KEYS.regions || value === COLOR_KEYS.blocks;
    }

    applyColoring(value: string): void {
        const list: GroupList = value === COLOR_KEYS.blocks ? 'blocks' : 'regions';
        const entries = this.inList(list).filter((g) => g.info.surfaceCells);
        if (entries.length === 0) return;
        this.coloring = list;
        const infos = entries.map((g) => g.info);
        const labels = groupLabels(infos, this.numCells, smallestFirst(infos));
        this.o.renderer.setLabelColoring(
            labels,
            entries.map((g) => g.color)
        );
        this.o.legend.showCategories(
            list === 'regions' ? 'Regions' : 'Cell blocks',
            entries.map((g) => ({ label: g.info.name, color: g.color }))
        );
    }

    /** Some other colour-by was chosen. */
    clearColoring(): void {
        this.coloring = null;
    }

    // --- panel ----------------------------------------------------------------------------------

    private publish(): void {
        this.o.onState({
            status: this.status,
            groups: this.groups.map((g) => ({
                id: g.id,
                name: g.info.name,
                list: g.list,
                cells: g.info.surfaceCells?.length ?? 0,
                visible: this.isolated[g.list] !== null ? this.isolated[g.list] === g.id : g.visible,
                user: g.user,
            })),
        });
    }

    private render(): void {
        const section = $('groups-section');
        show(section, this.groups.length > 0 || this.status === 'loading' || this.status === 'unavailable' || this.status === 'failed');
        $('groups-note').textContent = this.note;
        $('groups-note').hidden = !this.note;
        for (const list of ['regions', 'blocks'] as const) {
            const entries = this.inList(list);
            const head = $(`groups-${list}-head`);
            const ul = $(`groups-${list}`);
            head.hidden = entries.length === 0;
            ul.hidden = entries.length === 0;
            ul.replaceChildren(...entries.map((g) => this.row(g)));
        }
        $('groups-show-all').hidden = this.groups.length === 0;
        if (this.coloring) this.applyColoring(COLOR_KEYS[this.coloring]);
        this.publish();
    }

    private row(g: GroupEntry): HTMLElement {
        const li = document.createElement('li');
        li.dataset.id = String(g.id);
        const cells = g.info.surfaceCells?.length ?? 0;
        const shown = this.isolated[g.list] !== null ? this.isolated[g.list] === g.id : g.visible;

        const eye = document.createElement('input');
        eye.type = 'checkbox';
        eye.className = 'g-eye';
        eye.checked = shown;
        eye.setAttribute('aria-label', `Show ${g.info.name}`);
        eye.addEventListener('change', () => this.toggle(g));

        const color = document.createElement('input');
        color.type = 'color';
        color.className = 'g-color';
        color.value = toHex(g.color);
        color.title = 'Colour';
        color.addEventListener('input', () => {
            g.color = fromHex(color.value);
            if (this.coloring === g.list) this.applyColoring(COLOR_KEYS[g.list]);
        });

        const name = document.createElement('span');
        name.className = 'g-name';
        name.textContent = g.info.name;
        const tip = g.user
            ? `${cells.toLocaleString('en-US')} surface cells (saved selection)`
            : g.info.kind === 'point'
              ? `${g.info.numEntries.toLocaleString('en-US')} points in the file, ${(g.info.surfacePoints?.length ?? 0).toLocaleString('en-US')} on the surface`
              : `${g.info.numEntries.toLocaleString('en-US')} in the file, ${cells.toLocaleString('en-US')} surface cells${g.info.note ? ` — ${g.info.note}` : ''}`;
        name.title = tip;

        const count = document.createElement('span');
        count.className = 'g-count';
        count.textContent = g.info.kind === 'point' ? `${g.info.numEntries} pts` : cells.toLocaleString('en-US');

        li.append(eye, color, name, count);
        if (g.info.surfaceCells) {
            const btn = (icon: Parameters<typeof glyph>[0], title: string, run: () => void) => {
                const b = document.createElement('button');
                b.type = 'button';
                b.className = 'g-action';
                b.title = title;
                b.innerHTML = glyph(icon);
                b.addEventListener('click', run);
                return b;
            };
            li.append(
                btn('pointer', 'Select its cells', () => this.o.selectMask(groupMask(g.info, this.numCells), 'replace')),
                btn('isolate', 'Show only this group', () => this.isolate(g)),
                btn('maximize', 'Zoom to this group', () => this.o.renderer.frameCells(groupMask(g.info, this.numCells)))
            );
        } else {
            li.classList.add('g-disabled');
        }
        return li;
    }
}

/** Swatch helper shared with tests of the panel. */
export { cssColor as groupCss };
