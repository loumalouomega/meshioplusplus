/**
 * Page shell behaviour: the resizable sidebar, its show/hide toggle, and the
 * keyboard-shortcut sheet (generated from the action registry, so it cannot
 * drift from what the keys actually do).
 */
import { listActions } from '../app/actions';
import { $ } from '../ui/dom';

const WIDTH_KEY = 'meshioplusplus.viewer.sidebar';
const MIN_WIDTH = 220;
const MAX_WIDTH = 640;
const DEFAULT_WIDTH = 300;
const STEP = 16;

export function clampSidebarWidth(width: number, viewport = Infinity): number {
    // Never let the sidebar eat the whole viewport on a narrow window.
    const max = Math.min(MAX_WIDTH, Math.max(MIN_WIDTH, viewport - 240));
    return Math.min(max, Math.max(MIN_WIDTH, Math.round(width)));
}

function loadWidth(): number {
    try {
        const stored = Number(localStorage.getItem(WIDTH_KEY));
        if (Number.isFinite(stored) && stored > 0) return stored;
    } catch {
        // storage unavailable
    }
    return DEFAULT_WIDTH;
}

export function initSidebar(onResize: () => void): { toggle(): void } {
    const workspace = $('workspace');
    const sash = $('sidebar-resizer');

    const setWidth = (width: number, persist = true) => {
        const w = clampSidebarWidth(width, window.innerWidth);
        workspace.style.setProperty('--sidebar-w', `${w}px`);
        if (persist) {
            try {
                localStorage.setItem(WIDTH_KEY, String(w));
            } catch {
                // not worth failing over
            }
        }
        onResize();
        return w;
    };
    let current = setWidth(loadWidth(), false);

    sash.addEventListener('pointerdown', (e) => {
        e.preventDefault();
        sash.setPointerCapture(e.pointerId);
        sash.classList.add('dragging');
        const move = (ev: PointerEvent) => {
            const left = workspace.getBoundingClientRect().left;
            current = setWidth(ev.clientX - left, false);
        };
        const up = () => {
            sash.classList.remove('dragging');
            sash.removeEventListener('pointermove', move);
            sash.removeEventListener('pointerup', up);
            sash.removeEventListener('pointercancel', up);
            current = setWidth(current);
        };
        sash.addEventListener('pointermove', move);
        sash.addEventListener('pointerup', up);
        sash.addEventListener('pointercancel', up);
    });
    sash.addEventListener('keydown', (e) => {
        if (e.key === 'ArrowLeft') current = setWidth(current - STEP);
        else if (e.key === 'ArrowRight') current = setWidth(current + STEP);
        else if (e.key === 'Home') current = setWidth(MIN_WIDTH);
        else if (e.key === 'End') current = setWidth(MAX_WIDTH);
        else return;
        e.preventDefault();
    });
    window.addEventListener('resize', () => {
        current = setWidth(current, false);
    });

    return {
        toggle() {
            workspace.classList.toggle('sidebar-hidden');
            onResize();
        },
    };
}

/** Fill and open the shortcut sheet from the registered actions. */
export function showShortcuts(): void {
    const dialog = $<HTMLDialogElement>('shortcuts-dialog');
    const groups = new Map<string, { label: string; shortcut: string }[]>();
    for (const action of listActions()) {
        if (!action.shortcut) continue;
        const group = action.group ?? 'General';
        const rows = groups.get(group) ?? [];
        rows.push({ label: action.label, shortcut: action.shortcut });
        groups.set(group, rows);
    }
    const body = $('shortcuts-body');
    body.replaceChildren(
        ...[...groups].flatMap(([group, rows]) => {
            const heading = document.createElement('div');
            heading.className = 'shortcut-group';
            heading.textContent = group;
            return [
                heading,
                ...rows.map((row) => {
                    const el = document.createElement('div');
                    el.className = 'shortcut-row';
                    const label = document.createElement('span');
                    label.textContent = row.label;
                    const kbd = document.createElement('kbd');
                    kbd.textContent = row.shortcut;
                    el.append(label, kbd);
                    return el;
                }),
            ];
        })
    );
    if (!dialog.open) dialog.showModal();
}
