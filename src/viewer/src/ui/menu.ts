/**
 * Menubar dropdowns.
 *
 * One mechanism for every menu: a trigger `[data-menu="<name>"]` opens the
 * panel `#menu-<name>`. One menu is open at a time; an outside press or Escape
 * closes it (Escape returns focus to the trigger); hovering a sibling trigger
 * while one is open switches menus; the arrow keys, Home and End move between
 * items. Items carry `data-action="<action id>"` and run through the action
 * registry.
 */
import { runAction } from '../app/actions';

export function initMenus(root: ParentNode = document): void {
    const triggers = [...root.querySelectorAll<HTMLButtonElement>('[data-menu]')];
    const panelOf = (t: HTMLElement) =>
        document.getElementById(`menu-${t.dataset.menu ?? ''}`) as HTMLElement | null;

    const closeAll = (focusTrigger = false) => {
        for (const t of triggers) {
            const panel = panelOf(t);
            if (!panel || panel.hidden) continue;
            panel.hidden = true;
            t.setAttribute('aria-expanded', 'false');
            if (focusTrigger) t.focus();
        }
    };

    const open = (t: HTMLButtonElement) => {
        const panel = panelOf(t);
        if (!panel) return;
        closeAll();
        const rect = t.getBoundingClientRect();
        panel.style.left = `${Math.round(rect.left)}px`;
        panel.style.top = `${Math.round(rect.bottom + 2)}px`;
        panel.hidden = false;
        t.setAttribute('aria-expanded', 'true');
    };

    const isOpen = () => triggers.some((t) => panelOf(t)?.hidden === false);

    for (const t of triggers) {
        t.setAttribute('aria-haspopup', 'menu');
        t.setAttribute('aria-expanded', 'false');
        t.addEventListener('click', (e) => {
            e.stopPropagation();
            const panel = panelOf(t);
            if (panel && !panel.hidden) closeAll();
            else open(t);
        });
        t.addEventListener('pointerenter', () => {
            if (isOpen() && panelOf(t)?.hidden) open(t);
        });

        const panel = panelOf(t);
        panel?.addEventListener('click', (e) => {
            const item = (e.target as HTMLElement).closest<HTMLElement>('[data-action]');
            if (!item || item.hasAttribute('disabled')) return;
            const keepOpen = item.getAttribute('role') === 'menuitemcheckbox';
            if (!keepOpen) closeAll();
            runAction(item.dataset.action ?? '');
        });
        panel?.addEventListener('keydown', (e) => {
            const items = [
                ...panel.querySelectorAll<HTMLElement>('[role^="menuitem"]:not([disabled])'),
            ];
            const index = items.indexOf(document.activeElement as HTMLElement);
            let next = -1;
            if (e.key === 'ArrowDown') next = (index + 1) % items.length;
            else if (e.key === 'ArrowUp') next = (index - 1 + items.length) % items.length;
            else if (e.key === 'Home') next = 0;
            else if (e.key === 'End') next = items.length - 1;
            if (next >= 0) {
                e.preventDefault();
                items[next]?.focus();
            }
        });
        t.addEventListener('keydown', (e) => {
            if (e.key === 'ArrowDown' || e.key === 'Enter' || e.key === ' ') {
                e.preventDefault();
                open(t);
                panelOf(t)?.querySelector<HTMLElement>('[role^="menuitem"]')?.focus();
            }
        });
    }

    // Capture phase: the closing click must not reach the live 3D canvas.
    document.addEventListener(
        'pointerdown',
        (e) => {
            if (!isOpen()) return;
            const target = e.target as HTMLElement;
            if (target.closest('[data-menu]') || target.closest('.menu')) return;
            closeAll();
        },
        true
    );
    document.addEventListener('keydown', (e) => {
        if (e.key === 'Escape' && isOpen()) {
            e.preventDefault();
            closeAll(true);
        }
    });
}

/** Mirror a boolean into a checkbox menu item's `aria-checked`. */
export function setMenuChecked(action: string, checked: boolean): void {
    for (const item of document.querySelectorAll<HTMLElement>(`[data-action="${action}"]`)) {
        item.setAttribute('aria-checked', String(checked));
    }
}
