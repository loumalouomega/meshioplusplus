/**
 * Collapsible sidebar sections.
 *
 * Turns each `<section class="sb-section">` whose first child is an `<h2>` into
 * a header (chevron button + title) over a body. Done at runtime so the HTML
 * stays plain and the ids the specs and the offline page rely on do not move.
 * The collapsed state survives reloads (localStorage, best effort).
 */
import { glyph } from './glyphs';

const STORAGE_KEY = 'meshioplusplus.viewer.collapsed';

function loadCollapsed(): Set<string> {
    try {
        const raw = localStorage.getItem(STORAGE_KEY);
        const parsed: unknown = raw ? JSON.parse(raw) : [];
        return new Set(Array.isArray(parsed) ? parsed.filter((x) => typeof x === 'string') : []);
    } catch {
        return new Set();
    }
}

function saveCollapsed(ids: Set<string>): void {
    try {
        localStorage.setItem(STORAGE_KEY, JSON.stringify([...ids]));
    } catch {
        // not worth failing over
    }
}

export function setCollapsed(section: HTMLElement, collapsed: boolean): void {
    section.classList.toggle('collapsed', collapsed);
    const button = section.querySelector<HTMLButtonElement>(':scope > .sb-header > .sb-chevron');
    button?.setAttribute('aria-expanded', String(!collapsed));
}

export function initSections(root: ParentNode = document): void {
    const collapsed = loadCollapsed();
    for (const section of root.querySelectorAll<HTMLElement>('section.sb-section')) {
        const title = section.querySelector<HTMLElement>(':scope > h2');
        if (!title) continue;

        const header = document.createElement('div');
        header.className = 'sb-header';
        const button = document.createElement('button');
        button.type = 'button';
        button.className = 'sb-chevron';
        button.innerHTML = glyph('chevronDown');
        const bodyId = `${section.id || 'section'}-body`;
        button.setAttribute('aria-controls', bodyId);
        button.setAttribute('aria-label', `Toggle ${title.textContent?.trim() ?? 'section'}`);

        const body = document.createElement('div');
        body.className = 'sb-body';
        body.id = bodyId;
        while (title.nextSibling) body.append(title.nextSibling);

        title.replaceWith(header);
        header.append(button, title);
        section.append(body);

        const toggle = () => {
            const next = !section.classList.contains('collapsed');
            setCollapsed(section, next);
            if (section.id) {
                if (next) collapsed.add(section.id);
                else collapsed.delete(section.id);
                saveCollapsed(collapsed);
            }
        };
        button.addEventListener('click', toggle);
        title.addEventListener('click', toggle);
        setCollapsed(section, !!section.id && collapsed.has(section.id));
    }
}
