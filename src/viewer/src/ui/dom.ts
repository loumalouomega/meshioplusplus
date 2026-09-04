/**
 * DOM helpers.
 *
 * `$` throws with the id it was looking for rather than returning null. A
 * missing element is always a bug in `index.html`, and the alternative — an
 * optional chain that silently does nothing — once cost the sibling MDPA
 * viewer a whole webview that failed at module scope with no message.
 */

export function $<T extends HTMLElement = HTMLElement>(id: string): T {
    const el = document.getElementById(id);
    if (!el) throw new Error(`viewer: #${id} is missing from the page`);
    return el as T;
}

/**
 * Like {@link $}, for an SVG element. `getElementById` is typed
 * `HTMLElement | null`, which an `<svg>` is not, so this goes through
 * `querySelector<T>` -- properly typed for SVG, and no cast.
 */
export function svg$<T extends SVGElement = SVGSVGElement>(id: string): T {
    const el = document.querySelector<T>(`#${id}`);
    if (!el) throw new Error(`viewer: #${id} is missing from the page`);
    return el;
}

/** Like {@link $}, but for elements a given build may legitimately omit. */
export function maybe<T extends HTMLElement = HTMLElement>(id: string): T | null {
    return document.getElementById(id) as T | null;
}

export function show(el: HTMLElement, visible = true): void {
    el.hidden = !visible;
}

/** Replace a `<select>`'s options, preserving `selected` if it still exists. */
export function setOptions(
    select: HTMLSelectElement,
    options: { value: string; label: string }[],
    selected?: string
): void {
    const keep = selected ?? select.value;
    select.replaceChildren(
        ...options.map((o) => new Option(o.label, o.value, false, o.value === keep))
    );
    if (options.some((o) => o.value === keep)) select.value = keep;
}

/** A button with an inline SVG icon and a text label. */
export function iconButton(
    icon: string,
    label: string,
    title: string,
    onClick: () => void
): HTMLButtonElement {
    const button = document.createElement('button');
    button.type = 'button';
    button.title = title;
    button.innerHTML = `<span class="toolbar-icon">${icon}</span>`;
    if (label) button.append(document.createTextNode(label));
    button.addEventListener('click', onClick);
    return button;
}
