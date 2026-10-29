/** Save the display controls to localStorage and put them back on load. */
import { PREF_FIELDS, type Prefs, parsePrefs } from '../state/prefs';
import { $ } from '../ui/dom';

const KEY = 'meshioplusplus.viewer.prefs';

function load(): Prefs {
    try {
        const raw = localStorage.getItem(KEY);
        return raw ? parsePrefs(JSON.parse(raw) as unknown) : {};
    } catch {
        return {};
    }
}

function save(): void {
    const out: Prefs = {};
    for (const field of PREF_FIELDS) {
        const el = $<HTMLInputElement | HTMLSelectElement>(field.id);
        out[field.id] = field.kind === 'check' ? (el as HTMLInputElement).checked : el.value;
    }
    try {
        localStorage.setItem(KEY, JSON.stringify(out));
    } catch {
        // not worth failing over
    }
}

/** Call once every control is wired, so applying a value runs its listeners. */
export function initPrefs(): void {
    const prefs = load();
    for (const field of PREF_FIELDS) {
        const el = $<HTMLInputElement | HTMLSelectElement>(field.id);
        const v = prefs[field.id];
        if (v === undefined) continue;
        if (field.kind === 'check') {
            (el as HTMLInputElement).checked = v === true;
        } else if (el instanceof HTMLSelectElement) {
            if (![...el.options].some((o) => o.value === v)) continue;
            el.value = String(v);
        } else {
            el.value = String(v);
        }
        el.dispatchEvent(new Event(field.kind === 'range' ? 'input' : 'change', { bubbles: true }));
    }
    for (const field of PREF_FIELDS) {
        $(field.id).addEventListener(field.kind === 'range' ? 'input' : 'change', save);
    }
}
