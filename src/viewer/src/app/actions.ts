/**
 * The action registry.
 *
 * Every user-visible command (a menu entry, a toolbar button, a keyboard
 * shortcut, a command-palette row) is one action with an id, a label and an
 * optional shortcut, so the same thing is reachable the same way from each and
 * the shortcut sheet can be generated instead of written twice.
 */
export interface Action {
    id: string;
    label: string;
    /** Display form, e.g. `Ctrl+O` or `1`. Matching is done by {@link matchShortcut}. */
    shortcut?: string;
    /** Group heading in the shortcut sheet and the palette. */
    group?: string;
    /** Return false to hide/disable the action in the current state. */
    enabled?: () => boolean;
    run: () => void;
}

const actions = new Map<string, Action>();

export function registerAction(action: Action): void {
    actions.set(action.id, action);
}

export function getAction(id: string): Action | undefined {
    return actions.get(id);
}

export function listActions(): Action[] {
    return [...actions.values()];
}

export function runAction(id: string): boolean {
    const action = actions.get(id);
    if (!action || (action.enabled && !action.enabled())) return false;
    action.run();
    return true;
}

export interface KeyLike {
    key: string;
    ctrlKey: boolean;
    metaKey: boolean;
    altKey: boolean;
    shiftKey: boolean;
}

/**
 * Whether a key event matches a display shortcut such as `Ctrl+Shift+Z`, `Esc`,
 * `1` or `/`. `Ctrl` also accepts the Command key so one spelling serves both.
 */
export function matchShortcut(shortcut: string, e: KeyLike): boolean {
    const parts = shortcut.split('+').map((p) => p.trim());
    const key = parts.pop() ?? '';
    const want = new Set(parts.map((p) => p.toLowerCase()));
    const ctrl = e.ctrlKey || e.metaKey;
    if (want.has('ctrl') !== ctrl) return false;
    if (want.has('alt') !== e.altKey) return false;
    // Letters and digits must match Shift exactly (so `I` and `Shift+I` are
    // different shortcuts). Punctuation such as `?` is typed with Shift, so
    // Shift only has to be present when the shortcut names it.
    if (/^[a-z0-9]$/i.test(key)) {
        if (want.has('shift') !== e.shiftKey) return false;
    } else if (want.has('shift') && !e.shiftKey) {
        return false;
    }
    const k = e.key === 'Escape' ? 'Esc' : e.key;
    return k.length === 1 ? k.toLowerCase() === key.toLowerCase() : k === key;
}

/** The first enabled action whose shortcut matches, if any. */
export function findByShortcut(e: KeyLike): Action | undefined {
    for (const action of actions.values()) {
        if (!action.shortcut) continue;
        if (action.enabled && !action.enabled()) continue;
        if (matchShortcut(action.shortcut, e)) return action;
    }
    return undefined;
}
