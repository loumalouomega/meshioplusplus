/**
 * Colour scheme: auto (follow the OS), dark, light, or "scientific" (a white
 * scene for figures).
 *
 * The chrome follows `data-theme` on <html> (tokens.css). The 3D scene's
 * background and default body colour follow the same choice, with the values
 * the Kratos MDPA Preview extension uses (webview/themes.ts), so the viewer
 * and the VS Code preview agree.
 *
 * Everything that needs no DOM is a plain function, unit-tested in Node.
 */
import type { Vector3 } from '../types';

export type ThemeMode = 'auto' | 'dark' | 'light' | 'scientific';
export type ResolvedTheme = 'dark' | 'light' | 'scientific';

export const THEME_MODES: readonly ThemeMode[] = ['auto', 'dark', 'light', 'scientific'];

export interface SceneTheme {
    background: Vector3;
    /** Surface colour when nothing is coloured by data. */
    body: Vector3;
    /** Categorical palette (groups, selection sets), 8 entries. */
    palette: Vector3[];
}

const DARK_PALETTE: Vector3[] = [
    [0.26, 0.59, 0.98],
    [0.96, 0.62, 0.1],
    [0.2, 0.73, 0.4],
    [0.91, 0.3, 0.24],
    [0.61, 0.35, 0.71],
    [0.1, 0.74, 0.74],
    [0.85, 0.65, 0.13],
    [0.55, 0.55, 0.6],
];

const LIGHT_PALETTE: Vector3[] = [
    [0.15, 0.45, 0.8],
    [0.8, 0.45, 0.05],
    [0.1, 0.55, 0.25],
    [0.75, 0.15, 0.1],
    [0.45, 0.2, 0.58],
    [0.05, 0.55, 0.55],
    [0.65, 0.48, 0.05],
    [0.38, 0.38, 0.42],
];

const SCIENTIFIC_PALETTE: Vector3[] = [
    [0.4, 0.76, 0.65],
    [0.99, 0.55, 0.38],
    [0.55, 0.63, 0.8],
    [0.91, 0.54, 0.76],
    [0.65, 0.85, 0.33],
    [1.0, 0.85, 0.18],
    [0.9, 0.77, 0.58],
    [0.7, 0.7, 0.7],
];

export const SCENE_THEMES: Record<ResolvedTheme, SceneTheme> = {
    dark: { background: [0.118, 0.118, 0.118], body: [0.78, 0.82, 0.88], palette: DARK_PALETTE },
    light: { background: [0.941, 0.941, 0.941], body: [0.58, 0.64, 0.74], palette: LIGHT_PALETTE },
    scientific: { background: [1, 1, 1], body: [0.62, 0.68, 0.78], palette: SCIENTIFIC_PALETTE },
};

export function isThemeMode(value: unknown): value is ThemeMode {
    return typeof value === 'string' && (THEME_MODES as readonly string[]).includes(value);
}

/** What "auto" and the explicit modes resolve to, given the OS preference. */
export function resolveTheme(mode: ThemeMode, prefersDark: boolean): ResolvedTheme {
    if (mode === 'auto') return prefersDark ? 'dark' : 'light';
    return mode;
}

/** `rgb(r, g, b)` for a 0..1 colour, for swatches. */
export function cssColor(color: Vector3): string {
    const c = (v: number) => Math.round(Math.min(1, Math.max(0, v)) * 255);
    return `rgb(${c(color[0])}, ${c(color[1])}, ${c(color[2])})`;
}

// --- DOM half -------------------------------------------------------------- //

const STORAGE_KEY = 'meshioplusplus.viewer.theme';

export function loadThemeMode(): ThemeMode {
    try {
        const stored = localStorage.getItem(STORAGE_KEY);
        if (isThemeMode(stored)) return stored;
    } catch {
        // storage can be unavailable (privacy mode, file:// in some browsers)
    }
    return 'auto';
}

function saveThemeMode(mode: ThemeMode): void {
    try {
        localStorage.setItem(STORAGE_KEY, mode);
    } catch {
        // not worth failing over
    }
}

export interface ThemeController {
    readonly mode: ThemeMode;
    readonly resolved: ResolvedTheme;
    setMode(mode: ThemeMode): void;
    /** Called once now and again whenever the resolved scheme changes. */
    onChange(listener: (resolved: ResolvedTheme, scene: SceneTheme) => void): void;
}

export function initTheme(): ThemeController {
    const media = window.matchMedia('(prefers-color-scheme: dark)');
    let mode = loadThemeMode();
    let resolved = resolveTheme(mode, media.matches);
    const listeners: Array<(r: ResolvedTheme, s: SceneTheme) => void> = [];

    const apply = () => {
        const next = resolveTheme(mode, media.matches);
        resolved = next;
        document.documentElement.dataset.theme = next;
        document.documentElement.dataset.themeMode = mode;
        for (const listener of listeners) listener(next, SCENE_THEMES[next]);
    };
    media.addEventListener('change', () => {
        if (mode === 'auto') apply();
    });
    document.documentElement.dataset.theme = resolved;
    document.documentElement.dataset.themeMode = mode;

    return {
        get mode() {
            return mode;
        },
        get resolved() {
            return resolved;
        },
        setMode(next) {
            if (next === mode) return;
            mode = next;
            saveThemeMode(next);
            apply();
        },
        onChange(listener) {
            listeners.push(listener);
            listener(resolved, SCENE_THEMES[resolved]);
        },
    };
}
