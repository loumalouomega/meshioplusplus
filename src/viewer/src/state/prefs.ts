/**
 * Display preferences that survive a reload: colour map, scale, bands, draw
 * mode, opacity, edges, lighting and projection.
 *
 * Pure. What is stored is untrusted text (another version of the page, a
 * hand-edited value), so {@link parsePrefs} keeps only known controls with the
 * right type; applying a value to a `<select>` still checks the option exists.
 * Mesh-specific choices (the array to colour by, selections, clip) are
 * deliberately not here.
 */
export interface PrefField {
    id: string;
    kind: 'select' | 'range' | 'check';
}

/** Order matters: the colour-blind filter changes which colour maps exist. */
export const PREF_FIELDS: readonly PrefField[] = [
    { id: 'colormap-safe', kind: 'check' },
    { id: 'colormap', kind: 'select' },
    { id: 'color-scale', kind: 'select' },
    { id: 'color-bands', kind: 'select' },
    { id: 'representation', kind: 'select' },
    { id: 'opacity', kind: 'range' },
    { id: 'edges', kind: 'check' },
    { id: 'light-ambient', kind: 'range' },
    { id: 'light-diffuse', kind: 'range' },
    { id: 'light-specular', kind: 'range' },
    { id: 'light-cull', kind: 'check' },
    { id: 'view-projection', kind: 'select' },
];

export type Prefs = Record<string, string | boolean>;

export function parsePrefs(raw: unknown): Prefs {
    if (typeof raw !== 'object' || raw === null || Array.isArray(raw)) return {};
    const source = raw as Record<string, unknown>;
    const out: Prefs = {};
    for (const field of PREF_FIELDS) {
        const v = source[field.id];
        if (field.kind === 'check' && typeof v === 'boolean') out[field.id] = v;
        else if (field.kind === 'select' && typeof v === 'string' && v.length <= 60) out[field.id] = v;
        else if (field.kind === 'range' && typeof v === 'string' && v !== '' && Number.isFinite(Number(v))) out[field.id] = v;
    }
    return out;
}
