/**
 * The TypeScript twin of `src/python/meshioplusplus/_dataset.py`'s document
 * model — parse, validate, curate and serialize a DatasetManifest JSON.
 *
 * Two rules carried over verbatim from the Python side (doc/datasets.md):
 *
 * - **Strict parsing.** An unknown key anywhere is an error naming the
 *   offender, exactly like `_check_keys`. The UI refuses to load a document
 *   it does not fully understand rather than silently destroying the parts
 *   it doesn't — there is nothing legal to preserve, since Python itself
 *   would reject the same file.
 * - **Serialization parity.** `stringifyManifest` reproduces
 *   `DatasetManifest.save`'s bytes: insertion-ordered keys, `indent=2`,
 *   optional keys omitted when empty, non-ASCII escaped to `\uXXXX`
 *   (Python's `ensure_ascii` default) and a trailing newline — so hand
 *   edits, CLI edits and UI edits diff cleanly against one another. The one
 *   documented gap: JavaScript has no int/float distinction, so a
 *   Python-written `0.0` normalizes to `0` on the first UI save.
 *
 * Pure data — no DOM, no worker, no wasm — so the whole model is exercised
 * by `tests/viewer/unit/manifest.test.mjs`.
 */

export type TimeFrom = 'auto' | 'file' | 'filename' | 'index';

export interface SourceSpec {
    Pattern?: string;
    Path?: string;
    Paths?: string[];
    Format?: string;
    Times?: number[];
    TimeFrom?: TimeFrom;
    Sort?: boolean;
}

export interface ManifestEntry {
    id: string;
    source: SourceSpec;
    /** The paired coarse/fine series, when there is one. Null means
     *  self-supervised: one mesh supplies both sides. */
    target: SourceSpec | null;
    split: string | null;
    tags: string[];
    group: string | null;
    notes: string | null;
    metadata: Record<string, unknown>;
}

export interface Manifest {
    name: string | null;
    description: string | null;
    metadata: Record<string, unknown>;
    entries: ManifestEntry[];
}

const TOP_KEYS = ['Version', 'Name', 'Description', 'Metadata', 'Entries'];
const ENTRY_KEYS = ['Id', 'Source', 'Target', 'Split', 'Tags', 'Group', 'Notes', 'Metadata'];
const SOURCE_KEYS = ['Pattern', 'Path', 'Paths', 'Format', 'Times', 'TimeFrom', 'Sort'];
const TIME_FROM: TimeFrom[] = ['auto', 'file', 'filename', 'index'];

function err(message: string): Error {
    return new Error(`meshio++: dataset: ${message}`);
}

function checkKeys(obj: Record<string, unknown>, where: string, allowed: string[]): void {
    for (const key of Object.keys(obj)) {
        if (!allowed.includes(key)) {
            throw err(`unknown key '${key}' in ${where} (known: ${allowed.join(', ')})`);
        }
    }
}

function isObject(value: unknown): value is Record<string, unknown> {
    return typeof value === 'object' && value !== null && !Array.isArray(value);
}

function validateSource(raw: unknown, where: string, block = 'Source'): SourceSpec {
    where = `${where}.${block}`;
    if (!isObject(raw)) throw err(`${where} must be an object`);
    checkKeys(raw, where, SOURCE_KEYS);
    const kinds = (['Pattern', 'Path', 'Paths'] as const).filter((k) => k in raw);
    if (kinds.length !== 1) {
        throw err(
            `${where} needs exactly one of Pattern, Path or Paths ` +
                `(got ${kinds.length ? kinds.join(', ') : 'none'})`,
        );
    }
    const out: SourceSpec = {};
    const kind = kinds[0];
    if (kind === 'Paths') {
        const paths = raw.Paths;
        if (!Array.isArray(paths) || paths.length === 0) {
            throw err(`${where}.Paths must be a non-empty array of paths`);
        }
        if (!paths.every((p) => typeof p === 'string' && p.length > 0)) {
            throw err(`${where}.Paths entries must be non-empty strings`);
        }
        out.Paths = [...paths];
    } else {
        const value = raw[kind];
        if (typeof value !== 'string' || !value) {
            throw err(`${where}.${kind} must be a non-empty string`);
        }
        out[kind] = value;
    }
    if ('Format' in raw) {
        if (typeof raw.Format !== 'string' || !raw.Format) {
            throw err(`${where}.Format must be a non-empty string`);
        }
        out.Format = raw.Format;
    }
    if ('Times' in raw) {
        const times = raw.Times;
        if (!Array.isArray(times) || !times.every((t) => typeof t === 'number')) {
            throw err(`${where}.Times must be an array of numbers`);
        }
        out.Times = [...times];
    }
    if ('TimeFrom' in raw) {
        const tf = raw.TimeFrom;
        if (typeof tf !== 'string' || !(TIME_FROM as string[]).includes(tf)) {
            throw err(
                `${where}.TimeFrom must be one of ` +
                    TIME_FROM.map((v) => `'${v}'`).join(', '),
            );
        }
        out.TimeFrom = tf as TimeFrom;
    }
    if ('Sort' in raw) {
        if (typeof raw.Sort !== 'boolean') {
            throw err(`${where}.Sort must be true or false`);
        }
        if (kind !== 'Paths') {
            throw err(`${where}.Sort applies only with Paths`);
        }
        out.Sort = raw.Sort;
    }
    return out;
}

function validateMetadata(value: unknown, where: string): Record<string, unknown> {
    if (value === undefined || value === null) return {};
    if (!isObject(value)) throw err(`${where} must be an object`);
    return { ...value };
}

function parseEntry(raw: unknown, where: string): ManifestEntry {
    if (!isObject(raw)) throw err(`${where} must be an object`);
    checkKeys(raw, where, ENTRY_KEYS);
    const id = raw.Id;
    if (typeof id !== 'string' || !id) {
        throw err(`${where}.Id is required and must be a non-empty string`);
    }
    if (!('Source' in raw)) throw err(`${where}.Source is required`);
    const source = validateSource(raw.Source, where);
    const target =
        raw.Target === undefined || raw.Target === null
            ? null
            : validateSource(raw.Target, where, 'Target');
    const split = raw.Split;
    if (split !== undefined && (typeof split !== 'string' || !split)) {
        throw err(`${where}.Split must be a non-empty string`);
    }
    const tags = raw.Tags ?? [];
    if (!Array.isArray(tags) || !tags.every((t) => typeof t === 'string' && t.length > 0)) {
        throw err(`${where}.Tags must be an array of non-empty strings`);
    }
    const group = raw.Group;
    if (group !== undefined && (typeof group !== 'string' || !group)) {
        throw err(`${where}.Group must be a non-empty string`);
    }
    const notes = raw.Notes;
    if (notes !== undefined && typeof notes !== 'string') {
        throw err(`${where}.Notes must be a string`);
    }
    return {
        id,
        source,
        target,
        split: typeof split === 'string' ? split : null,
        tags: [...tags],
        group: typeof group === 'string' ? group : null,
        notes: typeof notes === 'string' ? notes : null,
        metadata: validateMetadata(raw.Metadata, `${where}.Metadata`),
    };
}

/** Parse a manifest document (JSON text or an already-parsed object). */
export function parseManifest(input: string | unknown): Manifest {
    let doc: unknown;
    if (typeof input === 'string') {
        try {
            doc = JSON.parse(input);
        } catch (e) {
            throw err(`the manifest is not valid JSON: ${(e as Error).message}`);
        }
    } else {
        doc = input;
    }
    if (!isObject(doc)) throw err('the manifest document must be an object');
    checkKeys(doc, 'the manifest document', TOP_KEYS);
    const version = doc.Version ?? 1;
    if (version !== 1) {
        throw err(`unsupported Version ${JSON.stringify(version)} (this build knows 1)`);
    }
    const name = doc.Name;
    if (name !== undefined && typeof name !== 'string') throw err('Name must be a string');
    const description = doc.Description;
    if (description !== undefined && typeof description !== 'string') {
        throw err('Description must be a string');
    }
    const rawEntries = doc.Entries ?? [];
    if (!Array.isArray(rawEntries)) throw err('Entries must be an array');
    const entries = rawEntries.map((e, i) => parseEntry(e, `Entries[${i}]`));
    const seen = new Set<string>();
    for (const entry of entries) {
        if (seen.has(entry.id)) throw err(`duplicate entry id '${entry.id}'`);
        seen.add(entry.id);
    }
    return {
        name: typeof name === 'string' ? name : null,
        description: typeof description === 'string' ? description : null,
        metadata: validateMetadata(doc.Metadata, 'Metadata'),
        entries,
    };
}

function entryToDoc(entry: ManifestEntry): Record<string, unknown> {
    const out: Record<string, unknown> = { Id: entry.id, Source: { ...entry.source } };
    if (entry.target !== null) out.Target = { ...entry.target };
    if (entry.split !== null) out.Split = entry.split;
    if (entry.tags.length) out.Tags = [...entry.tags];
    if (entry.group !== null) out.Group = entry.group;
    if (entry.notes !== null) out.Notes = entry.notes;
    if (Object.keys(entry.metadata).length) out.Metadata = { ...entry.metadata };
    return out;
}

/** The document object, key order matching `DatasetManifest.to_dict`. */
export function manifestToDoc(manifest: Manifest): Record<string, unknown> {
    const out: Record<string, unknown> = { Version: 1 };
    if (manifest.name !== null) out.Name = manifest.name;
    if (manifest.description !== null) out.Description = manifest.description;
    if (Object.keys(manifest.metadata).length) out.Metadata = { ...manifest.metadata };
    out.Entries = manifest.entries.map(entryToDoc);
    return out;
}

/**
 * Serialize with `DatasetManifest.save`'s exact conventions: `indent=2`
 * (JSON.stringify's layout matches Python's `json.dump(indent=2)` element
 * for element), non-ASCII escaped like Python's `ensure_ascii=True`, and a
 * trailing newline.
 */
export function stringifyManifest(manifest: Manifest): string {
    const text = JSON.stringify(manifestToDoc(manifest), null, 2).replace(
        /[\u0080-\uffff]/g,
        (ch) => `\\u${ch.charCodeAt(0).toString(16).padStart(4, '0')}`,
    );
    return `${text}\n`;
}

/** Mirror `_derive_id`: source stem, glob chars removed, `-_.` trimmed. */
export function deriveId(source: SourceSpec, existing: Iterable<string>): string {
    const text = source.Paths ? source.Paths[0] : (source.Pattern ?? source.Path ?? '');
    const base = text.split('/').pop() ?? '';
    const dot = base.lastIndexOf('.');
    let stem = dot > 0 ? base.slice(0, dot) : base;
    stem = stem.replace(/[*?]/g, '').replace(/[-_.]+$/, '');
    if (!stem) throw err(`cannot derive an id from '${text}'; pass one explicitly`);
    const taken = new Set(existing);
    if (taken.has(stem)) {
        throw err(`the derived id '${stem}' already exists; pass a distinct id`);
    }
    return stem;
}

/** Curation helpers — every mutation replaces entries wholesale, the
 * frozen-dataclass discipline of the Python side. */
export function addEntry(
    manifest: Manifest,
    source: SourceSpec,
    opts: {
        id?: string;
        target?: SourceSpec | null;
        split?: string | null;
        tags?: string[];
        group?: string | null;
        notes?: string | null;
        metadata?: Record<string, unknown>;
    } = {},
): ManifestEntry {
    const validated = validateSource(source, 'add()');
    const ids = manifest.entries.map((e) => e.id);
    const id = opts.id ?? deriveId(validated, ids);
    if (!id) throw err('the entry id must be a non-empty string');
    if (ids.includes(id)) throw err(`an entry '${id}' already exists; pass a distinct id`);
    const entry: ManifestEntry = {
        id,
        source: validated,
        target: opts.target ? validateSource(opts.target, 'add()', 'Target') : null,
        split: opts.split ?? null,
        tags: [...(opts.tags ?? [])],
        group: opts.group ?? null,
        notes: opts.notes ?? null,
        metadata: { ...(opts.metadata ?? {}) },
    };
    manifest.entries.push(entry);
    return entry;
}

export function removeEntry(manifest: Manifest, id: string): void {
    const index = manifest.entries.findIndex((e) => e.id === id);
    if (index < 0) throw err(`no entry '${id}'`);
    manifest.entries.splice(index, 1);
}

export function getEntry(manifest: Manifest, id: string): ManifestEntry {
    const entry = manifest.entries.find((e) => e.id === id);
    if (!entry) throw err(`no entry '${id}'`);
    return entry;
}

function replaceEntry(manifest: Manifest, id: string, changes: Partial<ManifestEntry>): void {
    const index = manifest.entries.findIndex((e) => e.id === id);
    if (index < 0) throw err(`no entry '${id}'`);
    manifest.entries[index] = { ...manifest.entries[index], ...changes };
}

export function setSplit(manifest: Manifest, ids: string[] | 'all', split: string | null): void {
    const targets = ids === 'all' ? manifest.entries.map((e) => e.id) : ids;
    for (const id of targets) getEntry(manifest, id); // fail before mutating
    for (const id of targets) replaceEntry(manifest, id, { split });
}

/** Order-preserving tag add/remove, duplicates never accumulating —
 * `DatasetManifest.tag`'s semantics. */
export function tagEntries(
    manifest: Manifest,
    ids: string[] | 'all',
    changes: { add?: string[]; remove?: string[] },
): void {
    const add = changes.add ?? [];
    const remove = new Set(changes.remove ?? []);
    const targets = ids === 'all' ? manifest.entries.map((e) => e.id) : ids;
    for (const id of targets) getEntry(manifest, id);
    for (const id of targets) {
        const entry = getEntry(manifest, id);
        const tags = entry.tags.filter((t) => !remove.has(t));
        for (const t of add) if (!tags.includes(t) && !remove.has(t)) tags.push(t);
        replaceEntry(manifest, id, { tags });
    }
}

/** `annotate`'s merge semantics: null leaves a field unchanged, metadata
 * merges key-wise, dropMetadata removes keys. */
export function annotateEntry(
    manifest: Manifest,
    id: string,
    changes: {
        notes?: string | null;
        group?: string | null;
        metadata?: Record<string, unknown>;
        dropMetadata?: string[];
    },
): void {
    const entry = getEntry(manifest, id);
    const patch: Partial<ManifestEntry> = {};
    if (changes.notes !== undefined && changes.notes !== null) patch.notes = changes.notes;
    if (changes.group !== undefined && changes.group !== null) patch.group = changes.group;
    if (changes.metadata || changes.dropMetadata?.length) {
        const merged = { ...entry.metadata, ...(changes.metadata ?? {}) };
        for (const key of changes.dropMetadata ?? []) delete merged[key];
        patch.metadata = merged;
    }
    if (Object.keys(patch).length) replaceEntry(manifest, id, patch);
}

/** `entries(split=, tags=, group=)`'s filter semantics: exact split,
 * all-of tags, group path-segment prefix. */
export function filterEntries(
    manifest: Manifest,
    filter: { split?: string; tags?: string[]; group?: string },
): ManifestEntry[] {
    const wanted = new Set(filter.tags ?? []);
    return manifest.entries.filter((entry) => {
        if (filter.split !== undefined && entry.split !== filter.split) return false;
        if (wanted.size && ![...wanted].every((t) => entry.tags.includes(t))) return false;
        if (filter.group !== undefined) {
            const have = entry.group ?? '';
            if (have !== filter.group && !have.startsWith(`${filter.group}/`)) return false;
        }
        return true;
    });
}

/** `{split: count}` with unassigned entries under the empty string. */
export function splitCounts(manifest: Manifest): Map<string, number> {
    const counts = new Map<string, number>();
    for (const entry of manifest.entries) {
        const key = entry.split ?? '';
        counts.set(key, (counts.get(key) ?? 0) + 1);
    }
    return counts;
}

export function emptyManifest(): Manifest {
    return { name: null, description: null, metadata: {}, entries: [] };
}
