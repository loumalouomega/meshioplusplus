/**
 * The dataset page's typed test hook and session-state types — the sibling
 * of `src/types.ts`'s `ViewerState`, with the same discipline: Playwright
 * specs assert on `window.__datasetState` (never pixels), the specs import
 * this type, and a new field with un-updated tests is a compile error.
 *
 * UI-session state (selection, dirtiness, scan cache) lives HERE and never
 * in the manifest document itself — the document's schema is strict, and
 * stashing UI state in it would make Python reject the file.
 */

export interface EntryScan {
    steps: number;
    numNan: number;
    numInf: number;
}

export interface DatasetState {
    status: 'booting' | 'idle' | 'loading' | 'ready' | 'error';
    /** How the workspace was opened: File System Access, or the
     * webkitdirectory + download fallback. Null before a pick. */
    backendMode: 'fsa' | 'fallback' | null;
    manifestName: string | null;
    entryIds: string[];
    selected: string | null;
    /** The staged entry's resolved plan length and current scrubber step. */
    planLength: number;
    step: number;
    /** Of the current preview render. */
    numPoints: number;
    numCells: number;
    summary: { name: string; numNan: number; numInf: number }[] | null;
    dirty: boolean;
    /** The last save — `text` is the exact serialization, the byte-parity
     * test seam (asserting on it needs no filesystem access). */
    lastSave: { mode: 'in-place' | 'download'; text: string } | null;
    error: string | null;
}
