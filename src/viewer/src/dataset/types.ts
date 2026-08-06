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
    /** NaN/Inf over the entry's DATA arrays only — `quality:*` arrays are
     * excluded, since a quality metric's NaN means "N/A for this cell type"
     * by design, not a bad value. */
    numNan: number;
    numInf: number;
    /** The quality lane: cells with `quality:inverted` set, and the worst
     * `quality:scaled_jacobian` (null when the metric was unavailable). */
    numInverted: number;
    minScaledJacobian: number | null;
}

export interface DatasetState {
    status: 'booting' | 'idle' | 'loading' | 'ready' | 'error';
    /** How the workspace was opened: File System Access, or the
     * webkitdirectory + download fallback. Null before a pick. */
    backendMode: 'fsa' | 'fallback' | null;
    /** A persisted directory handle exists — the "Reopen" button is shown
     * (restore itself runs behind its click; permission needs the gesture). */
    restoreAvailable: boolean;
    manifestName: string | null;
    entryIds: string[];
    selected: string | null;
    /** The staged entry's resolved plan length and current scrubber step. */
    planLength: number;
    step: number;
    /** Of the current preview render. */
    numPoints: number;
    numCells: number;
    summary:
        | { name: string; min: number; mean: number; numNan: number; numInf: number }[]
        | null;
    /** Per-entry scan results (filled by the Scan action). */
    scans: Record<string, EntryScan>;
    dirty: boolean;
    /** The last save — `text` is the exact serialization, the byte-parity
     * test seam (asserting on it needs no filesystem access). */
    lastSave: { mode: 'in-place' | 'download'; text: string } | null;
    error: string | null;
}
