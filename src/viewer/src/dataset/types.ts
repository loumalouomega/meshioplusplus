/**
 * The dataset page's state hook (`window.__datasetState`), typed like the
 * viewer's `ViewerState` so the Playwright spec's expectations are checked at
 * compile time.
 *
 * UI-session state lives here and NEVER in the manifest document: the
 * document's schema is strict (an unknown key is an error), so stashing
 * `selected`/`step`/scan results in it would make Python reject the file.
 */

import type { CheckpointInfo, JobSummary, MetricRow } from './jobs';
import type { NotificationState } from './notify';

export interface EntryScan {
    steps: number;
    /** NaN/Inf over the entry's DATA arrays only — `quality:*` arrays are
     * excluded, since a quality metric's NaN means "N/A for this cell type"
     * by design, not a bad value. */
    numNan: number;
    numInf: number;
    /** The quality lane: cells with `quality:inverted` / `quality:degenerate`
     * set, and the worst `quality:scaled_jacobian` (null when the metric was
     * unavailable). */
    numInverted: number;
    numDegenerate: number;
    minScaledJacobian: number | null;
    /** Every data array present at step 0, as `<location>:<name>` — the
     * "fields missing across entries" input of the health summary. Empty for
     * an entry that could not be read (`steps === 0`). */
    arrays: string[];
    /** Steps in the entry's paired `Target`, or null when it has none (which
     * means self-supervised, not broken). */
    targetSteps?: number | null;
    /** Why the Source/Target pair does not line up, when it does not. A pair
     * whose sides disagree trains on mismatched steps — a silent wrong answer
     * rather than a crash — so it counts as a bad entry. */
    pairingError?: string | null;
}

/** One split's share of a manifest; `''` is "unassigned". */
export interface SplitBalance {
    split: string;
    count: number;
    fraction: number;
}

/**
 * A manifest-level health summary aggregated from per-entry scans. The
 * browser (WASM) and, later, the companion server are two producers of the
 * same shape — `producer` says which filled it.
 */
export interface ManifestHealth {
    producer: 'browser' | 'server';
    /** Entries with a scan / entries in the manifest. */
    scanned: number;
    total: number;
    numNan: number;
    numInf: number;
    numInverted: number;
    numDegenerate: number;
    minScaledJacobian: number | null;
    splitBalance: SplitBalance[];
    /** entry id -> data arrays (as `<location>:<name>`) that every OTHER
     * scanned entry carries and this one lacks. Only entries missing
     * something are listed. */
    fieldsMissing: Record<string, string[]>;
    /** Entries with any bad signal (NaN/Inf, inverted/degenerate cells, a
     * negative worst scaled Jacobian, or unreadable). */
    badEntries: string[];
}

/** One manifest in the overview grid. */
export interface ManifestCard {
    /** Workspace-relative path (root-level `*.json`). */
    path: string;
    name: string | null;
    description: string | null;
    numEntries: number;
    /** `{split: count}` with unassigned entries under `''`. */
    splits: Record<string, number>;
    /** Sorted unions over the entries. */
    tags: string[];
    groups: string[];
    /** Epoch milliseconds; null when the backend cannot say. */
    lastModified: number | null;
    /** A strict-parse failure still yields a card, flagged with the message. */
    parseError: string | null;
    health: ManifestHealth | null;
    /** A data URL captured from a preview render, session-only. */
    thumbnail: string | null;
    /** SHA-256 of the file text (null when `crypto.subtle` is unavailable). */
    sha256: string | null;
    /** Unsaved edits are pending in the drill-down for this manifest. */
    dirty: boolean;
    /** The same manifest on the companion process (bound by content hash),
     * or null when no connected server holds it. */
    serverPath: string | null;
    /** Known only to the server — not in the picked directory, so it cannot
     * be opened here; health still works through the server. */
    serverOnly: boolean;
}

/** The companion-process connection (`api.ts`), null until first tried. */
export interface ServerState {
    url: string;
    connected: boolean;
    version: string | null;
    tools: string[];
    root: string | null;
    runsDir: string | null;
    error: string | null;
}

export type DatasetView = 'overview' | 'manifest';

/** A prediction fetched from the companion process and rendered. */
export interface PredictionPreview {
    jobId: string;
    entryId: string;
    outputPath: string;
    rmse: number | null;
    maxError: number | null;
}

/** The run the page is following (`train.ts` / `jobs.ts`). */
export interface ActiveJob {
    jobId: string;
    status: string;
    epoch: number;
    epochs: number;
    metrics: MetricRow[];
    checkpoints: CheckpointInfo[];
    bestCheckpoint: string | null;
    /** The last polling error, cleared by the next successful poll. */
    error: string | null;
}

/** The counts of the last manifest diff shown. */
export interface DiffSummary {
    a: string;
    b: string;
    headerChanged: number;
    added: number;
    removed: number;
    changed: number;
}

export interface DatasetState {
    status: 'booting' | 'idle' | 'loading' | 'ready' | 'error';
    backendMode: 'fsa' | 'fallback' | null;
    /** A persisted directory handle was found; "Reopen" is offered. */
    restoreAvailable: boolean;
    /** Which depth of the page is showing: the overview grid of every
     * manifest in the workspace, or one manifest's curation view. */
    view: DatasetView;
    manifests: ManifestCard[];
    server: ServerState | null;
    /** Training runs the companion process knows (newest first). */
    jobs: JobSummary[];
    activeJob: ActiveJob | null;
    /** The runs selected for comparison, in the order picked. */
    compare: string[];
    /** The last prediction rendered in the viewer, or null. */
    prediction: PredictionPreview | null;
    notifications: NotificationState;
    manifestName: string | null;
    entryIds: string[];
    selected: string | null;
    planLength: number;
    step: number;
    numPoints: number;
    numCells: number;
    summary:
        | { name: string; min: number; mean: number; numNan: number; numInf: number }[]
        | null;
    scans: Record<string, EntryScan>;
    diff: DiffSummary | null;
    dirty: boolean;
    lastSave: { mode: 'in-place' | 'download'; text: string } | null;
    error: string | null;
}
