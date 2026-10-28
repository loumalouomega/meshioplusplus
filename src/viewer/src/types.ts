/**
 * Shared types for the viewer.
 *
 * `ViewerState` is the documented test hook (`window.__viewerState`). Owning it
 * here rather than building it ad hoc in `main.ts` is the main reason this app
 * is TypeScript: the Playwright specs import this type, so adding a field and
 * forgetting to update a test becomes a compile error instead of a silent pass.
 */

export type Vector3 = [number, number, number];

/** Where a data array lives. */
export type DataLocation = 'point' | 'cell';

/** One colourable choice: an array, or one component of it. */
export interface ArrayEntry {
    location: DataLocation;
    name: string;
    /** Component index, or -1 for the magnitude of a multi-component array. */
    component: number;
    /** What the colour-by menu shows. */
    label: string;
}

/** Result of colouring by an array. */
export interface ScalarRange {
    min: number;
    max: number;
    /** Non-finite entries, excluded from the range and shown in the legend. */
    nanCount: number;
}

/** One data array's value at a picked element, all components. */
export interface PickValue {
    name: string;
    components: number[];
}

/** What clicking a cell reports. */
export interface PickInfo {
    /** Index into the *rendered surface*, not the original mesh. */
    cellId: number;
    cellType: string;
    pointId: number;
    position: Vector3;
    /** Where the ray actually hit the surface (not snapped to a vertex). */
    hit: Vector3;
    /**
     * The originating cell, from `surface:parent_cell` — meaningful only when
     * no operations have been applied, since after one it names a cell of the
     * post-operation mesh. Null when the array is absent.
     */
    parentCell: number | null;
    cellValues: PickValue[];
    pointValues: PickValue[];
}

/** The selection tool's state, for the test hook. */
export interface SelectionState {
    tool: 'off' | 'click' | 'box' | 'lasso';
    active: number | null;
    sets: { id: number; name: string; count: number }[];
    /** How many surface cells are currently hidden. */
    hidden: number;
}

/** The groups panel's state, for the test hook. */
export interface GroupsState {
    status: 'none' | 'loading' | 'ready' | 'unavailable' | 'failed';
    groups: { id: number; name: string; list: 'regions' | 'blocks'; cells: number; visible: boolean; user: boolean }[];
}

/** The viewing tools' state, for the test hook. */
export interface ViewState {
    parallel: boolean;
    /** Unit vector from the focal point toward the eye, rounded to 3 places. */
    direction: Vector3;
    logScale: boolean;
    bands: number;
    clip: boolean;
    /** Largest displacement drawn by the deformed shape (0 = none). */
    deformed: number;
    measuring: boolean;
    /** Distance between the two measure marks, or null. */
    measured: number | null;
    bookmarks: string[];
}

/** The table and statistics panels' state, for the test hook. */
export interface AnalysisState {
    table: { open: boolean; kind: 'cells' | 'points'; rows: number; sortedBy: string | null; direction: 'asc' | 'desc' };
    histogram: { array: string | null; counts: number[] };
}

export type ViewerStatus = 'booting' | 'idle' | 'loading' | 'ready' | 'error';

/**
 * The Playwright hook. Deliberately a documented part of the app: asserting on
 * this instead of on pixels keeps the end-to-end tests GPU-independent.
 */
export interface ViewerState {
    status: ViewerStatus;
    numPoints: number;
    numCells: number;
    /** Labels from {@link ArrayEntry}, in menu order. */
    arrays: string[];
    name: string;
    /** The most recent click-to-inspect result, or null. */
    picked: PickInfo | null;
    selection: SelectionState;
    groups: GroupsState;
    view: ViewState;
    analysis: AnalysisState;
}
