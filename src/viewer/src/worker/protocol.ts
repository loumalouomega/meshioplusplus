/**
 * The worker contract, imported by **both** halves.
 *
 * Before the TypeScript migration this was a convention held together by
 * comments in two files; now a mismatch between what the client sends and what
 * the worker handles is a compile error.
 *
 * ## Why operations are a *pipeline*, not a sequence of mutations
 *
 * Every mesh operation in the WASM API takes and returns a JS `Mesh`, whose
 * flat representation cannot carry multi-component (vector/tensor) arrays. So
 * applying an operation through that API would silently destroy exactly the
 * data the display path goes out of its way to preserve.
 *
 * Instead the worker keeps the **original file bytes** staged for the life of
 * the session and replays the whole pipeline in C++ on every change, through
 * one `convertSurfaceOps` call. Three consequences, all good:
 *
 * - nothing is ever lost, because no mesh crosses the JS boundary;
 * - **undo is exact** — pop the last op and replay from pristine bytes, with
 *   no inverse operations and no snapshots;
 * - the display path and the post-operation display path are literally the
 *   same call, so they cannot drift.
 */

import type { Vector3 } from '../types';

/** One operation in the pipeline. */
export type OpSpec =
    | {
          op: 'quality';
      }
    | {
          op: 'clean';
          weld: boolean;
          atol: number;
          removeOrphans: boolean;
          dropDegenerate: boolean;
          dropDuplicateCells: boolean;
      }
    | {
          op: 'smooth';
          method: 'laplacian' | 'taubin';
          iterations: number;
          lambda: number;
          mu: number;
          fixBoundary: boolean;
      }
    | {
          op: 'refine';
          levels: number;
          /**
           * Refine only the cells satisfying a threshold on a scalar cell_data
           * array, resolving the hanging nodes that leaves with `closure`.
           * Empty `array` means "every cell" — the uniform refinement.
           */
          array: string;
          /** `compare`, not `op`: `op` is this union's own discriminant. */
          compare: '<' | '>';
          value: number;
          closure: 'redgreen' | 'propagate' | 'balanced';
          recordLevels: boolean;
      }
    | {
          op: 'partition';
          nparts: number;
          method: 'auto' | 'sfc' | 'kahip';
      }
    | {
          /**
           * The planar cross-section through the mesh (slice): the volume is
           * replaced by the surface where the plane intersects it, one
           * dimension lower.
           *
           * Stored in **world** coordinates rather than as a normalized slider
           * position, so the section stays where the user put it even if an
           * earlier operation changes the mesh's extent. `mode` is carried for
           * backward compatibility but ignored by a cross-section (there is no
           * "keep side").
           */
          op: 'section';
          point: Vector3;
          normal: Vector3;
          mode: 'all' | 'any';
      }
    | {
          /**
           * The level set of a scalar `point_data` field (isosurface) — the
           * data-driven sibling of `section`, and like it a surface one
           * dimension lower. `component` is negative for the row magnitude.
           */
          op: 'isosurface';
          array: string;
          isovalue: number;
          component: number;
      }
    | {
          /**
           * The gradient, divergence or curl of a `point_data` field. Unlike
           * every other op here it changes no geometry -- it only attaches one
           * array -- so the pipeline carries the mesh straight through and the
           * new field shows up in the colour-by menu.
           */
          op: 'gradient';
          array: string;
          operator: 'gradient' | 'divergence' | 'curl';
          method: 'green-gauss' | 'least-squares';
          location: 'point' | 'cell';
          output: string;
      }
    | {
          /**
           * A regular hexahedron grid around the mesh. The ONLY op here that
           * replaces the geometry rather than transforming it: the lattice is
           * what you see afterwards, not the mesh that went in. Undo restores
           * it, since the pipeline is replayed against the original bytes.
           */
          op: 'voxelize';
          /**
           * Cell counts per axis. The panel offers one number and sends it three
           * times: the spec goes to `convertSurfaceOps` verbatim, so it has to be
           * the wire shape rather than a friendlier one.
           */
          resolution: number[];
          fill: 'all' | 'surface' | 'inside';
      }
    | {
          /**
           * A signed distance field: a grid over the mesh's surface, filled.
           * Like `voxelize` it REPLACES the geometry -- what you see afterwards
           * is the grid, coloured by `sdf:distance`. `octree` refines only near
           * the surface, so its output is 1-irregular (it has hanging nodes).
           */
          op: 'computeSdf';
          structure: 'voxel' | 'octree';
          /** Voxel only; the panel sends one number three times (the wire shape). */
          resolution: number[];
          /** Octree only: the root lattice, then how many halving passes. */
          rootResolution: number;
          maxDepth: number;
      };

export type OpName = OpSpec['op'];

/** Per-operation counters the C++ side reports back. */
export interface OpStepReport {
    op: OpName;
    [counter: string]: number | string;
}

export interface OpReport {
    steps: OpStepReport[];
    /** Caveats authored in C++ — never guessed at in the UI. */
    warnings: string[];
}

/** A cheap summary of the file, from `readMetadata`. */
export interface MeshMeta {
    numPoints: number;
    numCells: number;
    format: string;
    cellBlocks: { type: string; numCells: number }[];
    pointDataNames: string[];
    cellDataNames: string[];
    bboxMin?: number[];
    bboxMax?: number[];
    /** Per-step time values of a multi-step file — the dataset page's step
     * scrubber's count; `readMetadata` always returns it. */
    timeValues?: number[];
}

// --- dataset-manager requests (dataset.html) ------------------------------ //
// Strictly additive: the viewer's own requests and its single staged
// `session` are untouched, so the existing Playwright suite is the
// regression net for everything below.

/** One resolved step of a staged entry (`sequenceEntries`' shape). */
export interface PlanEntry {
    path: string;
    step: number;
    time: number;
    timeSource: 'explicit' | 'file' | 'filename' | 'index';
}

/** The `Source` plan inputs of a manifest entry, workspace-relative
 * (doc/datasets.md); the worker prefixes them into its staging directory. */
export interface StagedSource {
    Pattern?: string;
    Path?: string;
    Paths?: string[];
    Format?: string;
    Times?: number[];
    TimeFrom?: 'auto' | 'file' | 'filename' | 'index';
    Sort?: boolean;
}

/** The per-array facts the summary table shows — the "bad case" signals. */
export interface ArraySummary {
    location: string;
    name: string;
    dtype: string;
    numValues: number;
    numComponents: number;
    min: number;
    max: number;
    mean: number;
    numNan: number;
    numInf: number;
}

export interface StageResult {
    plan: PlanEntry[];
    meta: MeshMeta;
}

export interface SummaryResult {
    arrays: ArraySummary[];
}

/**
 * `Omit` collapses a union to its common keys, which would erase every
 * request's own payload. Distributing over the union keeps them.
 */
export type DistributiveOmit<T, K extends PropertyKey> = T extends unknown
    ? Omit<T, K>
    : never;

export type Request =
    | { id: number; type: 'init' }
    | { id: number; type: 'open'; name: string; bytes: ArrayBuffer }
    | { id: number; type: 'apply'; ops: OpSpec[] }
    | { id: number; type: 'convert'; outFormat: string; ops: OpSpec[] }
    | { id: number; type: 'close' }
    | {
          id: number;
          type: 'stageEntry';
          files: { relPath: string; bytes: ArrayBuffer }[];
          source: StagedSource;
      }
    | { id: number; type: 'previewStep'; step: number }
    | { id: number; type: 'summaryStep'; step: number }
    | { id: number; type: 'evictEntry' };

export interface InitResult {
    formats: { readers: string[]; writers: string[] };
    backend: string;
}

export interface RenderResult {
    vtp: ArrayBuffer;
    meta: MeshMeta;
    format: string;
    report: OpReport;
    /** World-space bounds of the *original* mesh, for the section slider. */
    bounds: { min: Vector3; max: Vector3 } | null;
}

export interface ConvertResult {
    bytes: ArrayBuffer;
    outFormat: string;
}

export type Response =
    | { id: number; type: 'progress'; stage: string }
    | ({ id: number; type: 'result' } & (
          | ({ kind: 'init' } & InitResult)
          | ({ kind: 'render' } & RenderResult)
          | ({ kind: 'convert' } & ConvertResult)
          | { kind: 'closed' }
          | ({ kind: 'staged' } & StageResult)
          | ({ kind: 'summary' } & SummaryResult)
          | { kind: 'evicted' }
      ))
    | { id: number; type: 'error'; message: string };

/** Sensible starting parameters for each operation, used by the panel. */
export const OP_DEFAULTS: { [K in OpName]: Extract<OpSpec, { op: K }> } = {
    quality: { op: 'quality' },
    clean: {
        op: 'clean',
        weld: false,
        atol: 1e-8,
        removeOrphans: true,
        dropDegenerate: true,
        dropDuplicateCells: true,
    },
    smooth: {
        op: 'smooth',
        method: 'taubin',
        iterations: 20,
        // Negative means "this method's own default" on the C++ side, which is
        // 0.33 for Taubin and 0.5 for Laplacian — genuinely different values,
        // so the sentinel is better than picking one here.
        lambda: -1,
        mu: -0.34,
        fixBoundary: true,
    },
    refine: {
        op: 'refine',
        levels: 1,
        array: '',
        compare: '<',
        value: 0,
        closure: 'redgreen',
        recordLevels: false,
    },
    partition: { op: 'partition', nparts: 4, method: 'auto' },
    section: {
        op: 'section',
        point: [0, 0, 0],
        normal: [0, 0, 1],
        mode: 'all',
    },
    isosurface: { op: 'isosurface', array: '', isovalue: 0, component: -1 },
    gradient: {
        op: 'gradient',
        array: '',
        operator: 'gradient',
        method: 'green-gauss',
        location: 'point',
        output: '',
    },
    voxelize: { op: 'voxelize', resolution: [32, 32, 32], fill: 'surface' },
    computeSdf: {
        op: 'computeSdf',
        structure: 'voxel',
        resolution: [32, 32, 32],
        rootResolution: 8,
        maxDepth: 3,
    },
};

/** Human label for a pipeline chip. */
export function describeOp(spec: OpSpec): string {
    switch (spec.op) {
        case 'quality':
            return 'quality';
        case 'clean':
            return spec.weld ? `clean · weld ${spec.atol}` : 'clean';
        case 'smooth':
            return `smooth · ${spec.method} · ${spec.iterations}`;
        case 'refine':
            return spec.array
                ? `refine · ${spec.array} ${spec.compare} ${spec.value}`
                : `refine · ${spec.levels}×`;
        case 'partition':
            return `partition · ${spec.nparts} · ${spec.method}`;
        case 'section': {
            const axis = spec.normal.findIndex((v) => Math.abs(v) > 0.5);
            const name = ['X', 'Y', 'Z'][axis] ?? '?';
            const sign = (spec.normal[axis] ?? 1) < 0 ? '−' : '+';
            return `section · ${sign}${name}`;
        }
        case 'isosurface':
            return `isosurface · ${spec.array} = ${spec.isovalue}`;
        case 'gradient':
            return `${spec.operator} · ${spec.array}`;
        case 'voxelize':
            return `voxelize · ${spec.fill} ${spec.resolution[0]}³`;
        case 'computeSdf':
            return spec.structure === 'octree'
                ? `sdf · octree ${spec.rootResolution}³ ×${spec.maxDepth}`
                : `sdf · ${spec.resolution[0]}³`;
    }
}
