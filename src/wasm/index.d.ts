// Ambient declarations for @meshioplusplus/wasm's hand-written wrapper
// (src/index.mjs). Mirrors the JS-facing mesh object shape produced/consumed
// by bindings/wasm/js_bindings.cpp's meshToVal/valToMesh (see doc/wasm.md for
// the full format-support table and known v1 limitations).

/** A single homogeneous group of cells, all the same meshio++ cell type. */
export interface CellBlock {
  /** meshio++ cell type name, e.g. "triangle", "tetra10", "hexahedron". */
  type: string;
  /** Flat, row-major connectivity: length === numCells * nodesPerCell. */
  data: Int32Array;
  nodesPerCell: number;
}

/**
 * A mesh as exchanged with the WASM boundary: every array is copied (there
 * is no zero-copy path across the JS/WASM memory boundary, unlike the
 * Python bindings' numpy views) and cell connectivity is always Int32Array
 * (down-cast from the C++ core's Int64, which is safe for any mesh size a
 * browser can reasonably hold). Ragged cell blocks (polygon/polyhedron with
 * varying node counts) are not representable in this shape and are rejected
 * by both readMesh (throws) and writeMesh (cannot be constructed).
 */
export interface Mesh {
  /** Flat, row-major point coordinates: length === numPoints * dim. */
  points: Float64Array;
  /** 2 or 3. */
  dim: number;
  cells: CellBlock[];
  /** name -> flat, row-major per-point data. */
  point_data?: Record<string, Float64Array>;
  /** name -> one flat array per cell block, same order as `cells`. */
  cell_data?: Record<string, Float64Array[]>;
  /** name -> scalar/small metadata arrays (e.g. material ids). */
  field_data?: Record<string, Float64Array>;
}

export interface ConvertOptions {
  /** Explicit input format key, or omit to infer from inPath's extension. */
  inFormat?: string;
  /** Explicit output format key, or omit to infer from outPath's extension. */
  outFormat?: string;
}

/** One cell block's shape, as reported by {@link MeshioPlusPlusModule.readMetadata}. */
export interface MeshMetadataCellBlock {
  /** meshio++ cell type name, e.g. `"triangle"`, `"tetra10"`. */
  type: string;
  numCells: number;
  /** 0 for a ragged block, whose rows have no single node count. */
  nodesPerCell: number;
  ragged: boolean;
}

/**
 * A file's shape without its heavy arrays -- the result of `readMetadata`.
 *
 * `bboxMin`/`bboxMax` are **omitted** rather than null when no bounding box was
 * computed, so "not computed" cannot be misread as a box at the origin.
 */
export interface MeshMetadata {
  numPoints: number;
  pointDim: number;
  /** Total across every block. */
  numCells: number;
  cellBlocks: MeshMetadataCellBlock[];
  pointDataNames: string[];
  cellDataNames: string[];
  fieldDataNames: string[];
  /** The format that was actually used, whether given or inferred/sniffed. */
  format: string;
  /**
   * True when the format has no header-only path and the file had to be read
   * whole. The summary is still correct, just not cheap.
   */
  fellBackToFullRead: boolean;
  bboxMin?: number[];
  bboxMax?: number[];
}

/**
 * One operation in a {@link MeshioPlusPlusModule.convertSurfaceOps} pipeline.
 *
 * Parameters are optional and fall back to the same defaults the Python API
 * uses; only `op` is required.
 */
export type OpSpec =
  | { op: 'quality' }
  | {
      op: 'clean';
      weld?: boolean;
      atol?: number;
      removeOrphans?: boolean;
      dropDegenerate?: boolean;
      dropDuplicateCells?: boolean;
    }
  | {
      op: 'smooth';
      method?: SmoothMethod;
      iterations?: number;
      /** Negative means "this method's own default" (0.5 Laplacian, 0.33 Taubin). */
      lambda?: number;
      mu?: number;
      fixBoundary?: boolean;
    }
  | { op: 'refine'; levels?: number }
  | {
      /** Attaches the assignment as `partition:part` cell data. */
      op: 'partition';
      nparts?: number;
      method?: PartitionMethod;
    }
  | {
      /**
       * The planar cross-section through the mesh (slice), in **world**
       * coordinates: the volume is replaced by the surface where the plane
       * intersects it, one dimension lower — a genuine, flat, correctly-coloured
       * section. `mode` is accepted for backward compatibility but ignored (a
       * cross-section has no "keep side").
       */
      op: 'section';
      point: number[];
      normal: number[];
      mode?: CropMode;
    };

/** Per-operation counters and caveats from a pipeline run. */
export interface OpReport {
  steps: ({ op: OpSpec['op'] } & Record<string, number | string>)[];
  warnings: string[];
}

/** One data array's location: `point_data`, `cell_data`, or `field_data`. */
export type DataLocation = 'point' | 'cell' | 'field';

/** `data_condition`'s transform. See doc/data_condition.md. */
export type ConditionMode = 'clamp' | 'normalize' | 'standardize';

/** `data_condition`'s scope: independent components, or by row magnitude. */
export type ConditionScope = 'component' | 'magnitude';

/** Weighting for `dataCellToPoint`. See doc/data_average.md. */
export type CellPointWeight = 'uniform' | 'measure';

/** What reaches the output for non-finite (NaN/inf) values. They are always
 *  excluded from reductions regardless of this setting. */
export type NanPolicy = 'ignore' | 'replace' | 'fail';

/** Cell-keeping rule for the crop operations: every node inside, or any. */
export type CropMode = 'all' | 'any';

/** Partitioning criterion for `split`. */
export type SplitBy = 'type' | 'component' | 'region' | 'tag';

/** Node-renumbering method for `reorder`. */
export type ReorderMethod = 'rcm' | 'morton' | 'hilbert';

/** How `merge` reconciles data arrays that are not present in every input. */
export type MergeDataPolicy = 'intersection' | 'fill';

/** Element-representation conversion performed by `convertCells`. */
export type ConvertCellsMode = 'linearize' | 'simplexify' | 'elevate';

/** Smoothing operator applied by `smooth`. `'taubin'` is shrink-free. */
export type SmoothMethod = 'laplacian' | 'taubin';

/** How `interpolate` draws a target sample's value from the source.
 *  `'barycentric'` simplexifies the source first (simplex-linear on quad/hex
 *  sources; triangles evaluated in the xy-plane) and is exact on a linear
 *  field; `'nearest'` copies the nearest source point's value bit-for-bit. */
export type InterpolateMethod = 'nearest' | 'barycentric';

/** What `interpolate` does when a transferred name already exists on the
 *  target: throw, replace, or write to `name + '_interp'`. */
export type InterpolateOnConflict = 'error' | 'overwrite' | 'suffix';

/** Partitioning backend: SFC is always available; KaHIP is never compiled
 *  into the WASM build, so `'kahip'` always throws and `'auto'` = `'sfc'`. */
export type PartitionMethod = 'sfc' | 'kahip' | 'auto';

/** KaHIP preconfiguration (ignored by the SFC method). */
export type PartitionMode = 'fast' | 'eco' | 'strong';

/** Read-only per-array summary returned by `dataInfo`. See doc/data_info.md. */
export interface DataArrayInfo {
  location: 'point_data' | 'cell_data' | 'field_data';
  name: string;
  /** numpy-style dtype string, e.g. "f8", "i4". */
  dtype: string;
  shape: number[];
  /** cell_data: number of cell blocks; 1 otherwise. */
  numBlocks: number;
  numEntries: number;
  numComponents: number;
  numValues: number;
  /** Over finite values only; NaN when there are none. */
  min: number;
  max: number;
  mean: number;
  minPerComponent: number[];
  maxPerComponent: number[];
  meanPerComponent: number[];
  numNan: number;
  numInf: number;
  numFinite: number;
  /** cell_data whose blocks disagree in component count. */
  inconsistentBlocks: boolean;
}

/**
 * The instantiated module returned by `loadMeshioPlusPlus()`. `FS` is
 * Emscripten's virtual filesystem (MEMFS by default) -- write the bytes of a
 * mesh file there before calling `readMesh`, and read them back out after
 * `writeMesh`/`convert`. See https://emscripten.org/docs/api_reference/Filesystem-API.html
 */
export interface MeshioPlusPlusModule {
  FS: {
    writeFile(path: string, data: string | ArrayBufferView, opts?: object): void;
    readFile(path: string, opts?: { encoding?: 'binary' | 'utf8' }): Uint8Array | string;
    unlink(path: string): void;
    mkdir(path: string): void;
    [key: string]: unknown;
  };

  /**
   * Read a mesh file from the virtual filesystem.
   * @param path virtual FS path.
   * @param format explicit format key (see doc/wasm.md's table), or omit to
   *   infer from `path`'s extension. `.msh` defaults to gmsh, `.inp` to
   *   abaqus -- pass `format` explicitly to select ansys/freefem/ansysinp.
   * @throws {Error} on an unknown/unsupported format or a malformed file.
   */
  readMesh(path: string, format?: string): Mesh;

  /**
   * Read only part of a file: geometry alone, or a named subset of data arrays.
   *
   * `arrays` omitted/null reads every array; `arrays: []` reads none. Readers
   * with a native selective path (vtu/vtp/xdmf/gmsh) skip the unwanted arrays
   * outright; the rest are read whole and filtered, so the result is identical
   * either way and only the cost differs.
   */
  readMeshSelective(
    path: string,
    options?: { format?: string; pointsOnly?: boolean; arrays?: string[] | null }
  ): Mesh;

  /**
   * Summarize a mesh file without loading its heavy arrays.
   *
   * `bboxMin`/`bboxMax` are present only when a bounding box was computed --
   * omitted rather than null, so "not computed" cannot read as a box at the
   * origin. `fellBackToFullRead` is true when the format has no header-only
   * path: the summary is still correct, just not cheap.
   */
  readMetadata(path: string, format?: string): MeshMetadata;

  /** Whether `format` has a native selective-read path. */
  readerSupportsOptions(format: string): boolean;

  /**
   * Write a mesh to the virtual filesystem.
   * @throws {Error} on an unknown/write-unsupported format or malformed input
   *   (e.g. a points/connectivity array length not divisible by its
   *   declared dim/nodesPerCell).
   */
  writeMesh(path: string, mesh: Mesh, format?: string): void;

  /**
   * Read `inPath` and write it to `outPath` directly (no intermediate JS
   * mesh object). Mirrors the CLI's `convert` subcommand.
   */
  convert(inPath: string, outPath: string, options?: ConvertOptions): void;

  /**
   * Like {@link convert}, but writes a *renderable surface*: a mesh with
   * skinnable 3D cells becomes its boundary, anything else passes through, and
   * the result is linearized (a triangle renderer has no concept of a mid-side
   * node, so `triangle6` connectivity drawn verbatim is visible garbage).
   *
   * Prefer this over `readMesh` -> `extractSkin` -> `writeMesh` for anything
   * headed to a renderer. It never materializes a JS {@link Mesh}, so
   * multi-component (vector/tensor) arrays survive -- the flat JS
   * representation cannot carry them.
   */
  convertSurface(inPath: string, outPath: string, options?: ConvertOptions): void;

  /**
   * Like {@link convertSurface}, but applies a pipeline of mesh operations
   * first — all inside C++.
   *
   * Chaining the individual operation bindings (`smooth`, `clean`, …) would
   * route the mesh through the flat JS {@link Mesh} on every step and so
   * destroy every multi-component array; nothing crosses the boundary here.
   *
   * An **empty** pipeline is exactly {@link convertSurface}. That is
   * deliberate: one call serves both the plain and the post-operation display,
   * so they cannot drift, and undo becomes a replay of a shortened pipeline
   * rather than a set of inverse operations.
   *
   * Operations run on the **full-dimensional** mesh, before the boundary is
   * extracted — smoothing a solid is not the same as smoothing its skin.
   *
   * @throws {Error} on an unknown operation name or an unreadable file.
   */
  convertSurfaceOps(
    inPath: string,
    outPath: string,
    ops?: OpSpec[],
    options?: ConvertOptions & {
      /**
       * Keep the `surface:parent_cell` provenance array in the output. A
       * picker needs it; a colour-by menu must filter it out.
       */
      keepProvenance?: boolean;
    }
  ): OpReport;

  /** The shared cell-type -> node-count table (e.g. `{triangle: 3, tetra: 4, ...}`). */
  numNodesPerCell(): Record<string, number>;

  /** The shared cell-type -> topological-dimension table (0-3). */
  topologicalDimension(): Record<string, number>;

  /** The mesh backend this build was compiled with ("meshio"/"native"/"kratos"). */
  meshBackend(): string;

  /**
   * The format names this build can actually read and write, both sorted.
   *
   * Prefer this over a hardcoded table when building a file-picker filter or a
   * "convert to" menu: the wasm build carries no HDF5/netCDF-backed formats,
   * and the two lists genuinely differ -- `openfoam` is read-only, `svg` and
   * `tikz` are write-only.
   */
  availableFormats(): { readers: string[]; writers: string[] };

  /**
   * Mesh operations: computations ON a mesh (not file formats). The index maps
   * the C++ core returns for `cropBbox`/`cropPlane`/`split`/`convertCells` are
   * not carried across the JS boundary -- use the `recordIds`/`recordParentIds`
   * flags to get the same provenance as data arrays instead.
   */

  /** Extract the boundary of a mesh (faces of a volume, edges of a surface). */
  extractSurface(mesh: Mesh, recordParentIds?: boolean): Mesh;

  /** Extract the skin of a volume mesh (the volume-only special case). */
  extractSkin(mesh: Mesh, linearize?: boolean): Mesh;

  /** Return a clone with the quality metrics attached as `quality:*` cell_data. */
  attachQuality(mesh: Mesh): Mesh;

  /** Guess a file's format from its leading bytes, or `""` if ambiguous. */
  sniffFormat(path: string): string;

  /** Renumber nodes/cells to reduce bandwidth or improve locality. */
  reorder(mesh: Mesh, method?: ReorderMethod): { mesh: Mesh; nodePermutation: Int32Array; cellPermutations: Int32Array[] };

  /** Connectivity bandwidth: max over cells of (max node index - min node index). */
  computeBandwidth(mesh: Mesh): number;

  /** Structured comparison of two meshes. */
  diff(a: Mesh, b: Mesh, atol?: number, rtol?: number, unordered?: boolean): object;

  /** Whether two meshes are equal within tolerance. */
  meshesEqual(a: Mesh, b: Mesh, atol?: number, rtol?: number, unordered?: boolean): boolean;

  /** Combine several meshes into one, optionally welding coincident points. */
  merge(
    meshes: Mesh[],
    weld?: boolean,
    atol?: number,
    sourceTag?: boolean,
    dataPolicy?: MergeDataPolicy,
    dropDuplicateCells?: boolean,
  ): Mesh;

  /** Apply a row-major 4x4 affine transform to the point coordinates. */
  transform(mesh: Mesh, matrix: number[], rotateVectorData?: boolean): Mesh;

  /** Weld / prune / de-duplicate in one pass. */
  clean(
    mesh: Mesh,
    weld?: boolean,
    atol?: number,
    removeOrphans?: boolean,
    dropDegenerate?: boolean,
    dropDuplicateCells?: boolean,
  ): {
    mesh: Mesh;
    pointsWelded: number;
    pointsRemovedOrphan: number;
    cellsDroppedDegenerate: number;
    cellsDroppedDuplicate: number;
  };

  /**
   * Relax point coordinates toward their edge-neighbour centroids, leaving
   * connectivity and every data array untouched. `"taubin"` (the default)
   * alternates a `+lambda` and a `-mu` pass per iteration and is shrink-free;
   * `"laplacian"` is stronger per pass but contracts the mesh. A **negative**
   * `lambda` means "this method's own default" (0.5 Laplacian, 0.33 Taubin).
   * Boundary and feature nodes are pinned by default, and `guardInversion`
   * rejects any move that would flip an incident cell.
   * @throws {Error} on an unknown `method`, a non-negative `lambda` outside
   *   `(0, 1)`, or a `"taubin"` `mu` that does not satisfy `mu < -lambda < 0`.
   */
  smooth(
    mesh: Mesh,
    method?: SmoothMethod,
    iterations?: number,
    lambda?: number,
    mu?: number,
    fixBoundary?: boolean,
    preserveFeatures?: boolean,
    featureAngle?: number,
    guardInversion?: boolean,
  ): { mesh: Mesh; numNodesMoved: number; maxDisplacement: number; numSkippedInversion: number };

  /**
   * Sample data arrays from `source` onto `target` (cross-mesh field
   * transfer). Returns a copy of the target — geometry, connectivity and its
   * own data preserved exactly — with the requested source arrays sampled onto
   * it: source point_data at the target's points, source cell_data at the
   * target's cell centroids (always by nearest source-cell centroid, whatever
   * the method). An empty `arrays` transfers every source point_data array;
   * cell_data transfers only when named. Under `'barycentric'` a target point
   * outside the source domain receives `defaultValue` unless `extrapolate`.
   * @throws {Error} on an unknown method/onConflict, an unknown array name, a
   *   name collision under `'error'`, or a barycentric source with no
   *   triangle/tetrahedron cells after simplexification.
   */
  interpolate(
    source: Mesh,
    target: Mesh,
    method?: InterpolateMethod,
    arrays?: string[],
    extrapolate?: boolean,
    defaultValue?: number,
    onConflict?: InterpolateOnConflict,
  ): Mesh;

  /** Subset a mesh to an axis-aligned bounding box. */
  cropBbox(mesh: Mesh, lo: number[], hi: number[], mode?: CropMode, recordIds?: boolean): Mesh;

  /** Subset a mesh to the half-space `(p - point) . normal >= 0`. */
  cropPlane(
    mesh: Mesh,
    point: number[],
    normal: number[],
    mode?: CropMode,
    recordIds?: boolean,
  ): Mesh;

  /**
   * Planar cross-section of a mesh (marching tetrahedra on a simplexified
   * input): a volume mesh yields a triangle/quad surface, a 2D surface mesh a
   * line mesh. Crossing points on shared edges are deduped so the section is
   * watertight; each section cell inherits its parent's cell_data.
   */
  slice(mesh: Mesh, origin: number[], normal: number[], recordParentIds?: boolean): Mesh;

  /** Partition a mesh into submeshes by type, connected component, or tag. */
  split(mesh: Mesh, by: SplitBy, tagName?: string): { key: string; mesh: Mesh }[];

  /**
   * Convert the element representation: drop higher-order nodes
   * (`"linearize"`), decompose into same-dimension simplices (`"simplexify"`),
   * or promote linear cells to serendipity quadratic (`"elevate"`).
   * @throws {Error} on a polyhedron block under `"simplexify"`, or a
   *   full-Lagrange target (quad9/hexahedron27) under `"elevate"`.
   */
  convertCells(mesh: Mesh, mode?: ConvertCellsMode, recordParentIds?: boolean): Mesh;

  /**
   * Uniformly refine a mesh, subdividing every cell into same-type children
   * (`line` → 2, `triangle` → 4, `quad` → 4, `tetra` → 8, `wedge` → 8,
   * `hexahedron` → 8). New nodes sit at edge / quad-face / body midpoints and
   * are shared between neighbouring cells, so the result has no hanging nodes.
   * `levels` applies the templates repeatedly; `0` returns an unchanged copy.
   * @throws {Error} on a higher-order cell, a `pyramid`, or a ragged
   *   polygon/polyhedron block — none has a same-type subdivision.
   */
  refine(mesh: Mesh, levels?: number, recordParentIds?: boolean): Mesh;

  /**
   * Decompose a mesh into exactly `nparts` balanced pieces for domain
   * decomposition (the count-driven complement to `split`). Pieces keep the
   * input's cell-block structure 1:1 (empty blocks included, unlike `split`),
   * so concatenating them reproduces the input. The index maps are not
   * carried across the JS boundary — use `recordIds` for the
   * `partition:original_*_id` arrays, or `partitionLabels` for the raw
   * assignment. `weightsKey` names a scalar `cell_data` array of per-cell
   * weights. `ghostLayers` is reserved and must be 0.
   * @throws {Error} on `method: 'kahip'` (KaHIP is never part of the WASM
   *   build; the message names `MESHIOPLUSPLUS_WITH_KAHIP`), `nparts < 1`,
   *   `ghostLayers != 0`, or a bad weights array.
   */
  partition(
    mesh: Mesh,
    nparts: number,
    method?: PartitionMethod,
    imbalance?: number,
    mode?: PartitionMode,
    seed?: number,
    recordIds?: boolean,
    ghostLayers?: number,
    weightsKey?: string,
  ): { partId: number; mesh: Mesh }[];

  /**
   * The per-cell part assignment only: one array per cell block
   * (block-aligned, like a `cell_data` entry), values in `[0, nparts)`.
   */
  partitionLabels(
    mesh: Mesh,
    nparts: number,
    method?: PartitionMethod,
    imbalance?: number,
    mode?: PartitionMode,
    seed?: number,
    weightsKey?: string,
  ): number[][];

  /** Read-only geometric statistics (bbox, areas, volumes, inverted count). */
  stats(mesh: Mesh): object;

  /**
   * Data operations: act on `point_data`/`cell_data`/`field_data` only -- the
   * geometry is never modified. See doc/data_operations.md. Note the mesh
   * object's flat, shapeless array shape means these only support scalar
   * (1-component) arrays; a vector/tensor field cannot be round-tripped
   * through the WASM boundary at all.
   */

  /** Drop the named arrays at `location`. Empty `names` is a no-op. */
  dataDrop(mesh: Mesh, location: DataLocation, names?: string[], ignoreMissing?: boolean): Mesh;

  /** Keep only the named arrays at `location`; other locations are untouched. */
  dataKeep(mesh: Mesh, location: DataLocation, names?: string[], ignoreMissing?: boolean): Mesh;

  /** Rename one array, preserving its values and dtype. */
  dataRename(mesh: Mesh, location: DataLocation, from: string, to: string): Mesh;

  /** Average point_data onto the cells (mean over each cell's own nodes). */
  dataPointToCell(mesh: Mesh, names?: string[], suffix?: string): Mesh;

  /** Average cell_data onto the points, optionally measure-weighted. */
  dataCellToPoint(mesh: Mesh, names?: string[], weight?: CellPointWeight, suffix?: string): Mesh;

  /**
   * Evaluate an elementwise expression (`+ - * /`, parens, numbers, array
   * names, `abs`/`sqrt`/`min`/`max`/`norm`) and store the result as a new
   * array at `location`.
   * @throws {Error} on any lexical, syntactic, or name-resolution error.
   */
  dataCalc(
    mesh: Mesh,
    expression: string,
    location: DataLocation,
    outputName: string,
    overwrite?: boolean,
  ): Mesh;

  /** Clamp / normalize / standardize the selected arrays' values. */
  dataCondition(
    mesh: Mesh,
    location: DataLocation,
    names?: string[],
    mode?: ConditionMode,
    lo?: number,
    hi?: number,
    scope?: ConditionScope,
    nanPolicy?: NanPolicy,
    nanReplacement?: number,
    suffix?: string,
  ): Mesh;

  /** Read-only per-array summary of every data array the mesh carries. */
  dataInfo(mesh: Mesh): DataArrayInfo[];
}

/**
 * Instantiate a fresh, independent meshio++ WASM module instance. Safe to
 * call more than once (e.g. one instance per Web Worker).
 * @param moduleOverrides forwarded as-is to the Emscripten module factory
 *   (e.g. `{ locateFile }` to relocate the `.wasm` binary for a bundler/CDN).
 */
export function loadMeshioPlusPlus(moduleOverrides?: object): Promise<MeshioPlusPlusModule>;

export default loadMeshioPlusPlus;
