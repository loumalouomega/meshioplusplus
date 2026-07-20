// Ambient declarations for @meshioplusplus/wasm's hand-written wrapper
// (src/index.mjs). Mirrors the JS-facing mesh object shape produced/consumed
// by bindings_js/js_bindings.cpp's meshToVal/valToMesh (see doc/wasm.md for
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

  /** The shared cell-type -> node-count table (e.g. `{triangle: 3, tetra: 4, ...}`). */
  numNodesPerCell(): Record<string, number>;

  /** The shared cell-type -> topological-dimension table (0-3). */
  topologicalDimension(): Record<string, number>;

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
