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
}

/**
 * Instantiate a fresh, independent meshio++ WASM module instance. Safe to
 * call more than once (e.g. one instance per Web Worker).
 * @param moduleOverrides forwarded as-is to the Emscripten module factory
 *   (e.g. `{ locateFile }` to relocate the `.wasm` binary for a bundler/CDN).
 */
export function loadMeshioPlusPlus(moduleOverrides?: object): Promise<MeshioPlusPlusModule>;

export default loadMeshioPlusPlus;
