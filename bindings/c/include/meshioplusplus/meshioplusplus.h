/*  ██████   ██████ ██████████  █████████  █████   █████ █████    ███████
 * ░░██████ ██████ ░░███░░░░░█ ███░░░░░███░░███   ░░███ ░░███   ███░░░░░███      ███         ███
 *  ░███░█████░███  ░███  █ ░ ░███    ░░░  ░███    ░███  ░███  ███     ░░███    ░███        ░███
 *  ░███░░███ ░███  ░██████   ░░█████████  ░███████████  ░███ ░███      ░███ ███████████ ███████████
 *  ░███ ░░░  ░███  ░███░░█    ░░░░░░░░███ ░███░░░░░███  ░███ ░███      ░███░░░░░███░░░ ░░░░░███░░░
 *  ░███      ░███  ░███ ░   █ ███    ░███ ░███    ░███  ░███ ░░███     ███     ░███        ░███
 *  █████     █████ ██████████░░█████████  █████   █████ █████ ░░░███████░      ░░░         ░░░
 * ░░░░░     ░░░░░ ░░░░░░░░░░  ░░░░░░░░░  ░░░░░   ░░░░░ ░░░░░    ░░░░░░░
 *
 *
 *  License:         MIT License
 *                   meshio++ default license: LICENSE
 *
 *  Main authors:    Vicente Mataix Ferrandiz
 */

/**
 * @file meshioplusplus.h
 * @brief The meshio++ C API: a stable, pure-C99 interface to the C++ mesh
 *        I/O core, shipped as the `libmeshioplusplus` shared library. This is
 *        the only installed header -- the C++ headers behind it make no ABI
 *        promise.
 *
 * Conventions (the whole contract in five rules):
 *  1. Every fallible function returns a `mio_status` (`MIO_OK == 0`) or, for
 *     pointer/count returns, `NULL`/`-1`; the failure message is retrievable
 *     via mio_last_error() (thread-local, valid until the next mio_* call on
 *     the same thread). No C++ exception ever crosses this ABI.
 *  2. Setters COPY caller memory; the caller's buffers can be freed as soon
 *     as the call returns.
 *  3. Getters are ZERO-COPY: returned data pointers alias mesh-owned memory
 *     and stay valid until the next mutating `mio_mesh_*` call on that mesh
 *     or mio_mesh_free() -- read-only accessors never invalidate them.
 *  4. Arrays are row-major (C order). Points are `(num_points, dim)`,
 *     connectivity `(num_cells, nodes_per_cell)` with 0-based node indices.
 *  5. String getters use the snprintf convention: they copy at most
 *     `buflen - 1` bytes plus a NUL into `buf` (when `buflen > 0`) and
 *     return the full length excluding the NUL, or -1 on error.
 *
 * File I/O (`mio_read`/`mio_write`/`mio_convert`) infers the format from the
 * path's extension when `format` is NULL or ""; ambiguous extensions default
 * to `.msh` -> gmsh and `.inp` -> abaqus (pass "ansys"/"freefem"/"ansysinp"
 * explicitly instead). Formats backed by optional dependencies (cgns, h5m,
 * hmf, med need HDF5; exodus needs netCDF) exist only in builds configured
 * with them -- probe with mio_format_readable()/mio_format_writable().
 */

#ifndef MESHIOPLUSPLUS_H
#define MESHIOPLUSPLUS_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#if defined(_WIN32) && defined(MIO_SHARED)
#ifdef MIO_BUILDING
#define MIO_API __declspec(dllexport)
#else
#define MIO_API __declspec(dllimport)
#endif
#elif defined(__GNUC__)
#define MIO_API __attribute__((visibility("default")))
#else
#define MIO_API
#endif

/** Opaque mesh handle. Create with mio_mesh_create()/mio_read(); destroy
 *  with mio_mesh_free(). Not thread-safe per handle (distinct handles may be
 *  used from distinct threads freely). */
typedef struct mio_mesh mio_mesh;

/** Opaque result of mio_reorder(): owns the renumbered mesh plus the applied
 *  node/cell permutations. Destroy with mio_reorder_result_free(). */
typedef struct mio_reorder_result mio_reorder_result;

/** Opaque result of mio_diff(): the structured mesh comparison report. Destroy
 *  with mio_diff_result_free(). */
typedef struct mio_diff_result mio_diff_result;

/** Opaque result of mio_split(): the ordered submesh pieces. Destroy with
 *  mio_split_result_free(). */
typedef struct mio_split_result mio_split_result;

/** Opaque result of mio_convert_cells(): the converted mesh plus the point/cell
 *  index maps. Destroy with mio_convert_cells_result_free(). */
typedef struct mio_convert_cells_result mio_convert_cells_result;

/** Opaque result of mio_refine(): the refined mesh plus the point/cell index
 *  maps. Destroy with mio_refine_result_free(). */
typedef struct mio_refine_result mio_refine_result;

/** Opaque result of mio_decimate(): the decimated mesh plus the point/cell
 *  index maps and the collapse summary. Destroy with
 *  mio_decimate_result_free(). */
typedef struct mio_decimate_result mio_decimate_result;

/** Opaque result of mio_partition(): the pieces plus their point/cell index
 *  maps. Destroy with mio_partition_result_free(). */
typedef struct mio_partition_result mio_partition_result;

typedef enum mio_status {
    MIO_OK = 0,           /**< success */
    MIO_ERR_READ = 1,     /**< file could not be parsed / read-side failure */
    MIO_ERR_WRITE = 2,    /**< mesh could not be serialized / write-side failure */
    MIO_ERR_INVALID_ARG = 3, /**< NULL handle/pointer, bad dtype, bad shape, ... */
    MIO_ERR_NOT_FOUND = 4,   /**< unknown format/data name/block index */
    MIO_ERR_UNSUPPORTED = 5, /**< valid but unsupported (e.g. ragged connectivity) */
    MIO_ERR_INTERNAL = 6     /**< unexpected failure; see mio_last_error() */
} mio_status;

/** Overall outcome of comparing two meshes (see mio_diff()), in increasing
 *  order of severity — mirrors meshioplusplus::DiffVerdict. */
typedef enum mio_diff_verdict {
    MIO_DIFF_IDENTICAL = 0,             /**< structurally equal, values bitwise-equal */
    MIO_DIFF_EQUAL_WITHIN_TOLERANCE = 1, /**< values drifted but within tolerance */
    MIO_DIFF_DIFFERENT = 2               /**< a structural or beyond-tolerance difference */
} mio_diff_verdict;

/** Element dtypes, matching numpy's fixed-width scalars. */
typedef enum mio_dtype {
    MIO_FLOAT32 = 0,
    MIO_FLOAT64 = 1,
    MIO_INT8 = 2,
    MIO_INT16 = 3,
    MIO_INT32 = 4,
    MIO_INT64 = 5,
    MIO_UINT8 = 6,
    MIO_UINT16 = 7,
    MIO_UINT32 = 8,
    MIO_UINT64 = 9
} mio_dtype;

/** Maximum rank of any data array crossing this API (out-param shape arrays
 *  must have at least this many elements). */
#define MIO_MAX_NDIM 8

/* X(EnumName) -- one entry per meshio++ cell type, in the exact order of the
 * C++ `meshioplusplus::CellType` enum (src/cpp/include/meshioplusplus/
 * cell_type.hpp). c_api.cpp static_asserts every entry and the terminal
 * MIO_CELL_Custom against the C++ enum, so any drift between the two lists
 * is a compile error, never a runtime mismatch. */
#define MIO_CELL_TYPES(X)                                                                        \
    X(Vertex) X(Line) X(Line3) X(Line4) X(Line5) X(Line6) X(Line7) X(Line8) X(Line9) X(Line10)  \
    X(Line11) X(Triangle) X(Triangle6) X(Triangle10) X(Triangle15) X(Triangle21) X(Triangle28)  \
    X(Triangle36) X(Triangle45) X(Triangle55) X(Triangle66) X(Quad) X(Quad8) X(Quad9) X(Quad16) \
    X(Quad25) X(Quad36) X(Quad49) X(Quad64) X(Quad81) X(Quad100) X(Quad121) X(Tetra) X(Tetra10) \
    X(Tetra20) X(Tetra35) X(Tetra56) X(Tetra84) X(Tetra120) X(Tetra165) X(Tetra220) X(Tetra286) \
    X(Hexahedron) X(Hexahedron20) X(Hexahedron24) X(Hexahedron27) X(Hexahedron64)               \
    X(Hexahedron125) X(Hexahedron216) X(Hexahedron343) X(Hexahedron512) X(Hexahedron729)        \
    X(Hexahedron1000) X(Hexahedron1331) X(Wedge) X(Wedge15) X(Wedge18) X(Wedge40) X(Wedge75)    \
    X(Wedge126) X(Wedge196) X(Wedge288) X(Wedge405) X(Wedge550) X(Pyramid) X(Pyramid13)         \
    X(Pyramid14) X(Polygon) X(Polyhedron) X(VtkLagrangeCurve) X(VtkLagrangeTriangle)            \
    X(VtkLagrangeQuadrilateral) X(VtkLagrangeTetrahedron) X(VtkLagrangeHexahedron)              \
    X(VtkLagrangeWedge) X(VtkLagrangePyramid)

/** Integer mirror of the meshio++ cell-type table. The string names (e.g.
 *  "tetra10") are the primary representation everywhere in this API; the
 *  enum exists for C-side switching and metadata lookup. `MIO_CELL_Custom`
 *  is the catch-all for names outside the table. */
typedef enum mio_cell_type {
#define MIO_CELL_TYPE_ENUM(Name) MIO_CELL_##Name,
    MIO_CELL_TYPES(MIO_CELL_TYPE_ENUM)
#undef MIO_CELL_TYPE_ENUM
        MIO_CELL_Custom
} mio_cell_type;

/**
 * What kind of entity a region groups. Mirrors meshioplusplus::RegionKind.
 * Declared here (rather than down with the rest of the Named regions section)
 * so mio_read_metadata's region accessors can use it too.
 */
typedef enum mio_region_kind {
    MIO_REGION_POINT = 0,
    MIO_REGION_CELL = 1,
    MIO_REGION_SIDE = 2
} mio_region_kind;

/** Fixed-size description of one region (see mio_regions_info). */
typedef struct mio_region_info {
    int32_t kind;         /**< a mio_region_kind */
    int32_t dim;          /**< topological dimension, or -1 if unspecified */
    int64_t tag;          /**< format-native integer id, or -1 if none */
    int64_t num_entries;  /**< number of grouped entities (rows) */
    int64_t stride;       /**< values per entry: 2 for SIDE, else 1 (0 where entries aren't
                            *   carried at all, as in mio_read_metadata_region_info) */
} mio_region_info;

/* ---------------------------------------------------------------------
 * Version / build introspection
 * --------------------------------------------------------------------- */

/*
 * The RELEASE version as compile-time macros, so a consumer can feature-detect
 * with the preprocessor:
 *
 *     #if MIO_VERSION_AT_LEAST(9, 5, 0)
 *         mio_refine_ex(mesh, &opts);      // new in 9.5.0
 *     #else
 *         mio_refine(mesh, 1, 0);
 *     #endif
 *
 * These describe the header you COMPILED against. mio_version() below reports
 * the library you actually LINKED against, and with a shared library the two can
 * genuinely differ -- which is why both exist. (A third question, "are my
 * headers binary-compatible with that library", is answered by the ABI version
 * and its link-time sentinel; see doc/abi.md.)
 *
 * MIO_VERSION is major*10000 + minor*100 + patch, so it orders exactly like the
 * release and can be compared directly when the macro is not expressive enough.
 * These are static_asserted against the C++ MESHIOPLUSPLUS_VERSION_* macros in
 * c_api.cpp, and CMake hard-fails at configure time if either disagrees with
 * project(... VERSION ...), so the copies cannot drift.
 */
#define MIO_VERSION_MAJOR 9
#define MIO_VERSION_MINOR 11
#define MIO_VERSION_PATCH 0
#define MIO_VERSION (MIO_VERSION_MAJOR * 10000 + MIO_VERSION_MINOR * 100 + MIO_VERSION_PATCH)

/** Whether the header being compiled against is at least major.minor.patch. */
#define MIO_VERSION_AT_LEAST(major, minor, patch) \
    (MIO_VERSION >= ((major) * 10000 + (minor) * 100 + (patch)))

/** Whether the header being compiled against is older than major.minor.patch. */
#define MIO_VERSION_BEFORE(major, minor, patch) (!MIO_VERSION_AT_LEAST(major, minor, patch))

/**
 * @return the meshio++ version string of the LINKED library, e.g. "9.6.0"
 *         (static storage). Compare with the MIO_VERSION_* macros above, which
 *         describe the header this call was compiled against.
 */
MIO_API const char* mio_version(void);

/** @return the compile-time mesh backend: "meshio", "native", or "kratos"
 *          (static storage). */
MIO_API const char* mio_mesh_backend(void);

/** @return 1 if `format` (e.g. "gmsh", "vtu", "med") is readable in this
 *          build, 0 otherwise. */
MIO_API int mio_format_readable(const char* format);

/** @return 1 if `format` is writable in this build, 0 otherwise (read-only
 *          formats like "openfoam" report 0). */
MIO_API int mio_format_writable(const char* format);

/** @return the failure message of the most recent failed mio_* call on this
 *          thread ("" if none). Valid until the next mio_* call on the same
 *          thread; never NULL. */
MIO_API const char* mio_last_error(void);

/* ---------------------------------------------------------------------
 * Cell-type metadata
 * --------------------------------------------------------------------- */

/** @return the meshio name of `t` (e.g. "tetra10"; "" for MIO_CELL_Custom or
 *          out-of-range values). Static storage. */
MIO_API const char* mio_cell_type_name(mio_cell_type t);

/** @return the enum value for a meshio cell-type name, or MIO_CELL_Custom if
 *          `name` is NULL or not in the table. */
MIO_API mio_cell_type mio_cell_type_from_name(const char* name);

/** @return the fixed nodes-per-cell of `t`, or -1 for variable-size types
 *          (polygon, polyhedron, VTK Lagrange) and MIO_CELL_Custom. */
MIO_API int mio_cell_type_num_nodes(mio_cell_type t);

/** @return the topological dimension (0-3) of `t`, or -1 for MIO_CELL_Custom. */
MIO_API int mio_cell_type_dimension(mio_cell_type t);

/* ---------------------------------------------------------------------
 * Lifecycle & file I/O
 * --------------------------------------------------------------------- */

/** @return a new empty mesh, or NULL on allocation failure. */
MIO_API mio_mesh* mio_mesh_create(void);

/** Destroy a mesh and every pointer previously handed out from it. NULL-safe. */
MIO_API void mio_mesh_free(mio_mesh* mesh);

/**
 * Read a mesh file.
 * @param path   filesystem path.
 * @param format explicit format name, or NULL/"" to infer from the extension.
 * @return the mesh, or NULL on failure (see mio_last_error()).
 */
MIO_API mio_mesh* mio_read(const char* path, const char* format);

/**
 * Options narrowing what a read materializes (selective / partial reads).
 *
 * ABI NOTE: this struct is part of the installed library's permanent ABI.
 * New fields may only be appended, replacing `reserved` capacity; never
 * reorder, resize or repurpose an existing field. Always zero-initialize
 * through mio_read_opts_init() rather than by hand, so fields added later
 * default sensibly in code compiled against an older header.
 */
typedef struct mio_read_opts {
    int points_only;  /**< read geometry only, skipping every data array */
    int metadata_only; /**< read the header/summary only (see mio_read_metadata) */
    /** NULL = every array; otherwise `num_arrays` names to keep. An array
     *  count of 0 with a non-NULL pointer means "no arrays" -- the distinction
     *  from NULL is deliberate and load-bearing. Names absent from the file
     *  are ignored. The pointed-to strings are copied during the call. */
    const char* const* arrays;
    int64_t num_arrays;
    int mmap_mode;      /**< 0 = auto, 1 = on, 2 = off */
    /** Which step of a multi-step file to materialize: 0 (the default) is the
     *  first, preserving the historical behaviour; negative counts from the end.
     *  Out of range fails the call rather than clamping. Honoured by formats
     *  carrying a time series (currently exodus); ignored by the rest.
     *  Takes one of the former `reserved` slots, keeping the struct's size and
     *  every preceding field's offset unchanged. */
    int64_t time_step;
    /** Nonzero downgrades "this reader cannot represent construct X" errors to
     *  a warning plus a skip, where the reader can do that and still return a
     *  correct mesh (currently mdpa's Table/Geometries/Mesh/Constraints blocks).
     *  It is NOT "ignore all errors": a malformed file, a truncated block or a
     *  bad node reference still fails the call. 0 (the default) preserves the
     *  historical behaviour. Takes one of the former `reserved` slots, keeping
     *  the struct's size and every preceding field's offset unchanged. */
    int64_t lenient;
    int64_t reserved[4]; /**< must be zero; room for additive growth */
} mio_read_opts;

/** Initialize `opts` to the defaults (read everything). Always use this. */
MIO_API void mio_read_opts_init(mio_read_opts* opts);

/**
 * Read a mesh file honouring `opts`. `mio_read` is exactly
 * `mio_read_ex(path, format, <defaults>)` and is unchanged.
 *
 * Formats without a native selective path are read in full; the options are
 * still honoured, just without the saving.
 * @return the mesh, or NULL on failure (see mio_last_error()).
 */
MIO_API mio_mesh* mio_read_ex(const char* path, const char* format, const mio_read_opts* opts);

/** Opaque file summary. Destroy with mio_read_metadata_free(). */
typedef struct mio_read_metadata mio_read_metadata;

/**
 * Summarize a mesh file without loading its heavy arrays.
 * @return a result handle (free with mio_read_metadata_free), or NULL on failure.
 */
MIO_API mio_read_metadata* mio_read_metadata_create(const char* path, const char* format);

/** @return number of points, or -1 on error. */
MIO_API int64_t mio_read_metadata_num_points(const mio_read_metadata* meta);
/** @return coordinate dimension, or -1 on error. */
MIO_API int64_t mio_read_metadata_point_dim(const mio_read_metadata* meta);
/** @return total cells over all blocks, or -1 on error. */
MIO_API int64_t mio_read_metadata_num_cells(const mio_read_metadata* meta);
/** @return number of cell blocks, or -1 on error. */
MIO_API int64_t mio_read_metadata_num_cell_blocks(const mio_read_metadata* meta);

/**
 * The `index`-th cell block's shape. Any out pointer may be NULL.
 * `num_nodes_per_cell` is 0 for a ragged block (`is_ragged` nonzero).
 */
MIO_API mio_status mio_read_metadata_cell_block(const mio_read_metadata* meta, int64_t index,
                                                int64_t* num_cells, int64_t* num_nodes_per_cell,
                                                int* is_ragged);

/**
 * Copy the `index`-th cell block's meshio type name into `buf`.
 * @return the required length excluding the NUL, or -1 on error.
 */
MIO_API int64_t mio_read_metadata_cell_block_type(const mio_read_metadata* meta, int64_t index,
                                                  char* buf, int64_t buflen);

/**
 * @return number of time steps the file records, or -1 on error.
 *
 * 0 for a format with no time concept. This is the count `mio_read_opts.time_step`
 * may name, so it is what makes a step request checkable before issuing it.
 */
MIO_API int64_t mio_read_metadata_num_time_values(const mio_read_metadata* meta);

/**
 * Copy the recorded time values into `out` (at most `count` of them).
 * @return the number written, or -1 on error.
 */
MIO_API int64_t mio_read_metadata_time_values(const mio_read_metadata* meta, double* out,
                                              int64_t count);

/**
 * @return number of named regions the file carries, or -1 on error.
 *
 * Populated whenever the summary came from an already-read mesh (every
 * fallback path, and any format -- like exodus -- that always falls back); 0
 * on a native metadata path that declined a full read, since none of the
 * formats with such a path (VTU/VTP/XDMF/Gmsh 4.1) currently map regions at
 * all, so this is never a wrong answer.
 */
MIO_API int64_t mio_read_metadata_num_regions(const mio_read_metadata* meta);

/**
 * Copy the `index`-th region's name into `buf`.
 * @return the required length excluding the NUL, or -1 on error.
 */
MIO_API int64_t mio_read_metadata_region_name(const mio_read_metadata* meta, int64_t index,
                                              char* buf, int64_t buflen);

/**
 * The `index`-th region's kind/dim/tag/entry count. Reuses `mio_region_info`'s
 * shape (see the Named regions section below); `stride` is meaningless here
 * (no entries are carried) and is always 0.
 */
MIO_API mio_status mio_read_metadata_region_info(const mio_read_metadata* meta, int64_t index,
                                                 mio_region_info* out);

/** @return number of data-array names at `location` (a mio_data_location), or -1. */
MIO_API int64_t mio_read_metadata_num_names(const mio_read_metadata* meta, int location);

/**
 * Copy the `index`-th data-array name at `location` into `buf`.
 * @return the required length excluding the NUL, or -1 on error.
 */
MIO_API int64_t mio_read_metadata_name(const mio_read_metadata* meta, int location, int64_t index,
                                       char* buf, int64_t buflen);

/**
 * The bounding box, when one was computed.
 * @return MIO_OK and fills `min`/`max` (3 doubles each, either may be NULL), or
 *         MIO_ERR_NOT_FOUND when no box is available -- which is the normal
 *         case for a native summary, since computing it would mean decoding the
 *         point coordinates and defeating the purpose.
 */
MIO_API mio_status mio_read_metadata_bbox(const mio_read_metadata* meta, double* min, double* max);

/**
 * Whether the whole file had to be read to produce this summary.
 * @return 1 when the format has no header-only path (the summary is still
 *         correct, just not cheap), 0 when a native path was used, -1 on error.
 */
MIO_API int mio_read_metadata_fell_back(const mio_read_metadata* meta);

/** Destroy a summary handle. Safe to call with NULL. */
MIO_API void mio_read_metadata_free(mio_read_metadata* meta);

/** @return 1 if `format` has a native selective-read path, 0 if not, -1 on error. */
MIO_API int mio_reader_supports_options(const char* format);

/** Write a mesh. `format` as in mio_read(). */
MIO_API mio_status mio_write(const char* path, const mio_mesh* mesh, const char* format);

/** ASCII/binary selection for mio_write_opts.encoding. */
typedef enum mio_write_encoding {
    MIO_ENCODING_DEFAULT = 0, /**< the format's own default (unchanged behaviour) */
    MIO_ENCODING_ASCII = 1,
    MIO_ENCODING_BINARY = 2
} mio_write_encoding;

/** Block compression codec for mio_write_opts.codec (vtu/vtp only). */
typedef enum mio_write_codec {
    MIO_CODEC_DEFAULT = 0, /**< leave the format's default in place */
    MIO_CODEC_NONE = 1,
    MIO_CODEC_ZLIB = 2,
    MIO_CODEC_LZ4 = 3,  /**< vtkLZ4DataCompressor; ParaView-readable */
    MIO_CODEC_ZSTD = 4  /**< a meshio++ extension; ParaView cannot read it */
} mio_write_codec;

/**
 * Options for how a mesh is written -- the symmetric counterpart of
 * mio_read_opts.
 *
 * ABI NOTE: identical discipline to mio_read_opts. New fields may only be
 * appended, replacing `reserved` capacity; never reorder, resize or repurpose
 * an existing field. Always zero-initialize through mio_write_opts_init().
 *
 * An option the target format cannot honour FAILS the call rather than being
 * ignored -- asking for zstd in a Gmsh file is a mistake worth reporting.
 */
typedef struct mio_write_opts {
    int encoding;      /**< a mio_write_encoding value */
    int codec;         /**< a mio_write_codec value; vtu/vtp only */
    /** printf-style float format for the ASCII writers that take one (e.g.
     *  ".16e"). NULL or empty keeps the writer's own default. Copied during
     *  the call. Currently honoured by flac3d. */
    const char* float_format;
    int64_t reserved[5]; /**< must be zero; room for additive growth */
} mio_write_opts;

/** Initialize `opts` to the defaults (exactly mio_write's behaviour). */
MIO_API void mio_write_opts_init(mio_write_opts* opts);

/**
 * Write a mesh honouring `opts`. `mio_write(path, mesh, format)` is exactly
 * `mio_write_ex(path, mesh, format, <defaults>)` and is unchanged.
 * @return MIO_OK, or an error code (see mio_last_error()).
 */
MIO_API mio_status mio_write_ex(const char* path, const mio_mesh* mesh, const char* format,
                                const mio_write_opts* opts);

/** Read `in_path` and immediately write it to `out_path` (the CLI's
 *  `convert`), without materializing a handle for the caller. */
MIO_API mio_status mio_convert(const char* in_path, const char* in_format, const char* out_path,
                               const char* out_format);

/* ---------------------------------------------------------------------
 * Mesh operations (computations on a mesh, not file formats)
 * --------------------------------------------------------------------- */

/**
 * Extract the boundary of a mesh's highest-dimension cells as a new mesh.
 * Volume cells -> boundary faces; a 2D surface mesh -> boundary edges.
 * @param mesh             input mesh.
 * @param record_parent_ids nonzero to attach an int64 cell_data
 *                          "surface:parent_cell" (owning input-cell index).
 * @return the boundary mesh, or NULL on failure (see mio_last_error()).
 */
MIO_API mio_mesh* mio_extract_surface(const mio_mesh* mesh, int record_parent_ids);

/**
 * Extract the boundary skin of a volume mesh as a new surface mesh.
 * @param mesh      input volume mesh.
 * @param linearize nonzero to emit only corner nodes (triangle/quad output).
 * @return the skin mesh, or NULL on failure (see mio_last_error()).
 */
MIO_API mio_mesh* mio_extract_skin(const mio_mesh* mesh, int linearize);

/**
 * Compute per-cell quality metrics and return a copy of the mesh with them
 * attached as cell_data (names "quality:<metric>"; read them back with the
 * cell-data accessors).
 * @param mesh input mesh.
 * @return the annotated mesh, or NULL on failure (see mio_last_error()).
 */
MIO_API mio_mesh* mio_attach_quality(const mio_mesh* mesh);

/**
 * Report aggregate quality counts. Any out-param may be NULL.
 * @param mesh            input mesh.
 * @param num_cells       total cells scored.
 * @param num_inverted    cells with negative signed volume/area.
 * @param num_degenerate  near-zero (degenerate) cells.
 */
MIO_API mio_status mio_quality_counts(const mio_mesh* mesh, int64_t* num_cells,
                                      int64_t* num_inverted, int64_t* num_degenerate);

/**
 * Guess a mesh file's format from its contents (magic-byte sniffing).
 * @param path   filesystem path to an existing, readable file.
 * @param buf    caller buffer for the format name (may be NULL to query length).
 * @param buflen size of `buf`.
 * @return the untruncated length of the format name (0 if undetermined), or -1
 *         on error; the name is written to `buf` NUL-terminated when it fits.
 */
MIO_API int64_t mio_sniff_format(const char* path, char* buf, int64_t buflen);

/* ---------------------------------------------------------------------
 * Mesh operations (computations on a mesh, not file formats)
 * --------------------------------------------------------------------- */

/**
 * Extract the boundary of a mesh's highest-dimension cells as a new mesh.
 * Volume cells -> boundary faces; a 2D surface mesh -> boundary edges.
 * @param mesh             input mesh.
 * @param record_parent_ids nonzero to attach an int64 cell_data
 *                          "surface:parent_cell" (owning input-cell index).
 * @return the boundary mesh, or NULL on failure (see mio_last_error()).
 */
MIO_API mio_mesh* mio_extract_surface(const mio_mesh* mesh, int record_parent_ids);

/**
 * Extract the boundary skin of a volume mesh as a new surface mesh.
 * @param mesh      input volume mesh.
 * @param linearize nonzero to emit only corner nodes (triangle/quad output).
 * @return the skin mesh, or NULL on failure (see mio_last_error()).
 */
MIO_API mio_mesh* mio_extract_skin(const mio_mesh* mesh, int linearize);

/**
 * Compute per-cell quality metrics and return a copy of the mesh with them
 * attached as cell_data (names "quality:<metric>"; read them back with the
 * cell-data accessors).
 * @param mesh input mesh.
 * @return the annotated mesh, or NULL on failure (see mio_last_error()).
 */
MIO_API mio_mesh* mio_attach_quality(const mio_mesh* mesh);

/**
 * Report aggregate quality counts. Any out-param may be NULL.
 * @param mesh            input mesh.
 * @param num_cells       total cells scored.
 * @param num_inverted    cells with negative signed volume/area.
 * @param num_degenerate  near-zero (degenerate) cells.
 */
MIO_API mio_status mio_quality_counts(const mio_mesh* mesh, int64_t* num_cells,
                                      int64_t* num_inverted, int64_t* num_degenerate);

/**
 * Guess a mesh file's format from its contents (magic-byte sniffing).
 * @param path   filesystem path to an existing, readable file.
 * @param buf    caller buffer for the format name (may be NULL to query length).
 * @param buflen size of `buf`.
 * @return the untruncated length of the format name (0 if undetermined), or -1
 *         on error; the name is written to `buf` NUL-terminated when it fits.
 */
MIO_API int64_t mio_sniff_format(const char* path, char* buf, int64_t buflen);

/**
 * Merge two or more meshes into one (concatenate, with optional welding of
 * coincident nodes). Cell blocks are merged by type; connectivity is offset so
 * indices stay valid; point_data/cell_data are combined per `data_policy`.
 * point_sets/cell_sets are not carried across the C ABI (as elsewhere).
 * @param meshes   array of `count` non-NULL input mesh handles.
 * @param count    number of input meshes (>= 1).
 * @param weld     nonzero to merge coincident nodes within `atol`.
 * @param atol     absolute coincidence tolerance for welding.
 * @param source_tag nonzero to add an int64 "source_mesh_id" cell_data.
 * @param data_policy 0 = intersection (drop partial keys), 1 = fill (NaN).
 * @param drop_duplicate_cells nonzero (with `weld`) to drop cells that become
 *                             identical after welding.
 * @return the merged mesh (free with mio_mesh_free), or NULL on failure.
 */
MIO_API mio_mesh* mio_merge(const mio_mesh* const* meshes, int64_t count, int weld, double atol,
                            int source_tag, int data_policy, int drop_duplicate_cells);

/**
 * Apply an affine transform to a mesh's point coordinates (translate / scale /
 * rotate / general matrix). Connectivity and data are carried through; only the
 * points change (plus optional vector/tensor point_data rotation).
 * @param mesh   input mesh.
 * @param matrix a row-major 4x4 affine matrix (16 doubles); point p maps to
 *               M * [p.x, p.y, p.z, 1].
 * @param rotate_vector_data nonzero to rotate vector (trailing dim 3) and tensor
 *               (trailing dim 9) point_data by the transform's 3x3 linear block.
 * @return the transformed mesh (free with mio_mesh_free), or NULL on failure.
 */
MIO_API mio_mesh* mio_transform(const mio_mesh* mesh, const double* matrix,
                                int rotate_vector_data);

/**
 * Clean a mesh in one pass: weld coincident points, drop degenerate/duplicate
 * cells, remove orphaned points. Each step is toggleable; the removal counts are
 * returned through the (optional) out-params. point_sets/cell_sets are not
 * carried across the C ABI.
 * @param mesh                    input mesh.
 * @param weld                    nonzero to fuse coincident points within `atol`.
 * @param atol                    absolute coincidence tolerance for welding.
 * @param remove_orphans          nonzero to drop points referenced by no cell.
 * @param drop_degenerate         nonzero to drop degenerate cells.
 * @param drop_duplicate_cells    nonzero to drop exact-duplicate cells.
 * @param points_welded           receives the number of points fused by welding.
 * @param points_removed_orphan   receives the number of orphaned points removed.
 * @param cells_dropped_degenerate receives the number of degenerate cells dropped.
 * @param cells_dropped_duplicate receives the number of duplicate cells dropped.
 * @return the cleaned mesh (free with mio_mesh_free), or NULL on failure.
 */
MIO_API mio_mesh* mio_clean(const mio_mesh* mesh, int weld, double atol, int remove_orphans,
                            int drop_degenerate, int drop_duplicate_cells, int64_t* points_welded,
                            int64_t* points_removed_orphan, int64_t* cells_dropped_degenerate,
                            int64_t* cells_dropped_duplicate);

/**
 * Smooth a mesh's point coordinates, leaving topology and data untouched. Only
 * the points move: connectivity, cell_data, field_data and point_data values are
 * carried through unchanged. The run summary is returned through the (optional)
 * out-params. The SmoothOptions `mFrozen` caller pin mask is not exposed across
 * the C ABI (documented flat-ABI gap).
 * @param mesh              input mesh.
 * @param method            "laplacian" or "taubin"; NULL means "taubin".
 * @param iterations        number of iterations (for Taubin one iteration is two
 *                          passes); zero or less returns an unchanged copy.
 * @param lambda            relaxation factor in (0, 1). **Negative means "this
 *                          method's own default"**: 0.5 for laplacian, 0.33 for
 *                          taubin.
 * @param mu                taubin only: the un-shrinking factor, which must
 *                          satisfy mu < -lambda < 0. Ignored by laplacian.
 * @param fix_boundary      nonzero to pin every node on a boundary facet.
 * @param preserve_features nonzero to additionally pin boundary nodes whose
 *                          incident boundary facets meet above feature_angle.
 * @param feature_angle     feature angle between boundary facet normals, in
 *                          degrees (0 = coplanar, 90 = the edge of a box).
 * @param guard_inversion   nonzero to reject any move that would flip the signed
 *                          measure of an incident cell, and count the rejection.
 * @param nodes_moved       receives the number of points whose net displacement
 *                          exceeded the relative move tolerance.
 * @param max_displacement  receives the largest net displacement, in mesh units.
 * @param skipped_inversion receives the number of (node, pass) events in which
 *                          the inversion guard rejected a move.
 * @return the smoothed mesh (free with mio_mesh_free), or NULL on failure.
 */
MIO_API mio_mesh* mio_smooth(const mio_mesh* mesh, const char* method, int iterations,
                             double lambda, double mu, int fix_boundary, int preserve_features,
                             double feature_angle, int guard_inversion, int64_t* nodes_moved,
                             double* max_displacement, int64_t* skipped_inversion);

/**
 * Sample data arrays from a source mesh onto a target mesh (cross-mesh field
 * transfer). Returns a copy of the target — geometry, connectivity and its own
 * data preserved exactly — with the requested source arrays sampled onto it:
 * source point_data at the target's points, source cell_data at the target's
 * cell centroids (always by nearest source-cell centroid, whatever the method).
 * Under "barycentric" the source is simplexified first, so on a quad/hex source
 * the result is the simplex-linear interpolant (triangles are evaluated in the
 * xy-plane); "nearest" copies the nearest source point's value bit-for-bit.
 * @param source        the mesh whose data is sampled.
 * @param target        the mesh receiving the samples.
 * @param method        "nearest" or "barycentric"; NULL means "nearest".
 * @param arrays        source array names to transfer, as an array of C strings
 *                      with an explicit count (the mio_data_* convention). NULL
 *                      or arrays_count <= 0 means every source point_data array;
 *                      cell_data transfers only when named explicitly.
 * @param arrays_count  number of entries in arrays (see above).
 * @param extrapolate   barycentric only: nonzero to give a target point outside
 *                      the source domain the nearest source point's value
 *                      instead of default_value.
 * @param default_value barycentric only: the fill value (every component) for a
 *                      target point outside the source domain when extrapolate
 *                      is zero.
 * @param on_conflict   what to do when a transferred name already exists on the
 *                      target: "error", "overwrite" or "suffix" (writes to
 *                      name + "_interp"); NULL means "error".
 * @return the target copy with the sampled arrays (free with mio_mesh_free), or
 *         NULL on failure.
 */
MIO_API mio_mesh* mio_interpolate(const mio_mesh* source, const mio_mesh* target,
                                  const char* method, const char* const* arrays,
                                  int64_t arrays_count, int extrapolate, double default_value,
                                  const char* on_conflict);

/**
 * Planar cross-section of a mesh: the intersection with a plane, one topological
 * dimension below the cut cells (a volume mesh yields a triangle/quad surface, a
 * 2D surface mesh yields a line mesh). The input is simplexified first
 * (marching tetrahedra); crossing points on shared edges are deduped so the
 * section is watertight, and each section cell inherits its parent's cell_data.
 * @param origin            a point on the cutting plane (3 doubles).
 * @param normal            the plane normal (3 doubles, non-zero; need not be
 *                          unit length).
 * @param record_parent_ids nonzero to attach an Int64 "slice:parent_cell"
 *                          cell_data array (per section cell, the global input
 *                          cell it was cut from).
 * @return the cross-section mesh (free with mio_mesh_free), empty when the plane
 *         misses the geometry, or NULL on failure.
 */
MIO_API mio_mesh* mio_slice(const mio_mesh* mesh, const double* origin, const double* normal,
                            int record_parent_ids);

/**
 * Isosurfaces / contours: the level set of a scalar point_data field, one
 * topological dimension below the cut cells (a volume mesh yields a
 * triangle/quad surface, a 2D surface mesh a line mesh). The data-driven
 * sibling of mio_slice, sharing its marching-tetrahedra cutter: crossing points
 * on shared edges are deduped so contours are watertight, faces are wound toward
 * increasing field, and a node exactly at the isovalue counts as being on the
 * positive side (so a plateau emits its boundary once, not twice).
 *
 * The field must be point_data: cell_data is piecewise constant and has no level
 * set, so naming one fails -- convert it with mio_data_cell_to_point first.
 * All contours land in the one returned mesh, tagged per cell with a Float64
 * "iso:value" and an Int64 "iso:index" (the ordinal in the ascending isovalue
 * list, which is the integer tag mio_split's "tag" criterion needs).
 * @param array_name        the point_data array to contour.
 * @param isovalues         the level values (n_isovalues doubles); sorted
 *                          ascending and de-duplicated before cutting.
 * @param n_isovalues       how many, at least one.
 * @param component         component of a multi-component array to contour;
 *                          negative for the row magnitude (in which case the
 *                          contoured value is approximate, not exact).
 * @param record_parent_ids nonzero to attach an Int64 "iso:parent_cell"
 *                          cell_data array (per contour cell, the global input
 *                          cell it was cut from).
 * @return the contour mesh (free with mio_mesh_free), empty when no isovalue
 *         crosses the field, or NULL on failure.
 */
MIO_API mio_mesh* mio_isosurface(const mio_mesh* mesh, const char* array_name,
                                 const double* isovalues, int n_isovalues, int component,
                                 int record_parent_ids);

/**
 * Differentiate a point_data field: its gradient, divergence or curl.
 *
 * The result is a copy of the input carrying one new array; geometry,
 * connectivity, regions and every existing array come through unchanged.
 *
 * Shapes: an nc-component input yields a gradient of 3*nc components, flat and
 * row-major as [component i][derivative j] at index i*3 + j, so a scalar gives
 * 3 and a 3-vector gives 9. Divergence gives 1 component and curl always 3;
 * both need an input of 2 or 3 components (a 2-component field reads as
 * (u, v, 0)). Output is always Float64.
 *
 * Cells that cannot be differentiated -- below the mesh's own topological
 * dimension, ragged, or a 3-D Lagrange type with no face table -- yield NaN and
 * are counted in num_skipped, never approximated. With "least-squares", a cell
 * whose neighbourhood is degenerate falls back to Green-Gauss and is counted in
 * num_fallback.
 *
 * @param mesh         input mesh.
 * @param array_name   name of the point_data array to differentiate. A
 *                     cell_data name is an error naming the fix, since a
 *                     piecewise-constant field has no derivative.
 * @param op           "gradient" (default), "divergence" or "curl"; NULL or ""
 *                     selects the default.
 * @param method       "green-gauss" (default) or "least-squares"; NULL or ""
 *                     selects the default.
 * @param location     "cell" (default) or "point"; NULL or "" selects the
 *                     default.
 * @param output_name  name for the produced array; NULL or "" uses
 *                     "<array_name>:<op>".
 * @param component    differentiate only this component of a multi-component
 *                     input (gradient only); negative means EVERY component.
 *                     Note this is the opposite of mio_isosurface's sentinel,
 *                     where negative means the row magnitude.
 * @param overwrite    nonzero to replace an existing array of the output name
 *                     instead of failing.
 * @param num_skipped  optional out: cells that produced a NaN row.
 * @param num_fallback optional out: least-squares cells that fell back.
 * @return the mesh carrying the produced array (free with mio_mesh_free), or
 *         NULL on failure.
 */
MIO_API mio_mesh* mio_gradient(const mio_mesh* mesh, const char* array_name, const char* op,
                               const char* method, const char* location, const char* output_name,
                               int component, int overwrite, int64_t* num_skipped,
                               int64_t* num_fallback);

/**
 * Crop a mesh to an axis-aligned bounding box (keep cells inside the box).
 * @param mesh       input mesh.
 * @param lo         box lower corner (3 doubles).
 * @param hi         box upper corner (3 doubles).
 * @param mode       0 = keep cells with ALL nodes inside, 1 = ANY node inside.
 * @param record_ids nonzero to attach crop:original_point_id / _cell_id data.
 * @return the cropped mesh (free with mio_mesh_free), or NULL on failure.
 */
MIO_API mio_mesh* mio_crop_bbox(const mio_mesh* mesh, const double* lo, const double* hi, int mode,
                                int record_ids);

/**
 * Crop a mesh to the half-space (p - point) . normal >= 0.
 * @param mesh       input mesh.
 * @param point      a point on the plane (3 doubles).
 * @param normal     the plane normal (3 doubles); the kept side is non-negative.
 * @param mode       0 = keep cells with ALL nodes inside, 1 = ANY node inside.
 * @param record_ids nonzero to attach crop:original_point_id / _cell_id data.
 * @return the cropped mesh (free with mio_mesh_free), or NULL on failure.
 */
MIO_API mio_mesh* mio_crop_plane(const mio_mesh* mesh, const double* point, const double* normal,
                                 int mode, int record_ids);

/**
 * Split a mesh into pieces by cell type, connected component, or integer
 * cell_data tag. point_sets/cell_sets are not carried across the C ABI.
 * @param mesh     input mesh.
 * @param by       "type", "component", or "region"/"tag".
 * @param tag_name for "region"/"tag", the integer cell_data name (NULL/"" =
 *                 auto-detect the first integer cell_data); ignored otherwise.
 * @return a result handle (free with mio_split_result_free), or NULL on failure.
 */
MIO_API mio_split_result* mio_split(const mio_mesh* mesh, const char* by, const char* tag_name);

/** Number of pieces in a split result, or -1 on error. */
MIO_API int64_t mio_split_result_count(const mio_split_result* result);

/**
 * Copy piece `index`'s key (type name / component index / tag value) into `buf`.
 * @return the untruncated key length, or -1 on error (string-getter protocol).
 */
MIO_API int64_t mio_split_result_key(const mio_split_result* result, int64_t index, char* buf,
                                     int64_t buflen);

/**
 * Borrow piece `index`'s mesh. Owned by the result: valid until
 * mio_split_result_free(result); do NOT pass it to mio_mesh_free().
 * @return the piece mesh, or NULL if out of range.
 */
MIO_API const mio_mesh* mio_split_result_mesh(const mio_split_result* result, int64_t index);

/**
 * Transfer ownership of piece `index`'s mesh out of the result.
 * @return a new owning mesh handle (free with mio_mesh_free), or NULL on error.
 */
MIO_API mio_mesh* mio_split_result_take_mesh(mio_split_result* result, int64_t index);

/** Free a split result (and any pieces it still owns). */
MIO_API void mio_split_result_free(mio_split_result* result);

/**
 * Convert the element representation of a mesh. point_sets/cell_sets are not
 * carried across the C ABI.
 * @param mode               "linearize" (higher-order cells -> their linear
 *                           base, pruning the orphaned nodes), "simplexify"
 *                           (decompose into same-dimension simplices), or
 *                           "elevate" (linear -> serendipity quadratic, adding
 *                           one node per unique edge).
 * @param record_parent_ids  nonzero to attach a convert:parent_cell cell_data
 *                           array of per-block source cell indices.
 * @return a result handle (free with mio_convert_cells_result_free), or NULL on
 *         failure — a polyhedron block under "simplexify" and the full-Lagrange
 *         targets (quad9/hexahedron27) under "elevate" both fail.
 */
MIO_API mio_convert_cells_result* mio_convert_cells(const mio_mesh* mesh, const char* mode,
                                                    int record_parent_ids);

/**
 * Borrow the converted mesh. Owned by the result: valid until
 * mio_convert_cells_result_free(result); do NOT pass it to mio_mesh_free().
 * @return the converted mesh, or NULL on error.
 */
MIO_API const mio_mesh* mio_convert_cells_result_mesh(const mio_convert_cells_result* result);

/**
 * Transfer ownership of the converted mesh out of the result.
 * @return a new owning mesh handle (free with mio_mesh_free), or NULL on error.
 */
MIO_API mio_mesh* mio_convert_cells_result_take_mesh(mio_convert_cells_result* result);

/**
 * Zero-copy borrow of the point map (int64, shape (num_points_in,), input point
 * index -> output index, -1 when the point was pruned). The pointer is valid
 * until the result is freed. Any out-param may be NULL.
 * @param result a convert-cells result.
 * @param data   receives a pointer into result-owned int64 memory.
 * @param dtype  receives MIO_INT64.
 * @param n      receives the input point count.
 */
MIO_API mio_status mio_convert_cells_result_point_map(const mio_convert_cells_result* result,
                                                      const void** data, mio_dtype* dtype,
                                                      int64_t* n);

/** Number of per-block cell maps in a convert-cells result, or -1 on error. */
MIO_API int64_t mio_convert_cells_result_num_cell_maps(const mio_convert_cells_result* result);

/**
 * Zero-copy borrow of cell block `block`'s cell map (int64, shape
 * (num_cells_in_block,), input cell -> the index of its FIRST child in the
 * corresponding output block). Under "linearize"/"elevate" there is one child
 * per parent, so the map is the identity; under "simplexify" a parent's
 * children occupy the contiguous range starting at this index. Any out-param
 * may be NULL.
 * @param result a convert-cells result.
 * @param block  cell-block index in [0, mio_convert_cells_result_num_cell_maps).
 * @param data   receives a pointer into result-owned int64 memory.
 * @param dtype  receives MIO_INT64.
 * @param n      receives the block's input cell count.
 */
MIO_API mio_status mio_convert_cells_result_cell_map(const mio_convert_cells_result* result,
                                                     int64_t block, const void** data,
                                                     mio_dtype* dtype, int64_t* n);

/** Free a convert-cells result (and the mesh it still owns). */
MIO_API void mio_convert_cells_result_free(mio_convert_cells_result* result);

/** How mio_refine_ex resolves the hanging nodes a partial refinement leaves. */
typedef enum mio_refine_closure {
    MIO_REFINE_CLOSURE_REDGREEN = 0,  /**< promote to the smallest admissible mask (local) */
    MIO_REFINE_CLOSURE_PROPAGATE = 1, /**< promote straight to a full split (spreads) */
    /**
     * Do not close at all: KEEP the hanging nodes and only enforce 2:1 balance.
     * The output is 1-irregular and NOT conforming; the constrained nodes are
     * reported in the refine:hanging point_data array.
     */
    MIO_REFINE_CLOSURE_BALANCED = 2
} mio_refine_closure;

/** The comparison a mio_refine_opts predicate selector applies. */
typedef enum mio_refine_compare {
    MIO_REFINE_LT = 0,
    MIO_REFINE_LE = 1,
    MIO_REFINE_GT = 2,
    MIO_REFINE_GE = 3,
    MIO_REFINE_EQ = 4,
    MIO_REFINE_NE = 5
} mio_refine_compare;

/**
 * Options for mio_refine_ex: which cells to refine and how to close up after.
 *
 * ABI NOTE: this struct is part of the installed library's permanent ABI. New
 * fields may only be appended, replacing `reserved` capacity; never reorder,
 * resize or repurpose an existing field. Always zero-initialize through
 * mio_refine_opts_init() rather than by hand, so fields added later default
 * sensibly in code compiled against an older header.
 *
 * All zero means "refine every cell once", i.e. exactly mio_refine(mesh, 1, 0).
 * At most ONE of `cells`, `region` and `predicate_array` may be given; two is an
 * error rather than a precedence rule.
 */
typedef struct mio_refine_opts {
    /** Global (block-major) indices of the cells to refine; NULL for none. */
    const int64_t* cells;
    /** Length of `cells`; ignored when `cells` is NULL. */
    int64_t num_cells;
    /**
     * Name of a region to refine, or NULL/"" for none. A cell region selects its
     * own cells; a point region selects every cell with ANY node in it; a side
     * region is an error, since a facet is not a cell.
     */
    const char* region;
    /**
     * Name of a scalar numeric cell_data array to threshold, or NULL/"" for
     * none. A non-finite cell value never matches.
     */
    const char* predicate_array;
    /** The predicate's right-hand side. */
    double predicate_value;
    /** How many times to apply the templates; 0 or less is an unchanged clone. */
    int32_t levels;
    /** Nonzero to attach the refine:parent_cell cell_data array. */
    int32_t record_parent_ids;
    /** Nonzero to attach the refine:level cell_data array. */
    int32_t record_levels;
    /** A mio_refine_closure. */
    int32_t closure;
    /** A mio_refine_compare. */
    int32_t predicate_op;
    int32_t reserved_pad; /**< must be zero; keeps the int64 tail aligned */
    int64_t reserved[6];  /**< must be zero; room for additive growth */
} mio_refine_opts;

/** Zero-initialize refine options (every field its default). */
MIO_API void mio_refine_opts_init(mio_refine_opts* opts);

/**
 * Refine a mesh, subdividing every cell into same-type children (line -> 2,
 * triangle -> 4, quad -> 4, tetra -> 8, wedge -> 8, hexahedron -> 8). New nodes
 * sit at edge / quad-face / body midpoints and are shared between neighbouring
 * cells, so the refined mesh has no hanging nodes.
 * point_sets/cell_sets are not carried across the C ABI.
 *
 * This is mio_refine_ex with every option at its default; use that to refine a
 * SUBSET of the cells instead.
 * @param levels             how many times to apply the templates; 0 or less
 *                           returns an unchanged clone.
 * @param record_parent_ids  nonzero to attach a refine:parent_cell cell_data
 *                           array naming each output cell's ORIGINAL ancestor
 *                           within its own block.
 * @return a result handle (free with mio_refine_result_free), or NULL on
 *         failure — higher-order cells, pyramids and ragged polygon/polyhedron
 *         blocks have no same-type subdivision and all fail.
 */
MIO_API mio_refine_result* mio_refine(const mio_mesh* mesh, int levels, int record_parent_ids);

/**
 * Refine a mesh with a cell selection and a conforming closure.
 *
 * With no selector set this is mio_refine. With one, only the selected cells get
 * the full template and the hanging nodes that leaves are resolved by the
 * closure, so the output is still a conforming mesh — see doc/refine.md for the
 * per-cell-type support matrix.
 * @param mesh  the mesh to refine (unchanged).
 * @param opts  the options, or NULL for every default.
 * @return a result handle (free with mio_refine_result_free), or NULL on
 *         failure — including more than one selector, an out-of-range cell
 *         index, an unknown region, a side region used as a selector, and a
 *         predicate array that is missing, does not cover every block or is not
 *         scalar.
 */
MIO_API mio_refine_result* mio_refine_ex(const mio_mesh* mesh, const mio_refine_opts* opts);

/**
 * Borrow the refined mesh. Owned by the result: valid until
 * mio_refine_result_free(result); do NOT pass it to mio_mesh_free().
 * @return the refined mesh, or NULL on error.
 */
MIO_API const mio_mesh* mio_refine_result_mesh(const mio_refine_result* result);

/**
 * Transfer ownership of the refined mesh out of the result.
 * @return a new owning mesh handle (free with mio_mesh_free), or NULL on error.
 */
MIO_API mio_mesh* mio_refine_result_take_mesh(mio_refine_result* result);

/**
 * Zero-copy borrow of the point map (int64, shape (num_points_in,), input point
 * index -> output index). Refinement never prunes, so this is the identity; it
 * is exposed in full so callers need not depend on that. The pointer is valid
 * until the result is freed. Any out-param may be NULL.
 * @param result a refine result.
 * @param data   receives a pointer into result-owned int64 memory.
 * @param dtype  receives MIO_INT64.
 * @param n      receives the input point count.
 */
MIO_API mio_status mio_refine_result_point_map(const mio_refine_result* result, const void** data,
                                               mio_dtype* dtype, int64_t* n);

/** Number of per-block cell maps in a refine result, or -1 on error. */
MIO_API int64_t mio_refine_result_num_cell_maps(const mio_refine_result* result);

/**
 * Zero-copy borrow of cell block `block`'s cell map (int64, shape
 * (num_cells_in_block,), input cell -> the index of its FIRST child in the
 * corresponding output block). A cell's children are contiguous, so parent c
 * owns [map[c], map[c+1]) (or the block's end for the last parent). Any
 * out-param may be NULL.
 * @param result a refine result.
 * @param block  cell-block index in [0, mio_refine_result_num_cell_maps).
 * @param data   receives a pointer into result-owned int64 memory.
 * @param dtype  receives MIO_INT64.
 * @param n      receives the block's input cell count.
 */
MIO_API mio_status mio_refine_result_cell_map(const mio_refine_result* result, int64_t block,
                                              const void** data, mio_dtype* dtype, int64_t* n);

/** Free a refine result (and the mesh it still owns). */
MIO_API void mio_refine_result_free(mio_refine_result* result);

/**
 * Decimate a SURFACE mesh by greedy quadric-error-metric (Garland-Heckbert)
 * edge collapse — the resolution-reducing inverse of mio_refine. The output is
 * all-triangle (quad/polygon blocks are triangulated first) with the block
 * structure kept 1:1. Boundary vertices (once-used-edge test) and feature
 * vertices (incident face normals differing by more than feature_angle
 * degrees) are pinned by default; the link condition and a normal-flip guard
 * reject any collapse that would change topology or fold the surface. The
 * caller frozen mask is not exposed across the C ABI (a documented flat-ABI
 * gap, like mio_smooth's).
 * @param target_ratio      fraction of the (triangulated) faces to KEEP, in
 *                          (0, 1]; negative = unset.
 * @param target_faces      absolute face count to stop at (the result lands
 *                          within one collapse of it); negative = unset.
 * @param max_error         collapse only while the cheapest candidate's
 *                          quadric error is at most this; negative = unset.
 *                          Exactly one of the three criteria must be set.
 * @param placement         "optimal" (quadric minimizer, midpoint when
 *                          ill-conditioned), "midpoint" or "endpoint"; NULL =
 *                          "optimal".
 * @param preserve_boundary nonzero to pin boundary vertices (default-on
 *                          behaviour: pass 1).
 * @param preserve_features nonzero to pin feature vertices.
 * @param feature_angle     the feature angle in degrees (30 is the
 *                          vtkFeatureEdges convention).
 * @return a result handle (free with mio_decimate_result_free), or NULL on
 *         failure — a 3D volume mesh (extract the surface first), higher-order
 *         or ragged blocks, line/vertex blocks, and a criterion count != 1 all
 *         fail by name.
 */
MIO_API mio_decimate_result* mio_decimate(const mio_mesh* mesh, double target_ratio,
                                          int64_t target_faces, double max_error,
                                          const char* placement, int preserve_boundary,
                                          int preserve_features, double feature_angle);

/**
 * Borrow the decimated mesh. Owned by the result: valid until
 * mio_decimate_result_free(result); do NOT pass it to mio_mesh_free().
 * @return the decimated mesh, or NULL on error.
 */
MIO_API const mio_mesh* mio_decimate_result_mesh(const mio_decimate_result* result);

/**
 * Transfer ownership of the decimated mesh out of the result.
 * @return a new owning mesh handle (free with mio_mesh_free), or NULL on error.
 */
MIO_API mio_mesh* mio_decimate_result_take_mesh(mio_decimate_result* result);

/**
 * Zero-copy borrow of the point map (int64, shape (num_points_in,), input
 * point index -> output index). A collapsed point maps to its SURVIVOR's
 * output index — which is what makes it usable for remapping external
 * per-point arrays — and is -1 only when the survivor itself ended up
 * unreferenced. The pointer is valid until the result is freed. Any out-param
 * may be NULL.
 * @param result a decimate result.
 * @param data   receives a pointer into result-owned int64 memory.
 * @param dtype  receives MIO_INT64.
 * @param n      receives the input point count.
 */
MIO_API mio_status mio_decimate_result_point_map(const mio_decimate_result* result,
                                                 const void** data, mio_dtype* dtype, int64_t* n);

/** Number of per-block cell maps in a decimate result, or -1 on error. */
MIO_API int64_t mio_decimate_result_num_cell_maps(const mio_decimate_result* result);

/**
 * Zero-copy borrow of INPUT cell block `block`'s cell map (int64, shape
 * (num_cells_in_block,), input cell -> the output index, within the
 * corresponding output block, of its first surviving triangle; -1 when none
 * survived). Any out-param may be NULL.
 * @param result a decimate result.
 * @param block  cell-block index in [0, mio_decimate_result_num_cell_maps).
 * @param data   receives a pointer into result-owned int64 memory.
 * @param dtype  receives MIO_INT64.
 * @param n      receives the block's input cell count.
 */
MIO_API mio_status mio_decimate_result_cell_map(const mio_decimate_result* result, int64_t block,
                                                const void** data, mio_dtype* dtype, int64_t* n);

/** Triangles removed (counted on the triangulated mesh), or -1 on error. */
MIO_API int64_t mio_decimate_result_faces_removed(const mio_decimate_result* result);

/** Points removed (collapsed plus pruned-unreferenced), or -1 on error. */
MIO_API int64_t mio_decimate_result_points_removed(const mio_decimate_result* result);

/** Guard-rejection events during the run, or -1 on error. */
MIO_API int64_t mio_decimate_result_collapses_rejected(const mio_decimate_result* result);

/** The largest committed collapse error (0.0 when nothing collapsed), or a
 *  negative value on error. */
MIO_API double mio_decimate_result_max_error_applied(const mio_decimate_result* result);

/** Free a decimate result (and the mesh it still owns). */
MIO_API void mio_decimate_result_free(mio_decimate_result* result);

/**
 * Decompose a mesh into `nparts` balanced pieces for domain decomposition (the
 * count-driven complement to mio_split). Methods: "sfc" (Hilbert space-filling
 * curve cut of the cell centroids — always available, deterministic), "kahip"
 * (KaHIP kaffpa() on the shared-face dual graph — only in builds with
 * -DMESHIOPLUSPLUS_WITH_KAHIP=ON; requesting it elsewhere fails with an error
 * naming that option), or "auto"/NULL/"" (kahip when built, else sfc). Every
 * piece keeps the input's cell-block structure 1:1 (empty blocks included,
 * unlike mio_split), so the cell maps index by input block and the pieces
 * recombine into the input. point_sets/cell_sets are not carried across the
 * C ABI.
 * @param nparts       number of pieces (>= 1); the result always has exactly
 *                     nparts pieces, possibly with empty meshes.
 * @param method       "sfc", "kahip", "auto" (NULL/"" = "auto").
 * @param imbalance    kahip only: allowed imbalance fraction in (0, 1)
 *                     (0.03 = 3%).
 * @param mode         kahip only: "fast", "eco", "strong" (NULL/"" = "eco").
 * @param seed         kahip only: random seed (deterministic per seed).
 * @param record_ids   nonzero to attach Int64 partition:original_point_id
 *                     point_data and partition:original_cell_id cell_data
 *                     (original input indices) to every piece.
 * @param ghost_layers grow each piece by this many shared-node BFS layers of
 *        other parts' cells (a halo), tagged with an Int64 `partition:ghost`
 *        cell_data (0 = owned, L = reached at layer L). 0 keeps the pieces
 *        disjoint; negative fails the call.
 * @param weights_key  name of a scalar numeric cell_data array of per-cell
 *                     weights, or NULL/"" for the unweighted rule.
 * @return a result handle (free with mio_partition_result_free), or NULL on
 *         failure (see mio_last_error()).
 */
MIO_API mio_partition_result* mio_partition(const mio_mesh* mesh, int nparts, const char* method,
                                            double imbalance, const char* mode, int seed,
                                            int record_ids, int ghost_layers,
                                            const char* weights_key);

/** Number of pieces in a partition result (== nparts), or -1 on error. */
MIO_API int64_t mio_partition_result_num_pieces(const mio_partition_result* result);

/** The part id of piece `index` (== index, in [0, nparts)), or -1 on error. */
MIO_API int mio_partition_result_part_id(const mio_partition_result* result, int64_t index);

/**
 * Borrow piece `index`'s mesh. Owned by the result: valid until
 * mio_partition_result_free(result); do NOT pass it to mio_mesh_free().
 * @return the piece mesh, or NULL on error.
 */
MIO_API const mio_mesh* mio_partition_result_mesh(const mio_partition_result* result,
                                                  int64_t index);

/**
 * Transfer ownership of piece `index`'s mesh out of the result.
 * @return a new owning mesh handle (free with mio_mesh_free), or NULL on error.
 */
MIO_API mio_mesh* mio_partition_result_take_mesh(mio_partition_result* result, int64_t index);

/**
 * Zero-copy borrow of piece `index`'s point map (int64, shape
 * (num_points_in,), input point index -> piece point index, -1 where the point
 * is not in the piece). The pointer is valid until the result is freed. Any
 * out-param may be NULL.
 * @param result a partition result.
 * @param index  piece index in [0, mio_partition_result_num_pieces).
 * @param data   receives a pointer into result-owned int64 memory.
 * @param dtype  receives MIO_INT64.
 * @param n      receives the input point count.
 */
MIO_API mio_status mio_partition_result_point_map(const mio_partition_result* result,
                                                  int64_t index, const void** data,
                                                  mio_dtype* dtype, int64_t* n);

/** Number of per-block cell maps in each piece (== the input's cell-block
 *  count), or -1 on error. */
MIO_API int64_t mio_partition_result_num_cell_maps(const mio_partition_result* result);

/**
 * Zero-copy borrow of piece `index`'s cell map for input block `block` (int64,
 * shape (num_cells_in_block,), input cell -> piece cell index within that
 * block, -1 where the cell is not in the piece). Any out-param may be NULL.
 * @param result a partition result.
 * @param index  piece index in [0, mio_partition_result_num_pieces).
 * @param block  cell-block index in [0, mio_partition_result_num_cell_maps).
 * @param data   receives a pointer into result-owned int64 memory.
 * @param dtype  receives MIO_INT64.
 * @param n      receives the block's input cell count.
 */
MIO_API mio_status mio_partition_result_cell_map(const mio_partition_result* result, int64_t index,
                                                 int64_t block, const void** data,
                                                 mio_dtype* dtype, int64_t* n);

/** Free a partition result (and every piece mesh it still owns). */
MIO_API void mio_partition_result_free(mio_partition_result* result);

/**
 * Compute the per-cell part assignment only, without building the pieces.
 * Fills a caller-allocated flat buffer in block-major global cell order (the
 * mesh's blocks concatenated in order — recover the block alignment from the
 * mesh's own per-block cell counts). Parameters as mio_partition (record_ids /
 * ghost_layers do not apply).
 * @param labels      caller buffer of labels_size int64 entries; receives
 *                    values in [0, nparts).
 * @param labels_size must equal the mesh's total cell count.
 * @return MIO_OK, or an error status (see mio_last_error()).
 */
MIO_API mio_status mio_partition_labels(const mio_mesh* mesh, int nparts, const char* method,
                                        double imbalance, const char* mode, int seed,
                                        const char* weights_key, int64_t* labels,
                                        int64_t labels_size);

/** Geometric statistics of a mesh (see mio_stats). Per-cell-type counts are not
 *  carried across the C ABI (use mio_num_cell_blocks / mio_cell_block_type). */
typedef struct mio_stats_report {
    int64_t num_points;     /**< number of points */
    int64_t num_cells;      /**< total number of cells */
    double bbox_min[3];     /**< bounding-box lower corner */
    double bbox_max[3];     /**< bounding-box upper corner */
    double extent[3];       /**< bounding-box extent (max - min) */
    double centroid[3];     /**< mean of the point coordinates */
    double total_area;      /**< area of 2D cells + boundary of 3D cells */
    double signed_volume;   /**< sum of signed volumes of 3D cells */
    double unsigned_volume; /**< sum of |volume| of 3D cells */
    int64_t num_inverted;   /**< 3D cells with negative signed volume */
} mio_stats_report;

/**
 * Compute geometric statistics of a mesh (read-only).
 * @param mesh input mesh.
 * @param out  receives the statistics (must be non-NULL).
 * @return MIO_OK, or an error status (see mio_last_error()).
 */
MIO_API mio_status mio_stats(const mio_mesh* mesh, mio_stats_report* out);

/* ------------------------------------------------------------------------- */
/* Data operations                                                           */
/*                                                                           */
/* These act on a mesh's point_data / cell_data / field_data arrays rather    */
/* than on its geometry, which they never modify. Each returns a NEW mesh     */
/* (free it with mio_mesh_free) or NULL on failure; see mio_last_error().     */
/*                                                                           */
/* Name lists cross as an array of C strings plus an explicit count (the same */
/* convention as mio_merge), not NULL-terminated -- easier for Fortran to     */
/* build. Passing count == 0 means "every array at that location" for the     */
/* averaging and conditioning entry points.                                   */
/*                                                                           */
/* Documented flat-ABI gap: the combined data_manage (keep + drop + rename in */
/* one call) is not exposed, since an array-of-triples-of-enum-and-two-string */
/* parameter is unwieldy and the three primitives below compose to the same   */
/* effect. Likewise point/cell SETS never reach the C++ core at all.           */
/* ------------------------------------------------------------------------- */

/** Which of a mesh's three data maps an array lives in. */
typedef enum mio_data_location {
    MIO_DATA_POINT = 0, /**< point_data -- one row per mesh point */
    MIO_DATA_CELL = 1,  /**< cell_data -- one array per cell block */
    MIO_DATA_FIELD = 2  /**< field_data -- mesh-global */
} mio_data_location;

/** Weighting for mio_data_cell_to_point. */
typedef enum mio_cell_point_weight {
    MIO_WEIGHT_UNIFORM = 0, /**< every incident cell contributes equally */
    MIO_WEIGHT_MEASURE = 1  /**< weight by |cell measure| (length/area/volume) */
} mio_cell_point_weight;

/** Conditioning transform for mio_data_condition. */
typedef enum mio_condition_mode {
    MIO_COND_CLAMP = 0,      /**< x -> min(max(x, lo), hi) */
    MIO_COND_NORMALIZE = 1,  /**< map the array's [min,max] onto [lo,hi] */
    MIO_COND_STANDARDIZE = 2 /**< (x - mean) / stddev */
} mio_condition_mode;

/** Whether conditioning treats components independently or by row magnitude. */
typedef enum mio_condition_scope {
    MIO_SCOPE_COMPONENT = 0, /**< each trailing component independently */
    MIO_SCOPE_MAGNITUDE = 1  /**< by row magnitude; direction preserved */
} mio_condition_scope;

/** What reaches the output for non-finite (NaN/inf) values. They are ALWAYS
 *  excluded from reductions regardless of this setting. */
typedef enum mio_nan_policy {
    MIO_NAN_IGNORE = 0,  /**< pass through unchanged */
    MIO_NAN_REPLACE = 1, /**< substitute the caller's replacement value */
    MIO_NAN_FAIL = 2     /**< fail with MIO_ERR_INVALID_ARG */
} mio_nan_policy;

/**
 * Drop the named data arrays at one location.
 * @param mesh           input mesh (unmodified).
 * @param location       which data map to drop from.
 * @param names          array of `count` names.
 * @param count          number of entries in `names`.
 * @param ignore_missing non-zero to skip names that do not exist.
 * @return a new mesh (free with mio_mesh_free), or NULL on failure.
 */
MIO_API mio_mesh* mio_data_drop(const mio_mesh* mesh, mio_data_location location,
                                const char* const* names, int64_t count, int ignore_missing);

/**
 * Keep only the named data arrays at one location, dropping the rest there.
 * The other two locations are left untouched; `count == 0` drops everything at
 * `location`.
 * @return a new mesh (free with mio_mesh_free), or NULL on failure.
 */
MIO_API mio_mesh* mio_data_keep(const mio_mesh* mesh, mio_data_location location,
                                const char* const* names, int64_t count, int ignore_missing);

/**
 * Rename one data array, preserving its values, dtype and shape.
 * @return a new mesh (free with mio_mesh_free), or NULL on failure.
 */
MIO_API mio_mesh* mio_data_rename(const mio_mesh* mesh, mio_data_location location,
                                  const char* from_name, const char* to_name);

/**
 * Average point_data onto the cells: each cell's value is the mean over its own
 * nodes. The output is always Float64. `count == 0` converts every point_data
 * array.
 * @return a new mesh (free with mio_mesh_free), or NULL on failure.
 */
MIO_API mio_mesh* mio_data_point_to_cell(const mio_mesh* mesh, const char* const* names,
                                         int64_t count, const char* suffix);

/**
 * Average cell_data onto the points: each point's value is the (optionally
 * measure-weighted) mean over the incident cells. A point touched by no cell
 * gets NaN. The output is always Float64.
 * @return a new mesh (free with mio_mesh_free), or NULL on failure.
 */
MIO_API mio_mesh* mio_data_cell_to_point(const mio_mesh* mesh, const char* const* names,
                                         int64_t count, mio_cell_point_weight weight,
                                         const char* suffix);

/**
 * Evaluate an elementwise expression over the arrays at `location` and store
 * the result there as `output_name`.
 *
 * The grammar accepts + - * / unary minus, parentheses, numeric literals,
 * array names, and the functions abs, sqrt, min, max and norm. Nothing else is
 * evaluated -- there is no arbitrary-code path.
 * @param mesh        input mesh (unmodified).
 * @param expression  the expression text.
 * @param location    where operands are looked up and the result is stored.
 * @param output_name name of the new array (required).
 * @param overwrite   non-zero to allow replacing an existing array.
 * @return a new mesh (free with mio_mesh_free), or NULL on failure.
 */
MIO_API mio_mesh* mio_data_calc(const mio_mesh* mesh, const char* expression,
                                mio_data_location location, const char* output_name,
                                int overwrite);

/**
 * Condition the values of the selected data arrays (clamp / normalize /
 * standardize). For cell_data the statistics are computed jointly across all
 * cell blocks. `count == 0` conditions every array at `location`.
 * @param lo,hi           clamp bounds, or the normalize target range.
 * @param nan_replacement used when `nan_policy` is MIO_NAN_REPLACE.
 * @param suffix          NULL or "" replaces in place; otherwise the result is
 *                        stored as name+suffix.
 * @return a new mesh (free with mio_mesh_free), or NULL on failure.
 */
MIO_API mio_mesh* mio_data_condition(const mio_mesh* mesh, mio_data_location location,
                                     const char* const* names, int64_t count,
                                     mio_condition_mode mode, double lo, double hi,
                                     mio_condition_scope scope, mio_nan_policy nan_policy,
                                     double nan_replacement, const char* suffix);

/** Opaque per-array data summary. Destroy with mio_data_info_free(). */
typedef struct mio_data_info mio_data_info;

/** Fixed-size summary of one data array (see mio_data_info_entry). Per-component
 *  statistics are retrieved separately with mio_data_info_component(). */
typedef struct mio_data_array_info {
    int location;           /**< a mio_data_location */
    int dtype;              /**< a mio_dtype, as stored */
    int64_t num_blocks;     /**< cell_data: number of cell blocks; else 1 */
    int64_t num_entries;    /**< rows: points / cells over all blocks / length */
    int64_t num_components; /**< product of the trailing dimensions */
    int64_t num_values;     /**< num_entries * num_components */
    double min;             /**< over finite values; NaN when there are none */
    double max;             /**< over finite values; NaN when there are none */
    double mean;            /**< over finite values; NaN when there are none */
    int64_t num_nan;        /**< count of NaN values */
    int64_t num_inf;        /**< count of +/-inf values */
    int64_t num_finite;     /**< count of finite values */
    int inconsistent_blocks; /**< cell_data blocks disagree in component count */
} mio_data_array_info;

/**
 * Summarize every data array a mesh carries (read-only).
 * @return a result handle (free with mio_data_info_free), or NULL on failure.
 */
MIO_API mio_data_info* mio_data_info_create(const mio_mesh* mesh);

/** @return the number of arrays described, or -1 on error. */
MIO_API int64_t mio_data_info_count(const mio_data_info* info);

/**
 * Copy the `index`-th array's name into `buf`.
 * @return the required length excluding the NUL (so a value >= buflen means the
 *         name was truncated), or -1 on error.
 */
MIO_API int64_t mio_data_info_name(const mio_data_info* info, int64_t index, char* buf,
                                   int64_t buflen);

/** Fill `out` with the `index`-th array's summary. */
MIO_API mio_status mio_data_info_entry(const mio_data_info* info, int64_t index,
                                       mio_data_array_info* out);

/**
 * Per-component statistics of the `index`-th array. Any out pointer may be NULL.
 * @param comp component index, in [0, num_components).
 */
MIO_API mio_status mio_data_info_component(const mio_data_info* info, int64_t index, int64_t comp,
                                           double* min, double* max, double* mean);

/** Destroy a summary handle. Safe to call with NULL. */
MIO_API void mio_data_info_free(mio_data_info* info);

/**
 * Renumber a mesh (RCM / Morton / Hilbert) as a pure permutation of node and
 * element indices, to reduce sparse-matrix bandwidth / improve cache locality.
 * @param mesh   input mesh.
 * @param method "rcm", "morton", or "hilbert".
 * @return a result handle (free with mio_reorder_result_free), or NULL on
 *         failure (see mio_last_error()).
 */
MIO_API mio_reorder_result* mio_reorder(const mio_mesh* mesh, const char* method);

/**
 * Borrow the renumbered mesh from a reorder result. The returned handle is
 * owned by the result: valid until mio_reorder_result_free(result); do NOT
 * pass it to mio_mesh_free().
 * @param result a reorder result.
 * @return the renumbered mesh, or NULL if `result` is NULL.
 */
MIO_API const mio_mesh* mio_reorder_result_mesh(const mio_reorder_result* result);

/**
 * Zero-copy borrow of the node permutation (int64, shape (num_points,),
 * old index -> new index). The pointer is valid until the result is freed.
 * Any out-param may be NULL.
 * @param result a reorder result.
 * @param data   receives a pointer into result-owned int64 memory.
 * @param dtype  receives MIO_INT64.
 * @param n      receives num_points.
 */
MIO_API mio_status mio_reorder_result_node_perm(const mio_reorder_result* result, const void** data,
                                                mio_dtype* dtype, int64_t* n);

/**
 * Number of per-block cell permutations in a reorder result (equal to the
 * number of cell blocks of the mesh).
 * @param result a reorder result.
 * @return the block count, or -1 on error (see mio_last_error()).
 */
MIO_API int64_t mio_reorder_result_num_cell_perms(const mio_reorder_result* result);

/**
 * Zero-copy borrow of cell block `block`'s cell permutation (int64, shape
 * (num_cells_in_block,), old index -> new index). Any out-param may be NULL.
 * @param result a reorder result.
 * @param block  cell-block index in [0, mio_reorder_result_num_cell_perms).
 * @param data   receives a pointer into result-owned int64 memory.
 * @param dtype  receives MIO_INT64.
 * @param n      receives the block's cell count.
 */
MIO_API mio_status mio_reorder_result_cell_perm(const mio_reorder_result* result, int64_t block,
                                                const void** data, mio_dtype* dtype, int64_t* n);

/**
 * Transfer ownership of the renumbered mesh out of the result. Afterwards the
 * result no longer owns a mesh (its borrowed accessor yields an empty mesh).
 * @param result a reorder result.
 * @return a new owning mesh handle (free it with mio_mesh_free), or NULL on error.
 */
MIO_API mio_mesh* mio_reorder_result_take_mesh(mio_reorder_result* result);

/** Free a reorder result (and its owned mesh and permutation arrays). */
MIO_API void mio_reorder_result_free(mio_reorder_result* result);

/**
 * Connectivity bandwidth: the maximum over all cells of (max node index - min
 * node index) within a cell. A before/after measure of a reordering.
 * @param mesh input mesh.
 * @return the bandwidth (0 for a node-less mesh), or -1 on error.
 */
MIO_API int64_t mio_compute_bandwidth(const mio_mesh* mesh);

/**
 * Compare two meshes and produce a structured report. Compares points, cell
 * blocks (types + counts + connectivity), and point/cell/field data within the
 * tolerance abs_err <= atol + rtol*|expected| (expected taken from mesh `a`).
 * Named point_sets/cell_sets are NOT compared (they never reach the C++ core).
 * @param a         first ("expected") mesh.
 * @param b         second mesh.
 * @param atol      absolute tolerance.
 * @param rtol      relative tolerance.
 * @param unordered nonzero to match points by spatial proximity (tolerant to a
 *                  shuffled node order) instead of by index.
 * @param out       receives a new result handle (free with mio_diff_result_free).
 * @return MIO_OK, or an error status (see mio_last_error()).
 */
MIO_API mio_status mio_diff(const mio_mesh* a, const mio_mesh* b, double atol, double rtol,
                            int unordered, mio_diff_result** out);

/**
 * Convenience boolean: are two meshes equal within tolerance? (verdict is not
 * MIO_DIFF_DIFFERENT). Named sets are not considered.
 * @param a,b       the meshes to compare.
 * @param atol,rtol tolerances.
 * @param unordered nonzero for proximity node matching.
 * @param out_equal receives 1 (equal) or 0 (different).
 * @return MIO_OK, or an error status (see mio_last_error()).
 */
MIO_API mio_status mio_meshes_equal(const mio_mesh* a, const mio_mesh* b, double atol, double rtol,
                                    int unordered, int* out_equal);

/** Overall verdict of a diff result (MIO_DIFF_DIFFERENT on a NULL result). */
MIO_API mio_diff_verdict mio_diff_result_verdict(const mio_diff_result* result);

/**
 * Points comparison summary from a diff result. Any out-param may be NULL.
 * @param result        a diff result.
 * @param max_abs_error largest absolute coordinate error.
 * @param max_rel_error largest relative coordinate error.
 * @param worst_index   flat index of the worst-offending coordinate, or -1.
 * @param num_exceeding number of coordinates beyond tolerance.
 * @return MIO_OK, or an error status.
 */
MIO_API mio_status mio_diff_result_point_summary(const mio_diff_result* result,
                                                 double* max_abs_error, double* max_rel_error,
                                                 int64_t* worst_index, int64_t* num_exceeding);

/** Number of per-block comparisons in a diff result (-1 on a NULL result). */
MIO_API int64_t mio_diff_result_num_block_diffs(const mio_diff_result* result);

/**
 * One block comparison from a diff result. Any out-param may be NULL.
 * @param result             a diff result.
 * @param block              block index in [0, mio_diff_result_num_block_diffs).
 * @param type_mismatch      set to 1 if the two blocks have different cell types.
 * @param count_mismatch     set to 1 if they have different cell counts.
 * @param conn_mismatch_count number of cells whose connectivity differs.
 * @return MIO_OK, MIO_ERR_NOT_FOUND if `block` is out of range, else an error.
 */
MIO_API mio_status mio_diff_result_block(const mio_diff_result* result, int64_t block,
                                         int* type_mismatch, int* count_mismatch,
                                         int64_t* conn_mismatch_count);

/** Free a diff result. */
MIO_API void mio_diff_result_free(mio_diff_result* result);

/* ---------------------------------------------------------------------
 * Building a mesh (setters -- all COPY caller memory)
 * --------------------------------------------------------------------- */

/**
 * Assign the point coordinates, replacing any previous ones.
 * @param dtype      MIO_FLOAT32 or MIO_FLOAT64; stored as given.
 * @param num_points number of points.
 * @param dim        coordinates per point (usually 2 or 3).
 * @param xyz        row-major `(num_points, dim)` buffer.
 */
MIO_API mio_status mio_mesh_set_points(mio_mesh* mesh, mio_dtype dtype, int64_t num_points,
                                       int64_t dim, const void* xyz);

/**
 * Append one homogeneous cell block.
 * @param cell_type      meshio type name (e.g. "triangle", "tetra10").
 * @param num_cells      number of cells in the block.
 * @param nodes_per_cell nodes per cell (must match the type's fixed count
 *                       when it has one).
 * @param dtype          MIO_INT32 or MIO_INT64; widened to int64 internally
 *                       (the core's connectivity type).
 * @param connectivity   row-major `(num_cells, nodes_per_cell)` buffer of
 *                       0-based point indices.
 */
MIO_API mio_status mio_mesh_add_cell_block(mio_mesh* mesh, const char* cell_type,
                                           int64_t num_cells, int64_t nodes_per_cell,
                                           mio_dtype dtype, const void* connectivity);

/**
 * Attach a named per-point data array. `shape[0]` must equal the number of
 * points (e.g. `{num_points}` for a scalar field, `{num_points, 3}` for a
 * vector field).
 */
MIO_API mio_status mio_mesh_add_point_data(mio_mesh* mesh, const char* name, mio_dtype dtype,
                                           int32_t ndim, const int64_t* shape, const void* data);

/**
 * Append one per-cell-block array to the named cell-data field: call once per
 * cell block, in block order, after adding the blocks (`shape[0]` must equal
 * that block's cell count).
 */
MIO_API mio_status mio_mesh_append_cell_data(mio_mesh* mesh, const char* name, mio_dtype dtype,
                                             int32_t ndim, const int64_t* shape, const void* data);

/** Attach a named mesh-level (field) data array of arbitrary shape. */
MIO_API mio_status mio_mesh_add_field_data(mio_mesh* mesh, const char* name, mio_dtype dtype,
                                           int32_t ndim, const int64_t* shape, const void* data);

/* ---------------------------------------------------------------------
 * Reading a mesh back (getters -- zero-copy, see rule 3)
 * --------------------------------------------------------------------- */

/** @return the number of points, or -1 if `mesh` is NULL. */
MIO_API int64_t mio_mesh_num_points(const mio_mesh* mesh);

/** @return coordinates per point (usually 2 or 3), or -1 if `mesh` is NULL. */
MIO_API int64_t mio_mesh_point_dim(const mio_mesh* mesh);

/** Borrow the `(num_points, point_dim)` row-major coordinate buffer. */
MIO_API mio_status mio_mesh_get_points(const mio_mesh* mesh, const void** data, mio_dtype* dtype);

/** @return the number of cell blocks, or -1 if `mesh` is NULL. */
MIO_API int64_t mio_mesh_num_cell_blocks(const mio_mesh* mesh);

/**
 * Describe cell block `block`. Any out-param may be NULL.
 * A ragged block (`is_ragged == 1`: polygons/polyhedra of varying size)
 * reports `nodes_per_cell == 0` and its connectivity is not accessible
 * through this API (v1 limitation).
 */
MIO_API mio_status mio_mesh_cell_block_info(const mio_mesh* mesh, int64_t block,
                                            int64_t* num_cells, int64_t* nodes_per_cell,
                                            int32_t* is_ragged);

/** Copy cell block `block`'s meshio type name into `buf` (string rule 5). */
MIO_API int64_t mio_mesh_cell_block_type(const mio_mesh* mesh, int64_t block, char* buf,
                                         int64_t buflen);

/** Borrow cell block `block`'s row-major `(num_cells, nodes_per_cell)`
 *  0-based connectivity. Fails with MIO_ERR_UNSUPPORTED on a ragged block. */
MIO_API mio_status mio_mesh_cell_block_conn(const mio_mesh* mesh, int64_t block, const void** conn,
                                            mio_dtype* dtype);

/** @return the number of named point-data arrays, or -1 if `mesh` is NULL. */
MIO_API int64_t mio_mesh_num_point_data(const mio_mesh* mesh);

/** Copy the `index`-th point-data name (ascending lexicographic order --
 *  identical on every mesh backend) into `buf` (string rule 5). */
MIO_API int64_t mio_mesh_point_data_name(const mio_mesh* mesh, int64_t index, char* buf,
                                         int64_t buflen);

/** Borrow the named point-data array. `shape` (when non-NULL) must hold
 *  MIO_MAX_NDIM elements; any out-param may be NULL. */
MIO_API mio_status mio_mesh_get_point_data(const mio_mesh* mesh, const char* name,
                                           const void** data, mio_dtype* dtype, int32_t* ndim,
                                           int64_t* shape);

/** @return the number of named cell-data fields, or -1 if `mesh` is NULL. */
MIO_API int64_t mio_mesh_num_cell_data(const mio_mesh* mesh);

/** Copy the `index`-th cell-data name (sorted, as above) into `buf`. */
MIO_API int64_t mio_mesh_cell_data_name(const mio_mesh* mesh, int64_t index, char* buf,
                                        int64_t buflen);

/** @return how many per-block arrays the named cell-data field has (normally
 *          the number of cell blocks), or -1 on error. */
MIO_API int64_t mio_mesh_cell_data_num_blocks(const mio_mesh* mesh, const char* name);

/** Borrow the named cell-data field's array for cell block `block`. */
MIO_API mio_status mio_mesh_get_cell_data(const mio_mesh* mesh, const char* name, int64_t block,
                                          const void** data, mio_dtype* dtype, int32_t* ndim,
                                          int64_t* shape);

/** @return the number of named field-data arrays, or -1 if `mesh` is NULL. */
MIO_API int64_t mio_mesh_num_field_data(const mio_mesh* mesh);

/** Copy the `index`-th field-data name (sorted, as above) into `buf`. */
MIO_API int64_t mio_mesh_field_data_name(const mio_mesh* mesh, int64_t index, char* buf,
                                         int64_t buflen);

/** Borrow the named field-data array. */
MIO_API mio_status mio_mesh_get_field_data(const mio_mesh* mesh, const char* name,
                                           const void** data, mio_dtype* dtype, int32_t* ndim,
                                           int64_t* shape);

/* ---------------------------------------------------------------------
 * Named regions
 *
 * A region is a named group of mesh entities: a gmsh physical group, an Exodus
 * block / node set / side set, an Abaqus *NSET / *ELSET / *SURFACE, a MED
 * family, a Kratos SubModelPart. Before meshio++ 8.1 these never left the
 * Python layer at all -- this section is what closed that gap. See
 * doc/regions.md.
 *
 * Entry conventions, by kind:
 *   MIO_REGION_POINT  point indices,                     stride 1
 *   MIO_REGION_CELL   GLOBAL cell indices, block-major   stride 1
 *                     (block 0's cells first, then block 1's, ...)
 *   MIO_REGION_SIDE   (global cell index, local facet)   stride 2
 *
 * Facets are numbered as meshio++ numbers the faces of a 3-D cell and the edges
 * of a 2-D one. Entries are stored ascending and de-duplicated.
 * --------------------------------------------------------------------- */

/** Opaque snapshot of a mesh's regions. Destroy with mio_regions_free(). */
typedef struct mio_regions mio_regions;

/**
 * Snapshot every named region a mesh carries (read-only).
 * @return a handle (free with mio_regions_free), or NULL on failure.
 */
MIO_API mio_regions* mio_regions_create(const mio_mesh* mesh);

/** @return the number of regions described, or -1 on error. */
MIO_API int64_t mio_regions_count(const mio_regions* regions);

/**
 * Copy the `index`-th region's name into `buf`.
 * @return the required length excluding the NUL (so a value >= buflen means the
 *         name was truncated), or -1 on error.
 */
MIO_API int64_t mio_regions_name(const mio_regions* regions, int64_t index, char* buf,
                                 int64_t buflen);

/** Fill `out` with the `index`-th region's kind / dim / tag / sizes. */
MIO_API mio_status mio_regions_info(const mio_regions* regions, int64_t index,
                                    mio_region_info* out);

/**
 * Borrow the `index`-th region's entries: `num_entries * stride` int64 values,
 * row-major. Valid until mio_regions_free().
 * @param count receives `num_entries * stride` (may be NULL).
 * @return the entry buffer, or NULL on error / when the region is empty.
 */
MIO_API const int64_t* mio_regions_entries(const mio_regions* regions, int64_t index,
                                           int64_t* count);

/** Destroy a regions handle. Safe to call with NULL. */
MIO_API void mio_regions_free(mio_regions* regions);

/**
 * Add a region to a mesh, replacing any existing one with the same
 * (kind, name, dim, tag). Entries are copied and canonicalized (sorted,
 * de-duplicated).
 * @param entries `count` int64 values: `count` indices for POINT/CELL, or
 *                `count / 2` (cell, facet) pairs for SIDE.
 */
MIO_API mio_status mio_mesh_add_region(mio_mesh* mesh, const char* name, mio_region_kind kind,
                                       int32_t dim, int64_t tag, const int64_t* entries,
                                       int64_t count);

/* ---------------------------------------------------------------------
 * Transient (time-series) XDMF writing
 *
 * The write half of what mio_read_opts.time_step / mio_read_metadata_time_value
 * expose on the read side. A solver writes the mesh ONCE and then one cheap step
 * per solve, so this is a stateful handle rather than a (path, mesh) call -- the
 * one writer here that mio_write() cannot express.
 *
 *   mio_xdmf_series* s = mio_xdmf_series_create("out.xdmf", "HDF", -1);
 *   mio_xdmf_series_write_points_cells(s, mesh);
 *   for (k = 0; k < nsteps; ++k) { solve(); mio_xdmf_series_write_data(s, k*dt, mesh); }
 *   mio_xdmf_series_free(s);          // finalizes if mio_xdmf_series_finalize wasn't called
 *
 * The .xdmf light data is buffered and written once, at finalize; a series is
 * therefore only readable after the handle is finalized or freed. Heavy data for
 * "HDF" goes to a <path minus extension>.h5 SIBLING of the .xdmf.
 * --------------------------------------------------------------------- */

/** Opaque transient XDMF writer. Destroy with mio_xdmf_series_free(). */
typedef struct mio_xdmf_series mio_xdmf_series;

/**
 * Open a transient XDMF series for writing. Nothing is written yet.
 * @param path        the .xdmf/.xmf light-data file to write.
 * @param data_format "HDF" (NULL means "HDF"; needs an HDF5-enabled build),
 *                    "XML" or "Binary".
 * @param gzip_level  gzip level for "HDF" datasets, negative for none.
 * @return a handle (free with mio_xdmf_series_free), or NULL on failure -- an
 *         unknown data_format, or "HDF" without HDF5 support.
 */
MIO_API mio_xdmf_series* mio_xdmf_series_create(const char* path, const char* data_format,
                                                int32_t gzip_level);

/** How a series relates to a file already at its path (see mio_xdmf_series_opts). */
typedef enum mio_xdmf_series_mode {
    MIO_XDMF_SERIES_TRUNCATE = 0, /**< start a fresh series (the default) */
    MIO_XDMF_SERIES_APPEND = 1    /**< continue the existing one, if any */
} mio_xdmf_series_mode;

/**
 * Options for mio_xdmf_series_create_ex().
 *
 * ABI NOTE: as for mio_read_opts/mio_write_opts, new fields may only be
 * appended, replacing `reserved` capacity; never reorder, resize or repurpose
 * an existing field. Always zero-initialize through mio_xdmf_series_opts_init().
 */
typedef struct mio_xdmf_series_opts {
    /** "HDF" (NULL means "HDF"), "XML" or "Binary". */
    const char* data_format;
    int32_t gzip_level; /**< gzip level for "HDF" datasets, negative for none */
    int32_t mode;       /**< a mio_xdmf_series_mode */
    /** Nonzero flushes the light data after every write_data. Off by default:
     *  a flush re-serializes the whole document, so per-step flushing is
     *  quadratic in the step count. Drive mio_xdmf_series_flush() yourself for
     *  any other cadence. */
    int32_t auto_flush;
    int32_t reserved_pad; /**< must be zero; keeps the int64 tail aligned */
    int64_t reserved[6];  /**< must be zero; room for additive growth */
} mio_xdmf_series_opts;

/** Initialize `opts` to the defaults (HDF, no gzip, truncate, no auto-flush). */
MIO_API void mio_xdmf_series_opts_init(mio_xdmf_series_opts* opts);

/**
 * Open a transient XDMF series honouring `opts`. mio_xdmf_series_create() is
 * exactly mio_xdmf_series_create_ex(path, <defaults with that data_format and
 * gzip_level>) and is unchanged.
 * @return a handle (free with mio_xdmf_series_free), or NULL on failure.
 */
MIO_API mio_xdmf_series* mio_xdmf_series_create_ex(const char* path,
                                                   const mio_xdmf_series_opts* opts);

/**
 * Write the static grid every step shares. Call once, before the first
 * mio_xdmf_series_write_data(). Only the mesh's points and cells are used.
 */
MIO_API mio_status mio_xdmf_series_write_points_cells(mio_xdmf_series* series,
                                                      const mio_mesh* mesh);

/**
 * Write one time step's point_data and cell_data; the mesh's geometry is
 * ignored, and its cell blocks must match those of the static grid.
 */
MIO_API mio_status mio_xdmf_series_write_data(mio_xdmf_series* series, double time,
                                              const mio_mesh* mesh);

/**
 * Write the .xdmf and close the heavy-data container. Idempotent, and done by
 * mio_xdmf_series_free() too -- call it explicitly to see a write failure.
 */
MIO_API mio_status mio_xdmf_series_finalize(mio_xdmf_series* series);

/** One solver array for mio_xdmf_series_write_data_arrays(). Nothing is owned:
 *  `values` is read during the call and may be freed afterwards. */
typedef struct mio_named_array {
    const char* name;
    int64_t num_components; /**< values per entity; 1 for a scalar field */
    const double* values;   /**< row-major, entities * num_components doubles */
    int64_t num_values;     /**< length of `values`, validated against the grid */
} mio_named_array;

/**
 * Write one step from raw solver arrays instead of a whole mesh -- the
 * granularity a solver has once mio_xdmf_series_write_points_cells() has fixed
 * the geometry. Arrays are emitted in the order given.
 * @param point_arrays nodal arrays (may be NULL when num_point_arrays is 0).
 * @param cell_arrays  cell arrays (may be NULL when num_cell_arrays is 0).
 * @return MIO_OK, or an error if an array's length does not match the grid.
 */
MIO_API mio_status mio_xdmf_series_write_data_arrays(mio_xdmf_series* series, double time,
                                                     const mio_named_array* point_arrays,
                                                     int64_t num_point_arrays,
                                                     const mio_named_array* cell_arrays,
                                                     int64_t num_cell_arrays);

/**
 * Write the .xdmf as it currently stands, without finalizing, so a run that is
 * killed or still going leaves a readable file covering every flushed step.
 * The document is written to a temp file and renamed, so a crash during a flush
 * cannot truncate the previous one. Safe to call repeatedly; a no-op once
 * finalized.
 */
MIO_API mio_status mio_xdmf_series_flush(mio_xdmf_series* series);

/** @return the number of steps written so far, or -1 on error. */
MIO_API int64_t mio_xdmf_series_num_steps(const mio_xdmf_series* series);

/** @return 1 if the series has been finalized, 0 if not, -1 on error. */
MIO_API int32_t mio_xdmf_series_finalized(const mio_xdmf_series* series);

/** Finalize (if needed) and destroy a series handle. Safe to call with NULL. */
MIO_API void mio_xdmf_series_free(mio_xdmf_series* series);

/* --------------------------------------------------------------------------
 * Settings pipeline (v9.11.0). A settings.json document describes read ->
 * operation chain -> write (PascalCase ops/keys; see doc/pipeline.md), and
 * these run it whole -- the flat ABI carries JSON text only, no step struct.
 *
 * Both need a build with the JSON parser (-DMESHIOPLUSPLUS_WITH_JSON=ON, the
 * default when the src/cpp/third_party/json submodule is checked out);
 * otherwise they return MIO_ERR_INTERNAL and mio_last_error() names the flag.
 * -------------------------------------------------------------------------- */

/**
 * Run the settings pipeline stored in the file at `settings_path`.
 * @return MIO_OK, or an error (schema errors name the offending op/key).
 */
MIO_API mio_status mio_pipeline_run_file(const char* settings_path);

/**
 * Run the settings pipeline given as JSON text.
 * @return MIO_OK, or an error (schema errors name the offending op/key).
 */
MIO_API mio_status mio_pipeline_run_json(const char* json_text);

/** @return 1 when this build carries the JSON pipeline parser, else 0. */
MIO_API int32_t mio_pipeline_has_json(void);

#ifdef __cplusplus
}
#endif

#endif /* MESHIOPLUSPLUS_H */
