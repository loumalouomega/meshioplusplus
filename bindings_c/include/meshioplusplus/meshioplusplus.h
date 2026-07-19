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

typedef enum mio_status {
    MIO_OK = 0,           /**< success */
    MIO_ERR_READ = 1,     /**< file could not be parsed / read-side failure */
    MIO_ERR_WRITE = 2,    /**< mesh could not be serialized / write-side failure */
    MIO_ERR_INVALID_ARG = 3, /**< NULL handle/pointer, bad dtype, bad shape, ... */
    MIO_ERR_NOT_FOUND = 4,   /**< unknown format/data name/block index */
    MIO_ERR_UNSUPPORTED = 5, /**< valid but unsupported (e.g. ragged connectivity) */
    MIO_ERR_INTERNAL = 6     /**< unexpected failure; see mio_last_error() */
} mio_status;

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
 * C++ `meshioplusplus::CellType` enum (cpp/include/meshioplusplus/
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

/* ---------------------------------------------------------------------
 * Version / build introspection
 * --------------------------------------------------------------------- */

/** @return the meshio++ version string, e.g. "6.1.0" (static storage). */
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

/** Write a mesh. `format` as in mio_read(). */
MIO_API mio_status mio_write(const char* path, const mio_mesh* mesh, const char* format);

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

#ifdef __cplusplus
}
#endif

#endif /* MESHIOPLUSPLUS_H */
