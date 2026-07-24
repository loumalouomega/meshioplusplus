## R interface to the meshio++ C library.
##
## Two facts govern the whole file, and both are documented rather than hidden:
##
##  * COLUMN-MAJOR. R matrices are column-major and the C core row-major, so an
##    R (dim x num_points) matrix is the SAME layout as the C API's row-major
##    (num_points, dim). Nothing is transposed anywhere.
##  * COPY-ONLY. R vectors are R-managed, so the C API's zero-copy borrow does
##    not survive into R without ALTREP machinery that is out of scope. EVERY
##    accessor here copies. There is deliberately no `_ptr` accessor, unlike the
##    Julia and Fortran bindings; the 0-based reader is named
##    `mio_connectivity_raw` so nobody reads it as a borrow.
##
## R also has no native 64-bit integer type, so the int64 arrays the C API
## returns (connectivity, region entries, index maps, permutations) come back as
## `double`. That is exact to 2^53 -- far beyond any mesh -- and avoids a hard
## dependency on `bit64` for a theoretical case.

# --- small helpers ------------------------------------------------------------

.mio_location <- function(location) {
  switch(location,
    point = 0L, cell = 1L, field = 2L,
    stop("`location` must be \"point\", \"cell\" or \"field\", not \"", location, "\"")
  )
}

.mio_mode <- function(mode) {
  switch(mode,
    all = 0L, any = 1L,
    stop("`mode` must be \"all\" or \"any\", not \"", mode, "\"")
  )
}

.mio_mmap <- function(mmap) {
  switch(mmap,
    auto = 0L, on = 1L, off = 2L,
    stop("`mmap` must be \"auto\", \"on\" or \"off\", not \"", mmap, "\"")
  )
}

.mio_block <- function(block) {
  block <- as.numeric(block)
  if (length(block) != 1L || is.na(block) || block < 1) {
    stop("`block` must be a single 1-based cell-block index")
  }
  block - 1  # 1-based in R, 0-based at the ABI
}

# --- introspection ------------------------------------------------------------

#' meshio++ build information
#'
#' @return `mio_version()` and `mio_mesh_backend()` return a single string;
#'   `mio_last_error()` returns the message of the most recent failed call on
#'   this thread, or `""`.
#' @export
mio_version <- function() .Call(R_mio_version)

#' @rdname mio_version
#' @export
mio_mesh_backend <- function() .Call(R_mio_mesh_backend)

#' @rdname mio_version
#' @export
mio_last_error <- function() .Call(R_mio_last_error)

#' Format support in this build
#'
#' @param format A format name, e.g. `"gmsh"`, `"vtu"`, `"med"`.
#' @return A single logical: whether the format can be read / written by the
#'   installed library. Formats backed by optional dependencies (HDF5, netCDF)
#'   exist only in builds configured with them.
#' @export
mio_format_readable <- function(format) .Call(R_mio_format_readable, as.character(format))

#' @rdname mio_format_readable
#' @export
mio_format_writable <- function(format) .Call(R_mio_format_writable, as.character(format))

#' Guess a mesh file's format from its contents
#'
#' @param path Path to an existing, readable file.
#' @return The format name, or `""` when the leading bytes are not a confident
#'   match. Ambiguous magics are never guessed at.
#' @export
mio_sniff_format <- function(path) .Call(R_mio_sniff_format, as.character(path))

#' Cell-type metadata
#'
#' @param name A meshio cell-type name, e.g. `"triangle"`, `"tetra10"`.
#' @return `mio_cell_type_num_nodes()` gives the fixed nodes per cell, or `-1`
#'   for the variable-size types (polygon, polyhedron, VTK Lagrange) and for
#'   unknown names. `mio_cell_type_dimension()` gives the topological dimension
#'   0-3, or `-1` for an unknown name.
#' @export
mio_cell_type_num_nodes <- function(name) {
  .Call(R_mio_cell_type_num_nodes, as.character(name))
}

#' @rdname mio_cell_type_num_nodes
#' @export
mio_cell_type_dimension <- function(name) {
  .Call(R_mio_cell_type_dimension, as.character(name))
}

#' @rdname mio_read
#' @export
mio_reader_supports_options <- function(format) {
  .Call(R_mio_reader_supports_options, as.character(format))
}

# --- lifecycle ----------------------------------------------------------------

#' A meshio++ mesh handle
#'
#' `mio_mesh()` creates an empty mesh. The handle is an external pointer with a
#' registered finalizer, so it is released automatically when garbage-collected;
#' `mio_release()` frees it immediately and is idempotent.
#'
#' @param mesh A `mio_mesh` object.
#' @param x A `mio_mesh` object.
#' @param ... Ignored.
#' @return `mio_mesh()` returns a new `mio_mesh`. `mio_release()` returns
#'   `NULL` invisibly. `mio_is_open()` returns a single logical.
#'   `print.mio_mesh()` returns its argument invisibly.
#' @examples
#' m <- mio_mesh()
#' mio_is_open(m)
#' mio_release(m)
#' @export
mio_mesh <- function() .Call(R_mio_mesh_create)

#' @rdname mio_mesh
#' @export
mio_release <- function(mesh) invisible(.Call(R_mio_mesh_release, mesh))

#' @rdname mio_mesh
#' @export
mio_is_open <- function(mesh) .Call(R_mio_mesh_is_open, mesh)

#' @rdname mio_mesh
#' @export
print.mio_mesh <- function(x, ...) {
  if (!mio_is_open(x)) {
    cat("<mio_mesh: released>\n")
  } else {
    cat(sprintf(
      "<mio_mesh: %g points, %g cells in %g block(s)>\n",
      mio_num_points(x), mio_num_cells(x), mio_num_cell_blocks(x)
    ))
  }
  invisible(x)
}

# --- file I/O -----------------------------------------------------------------

#' Read a mesh file
#'
#' @param path Path to the file.
#' @param format Explicit format name, or `NULL` to infer it from the extension.
#' @param points_only Read geometry only, skipping every data array.
#' @param metadata_only Read the header/summary only.
#' @param arrays `NULL` keeps every data array; a character vector keeps only
#'   the named ones. `character(0)` means *no* arrays, which is deliberately a
#'   different request from `NULL`.
#' @param mmap One of `"auto"`, `"on"`, `"off"`.
#' @param time_step Which step of a multi-step file to materialize. `0` (the
#'   default) is the first, preserving the historical behaviour; negative counts
#'   from the end. Out of range is an error naming the available count, never a
#'   silent clamp. Honoured by formats carrying a time series (currently
#'   `exodus`); `mio_read_metadata()$time_values` reports how many there are.
#' @return A `mio_mesh` object.
#' @examples
#' \dontrun{
#' m <- mio_read("bracket.msh")
#' m <- mio_read("bracket.vtu", points_only = TRUE)
#' m <- mio_read("run.exo", time_step = -1) # the last step
#' }
#' @export
mio_read <- function(path, format = NULL, points_only = FALSE, metadata_only = FALSE,
                     arrays = NULL, mmap = "auto", time_step = 0) {
  .Call(
    R_mio_read, as.character(path), format, isTRUE(points_only),
    isTRUE(metadata_only), if (is.null(arrays)) NULL else as.character(arrays),
    .mio_mmap(mmap), as.integer(time_step)
  )
}

#' Write a mesh file
#'
#' @param mesh A `mio_mesh` object.
#' @param path Destination path.
#' @param format Explicit format name, or `NULL` to infer it from the extension.
#' @return `NULL`, invisibly.
#' @export
mio_write <- function(mesh, path, format = NULL) {
  invisible(.Call(R_mio_write, mesh, as.character(path), format))
}

#' Convert a mesh file without materializing it
#'
#' Reads `in_path` and immediately writes `out_path`, the CLI's `convert` verb.
#'
#' @param in_path,out_path Source and destination paths.
#' @param in_format,out_format Explicit format names, or `NULL` to infer them.
#' @return `NULL`, invisibly.
#' @export
mio_convert <- function(in_path, out_path, in_format = NULL, out_format = NULL) {
  invisible(.Call(
    R_mio_convert, as.character(in_path), in_format,
    as.character(out_path), out_format
  ))
}

#' Summarize a mesh file without loading its arrays
#'
#' @param path Path to the file.
#' @param format Explicit format name, or `NULL` to infer it.
#' @return A named list with `num_points`, `point_dim`, `num_cells`,
#'   `num_cell_blocks`, the per-block vectors `cell_block_types`,
#'   `cell_block_num_cells`, `cell_block_nodes_per_cell`,
#'   `cell_block_is_ragged`, the three data-name vectors, `bbox` (a 3x2 matrix,
#'   or `NULL` -- absent is the normal case for a native summary, since
#'   computing it would mean decoding the coordinates), `fell_back`
#'   (`TRUE` when the whole file had to be read because the format has no
#'   header-only path), `time_values` (the file's recorded time-series
#'   values; length 0 for a format with no time concept -- this is the count
#'   `mio_read(time_step = ...)` may name), and `regions` (a list of
#'   `name`/`kind`/`dim`/`tag`/`num_entries` -- the `mio_regions()` shape minus
#'   the entries themselves; empty on a native metadata path, since none of
#'   those formats currently map regions).
#' @export
mio_read_metadata <- function(path, format = NULL) {
  .Call(R_mio_read_metadata, as.character(path), format)
}

# --- inspecting ---------------------------------------------------------------

#' Mesh size and cell-block structure
#'
#' @param mesh A `mio_mesh` object.
#' @param block A **1-based** cell-block index.
#' @return `mio_num_points()`, `mio_point_dim()`, `mio_num_cells()` and
#'   `mio_num_cell_blocks()` return a single number. `mio_cell_block_info()`
#'   returns a list with `num_cells`, `nodes_per_cell` and `is_ragged`.
#'   `mio_cell_block_type()` returns the meshio type name;
#'   `mio_cell_block_types()` returns one per block.
#' @export
mio_num_points <- function(mesh) .Call(R_mio_num_points, mesh)

#' @rdname mio_num_points
#' @export
mio_point_dim <- function(mesh) .Call(R_mio_point_dim, mesh)

#' @rdname mio_num_points
#' @export
mio_num_cell_blocks <- function(mesh) .Call(R_mio_num_cell_blocks, mesh)

#' @rdname mio_num_points
#' @export
mio_num_cells <- function(mesh) {
  # The C API has no single total-cell accessor (a documented gap), so this is
  # summed over the blocks.
  n <- mio_num_cell_blocks(mesh)
  if (n < 1) {
    return(0)
  }
  sum(vapply(seq_len(n), function(i) mio_cell_block_info(mesh, i)$num_cells, numeric(1)))
}

#' @rdname mio_num_points
#' @export
mio_cell_block_info <- function(mesh, block) {
  .Call(R_mio_cell_block_info, mesh, .mio_block(block))
}

#' @rdname mio_num_points
#' @export
mio_cell_block_type <- function(mesh, block) {
  .Call(R_mio_cell_block_type, mesh, .mio_block(block))
}

#' @rdname mio_num_points
#' @export
mio_cell_block_types <- function(mesh) {
  n <- mio_num_cell_blocks(mesh)
  if (n < 1) {
    return(character(0))
  }
  vapply(seq_len(n), function(i) mio_cell_block_type(mesh, i), character(1))
}

#' Point coordinates and cell connectivity
#'
#' `mio_points()` returns a `dim x num_points` matrix. Because R matrices are
#' column-major and the C core is row-major, that is the *same layout* as the C
#' API's `(num_points, dim)` -- nothing is transposed.
#'
#' `mio_connectivity()` returns **1-based** node indices;
#' `mio_connectivity_raw()` returns the ABI's own **0-based** ones. Both copy:
#' unlike the Julia and Fortran bindings, R has no zero-copy borrow, so there is
#' deliberately no `_ptr` accessor to imply one.
#'
#' Values arrive as `double` because R has no native 64-bit integer type; that
#' is exact to 2^53.
#'
#' @param mesh A `mio_mesh` object.
#' @param block A **1-based** cell-block index.
#' @return A numeric matrix: `dim x num_points` for `mio_points()`,
#'   `nodes_per_cell x num_cells` for the connectivity accessors.
#' @export
mio_points <- function(mesh) .Call(R_mio_points, mesh)

#' @rdname mio_points
#' @export
mio_connectivity <- function(mesh, block = 1) {
  .Call(R_mio_connectivity, mesh, .mio_block(block))
}

#' @rdname mio_points
#' @export
mio_connectivity_raw <- function(mesh, block = 1) {
  .Call(R_mio_connectivity_raw, mesh, .mio_block(block))
}

#' Named data arrays
#'
#' A data array shaped `(rows, components...)` in the C API comes back with the
#' dimensions reversed, so a scalar per-point field is a plain vector and a
#' 3-component vector field is a `3 x num_points` matrix -- the same memory,
#' read column-major. The dtype it was stored as is attached as a `"dtype"`
#' attribute; the values themselves are always `double`.
#'
#' @param mesh A `mio_mesh` object.
#' @param name The array name.
#' @param block A **1-based** cell-block index.
#' @return The `*_names()` functions return a character vector in ascending
#'   lexicographic order (identical on every mesh backend).
#'   `mio_cell_data_num_blocks()` returns a single number. The accessors return
#'   a numeric vector or array.
#' @export
mio_point_data_names <- function(mesh) .Call(R_mio_point_data_names, mesh)

#' @rdname mio_point_data_names
#' @export
mio_cell_data_names <- function(mesh) .Call(R_mio_cell_data_names, mesh)

#' @rdname mio_point_data_names
#' @export
mio_field_data_names <- function(mesh) .Call(R_mio_field_data_names, mesh)

#' @rdname mio_point_data_names
#' @export
mio_cell_data_num_blocks <- function(mesh, name) {
  .Call(R_mio_cell_data_num_blocks, mesh, as.character(name))
}

#' @rdname mio_point_data_names
#' @export
mio_point_data <- function(mesh, name) .Call(R_mio_point_data, mesh, as.character(name))

#' @rdname mio_point_data_names
#' @export
mio_cell_data <- function(mesh, name, block = 1) {
  .Call(R_mio_cell_data, mesh, as.character(name), .mio_block(block))
}

#' @rdname mio_point_data_names
#' @export
mio_field_data <- function(mesh, name) .Call(R_mio_field_data, mesh, as.character(name))

# --- building -----------------------------------------------------------------

#' Build a mesh from raw arrays
#'
#' All of these **copy** the caller's memory, so the R objects can be reclaimed
#' as soon as the call returns.
#'
#' `mio_set_points()` takes a `dim x num_points` matrix and
#' `mio_add_cell_block()` a `nodes_per_cell x num_cells` matrix of **1-based**
#' node indices; the shift to the core's 0-based indexing happens in the copy.
#'
#' `mio_append_cell_data()` must be called once per cell block, in block order,
#' after the blocks have been added.
#'
#' @param mesh A `mio_mesh` object.
#' @param points A `dim x num_points` numeric matrix.
#' @param cell_type A meshio cell-type name, e.g. `"tetra"`.
#' @param conn A `nodes_per_cell x num_cells` matrix of 1-based node indices.
#' @param name The array name.
#' @param data A numeric vector or array shaped `(components..., rows)`.
#' @return `NULL`, invisibly.
#' @examples
#' m <- mio_mesh()
#' mio_set_points(m, matrix(c(0, 0, 0, 1, 0, 0, 0, 1, 0), nrow = 3))
#' mio_add_cell_block(m, "triangle", matrix(c(1, 2, 3), nrow = 3))
#' mio_num_cells(m)
#' mio_release(m)
#' @export
mio_set_points <- function(mesh, points) {
  invisible(.Call(R_mio_set_points, mesh, as.matrix(points)))
}

#' @rdname mio_set_points
#' @export
mio_add_cell_block <- function(mesh, cell_type, conn) {
  invisible(.Call(R_mio_add_cell_block, mesh, as.character(cell_type), as.matrix(conn)))
}

#' @rdname mio_set_points
#' @export
mio_add_point_data <- function(mesh, name, data) {
  invisible(.Call(R_mio_add_point_data, mesh, as.character(name), data))
}

#' @rdname mio_set_points
#' @export
mio_append_cell_data <- function(mesh, name, data) {
  invisible(.Call(R_mio_append_cell_data, mesh, as.character(name), data))
}

#' @rdname mio_set_points
#' @export
mio_add_field_data <- function(mesh, name, data) {
  invisible(.Call(R_mio_add_field_data, mesh, as.character(name), data))
}

# --- named regions ------------------------------------------------------------

#' Named regions
#'
#' A region is a named group of mesh entities: a Gmsh physical group, an Exodus
#' block or node/side set, an Abaqus `*NSET`/`*ELSET`/`*SURFACE`, a MED family,
#' a Kratos SubModelPart.
#'
#' Entries are **1-based** here, with one exception: for a `"side"` region the
#' second row is a *facet ordinal within the cell type*, not a mesh index, so it
#' is passed through unshifted. Cell indices are **global and block-major**
#' (block 1's cells first, then block 2's).
#'
#' @param mesh A `mio_mesh` object.
#' @param name The region name.
#' @param kind One of `"point"`, `"cell"`, `"side"`.
#' @param entries A numeric vector of indices, or a `stride x n` matrix
#'   (`stride` is 2 for `"side"`, 1 otherwise).
#' @param dim Topological dimension, or `-1` when unspecified.
#' @param tag Format-native integer id, or `-1` when there is none.
#' @return `mio_regions()` returns a list of lists with `name`, `kind`, `dim`,
#'   `tag` and an `entries` matrix. `mio_add_region()` returns `NULL`
#'   invisibly.
#' @export
mio_regions <- function(mesh) .Call(R_mio_regions, mesh)

#' @rdname mio_regions
#' @export
mio_add_region <- function(mesh, name, kind, entries, dim = -1L, tag = -1L) {
  invisible(.Call(
    R_mio_add_region, mesh, as.character(name), as.character(kind),
    entries, as.integer(dim), as.numeric(tag)
  ))
}

# --- operations ---------------------------------------------------------------

#' Mesh operations
#'
#' Computations on a mesh rather than file formats. Each returns a new mesh (or
#' a list containing one); the input is never modified.
#'
#' @param mesh,source,target,a,b `mio_mesh` objects.
#' @param meshes A list of `mio_mesh` objects.
#' @param record_parent_ids Attach a cell-data array naming each output cell's
#'   source cell.
#' @param record_ids Attach `crop:`/`partition:` original-id arrays.
#' @param linearize Emit corner nodes only (triangle/quad output).
#' @param matrix A 4x4 affine matrix; point `p` maps to `M %*% c(p, 1)`.
#' @param rotate_vector_data Also rotate vector (trailing dim 3) and tensor
#'   (trailing dim 9) point data by the transform's 3x3 linear block.
#' @param weld Fuse coincident points within `atol`.
#' @param atol Absolute coincidence tolerance.
#' @param remove_orphans Drop points referenced by no cell.
#' @param drop_degenerate Drop degenerate cells.
#' @param drop_duplicate_cells Drop exact-duplicate cells.
#' @param method For `mio_smooth()`, `"taubin"` or `"laplacian"`; for
#'   `mio_reorder()`, `"rcm"`, `"morton"` or `"hilbert"`; for
#'   `mio_partition()`, `"sfc"`, `"kahip"` or `"auto"`; for
#'   `mio_interpolate()`, `"nearest"` or `"barycentric"`.
#' @param iterations Smoothing iterations.
#' @param lambda Relaxation factor in (0, 1). **Negative means "this method's
#'   own default"**: 0.5 for laplacian, 0.33 for taubin.
#' @param mu Taubin's un-shrinking factor, which must satisfy `mu < -lambda < 0`.
#' @param fix_boundary Pin every node on a boundary facet.
#' @param preserve_features Pin boundary nodes whose incident facets meet above
#'   `feature_angle`.
#' @param preserve_boundary Pin boundary vertices during decimation.
#' @param feature_angle Feature angle in degrees (30 is the `vtkFeatureEdges`
#'   convention).
#' @param guard_inversion Reject any move that would flip a cell's signed
#'   measure.
#' @param lo,hi Bounding-box corners (3 numbers each).
#' @param point,normal,origin Plane definition (3 numbers each).
#' @param mode `"all"` (every node inside) or `"any"` (at least one).
#' @param array The point-data array to contour.
#' @param isovalues One or more level values.
#' @param component Component of a multi-component array; negative for the row
#'   magnitude.
#' @param source_tag Add an integer `source_mesh_id` cell-data array.
#' @param data_policy `"intersection"` (keep keys present in every input) or
#'   `"fill"` (union, NaN for missing rows).
#' @param arrays Source array names to transfer; `NULL` means every source
#'   point-data array.
#' @param extrapolate Give a target point outside the source domain the nearest
#'   source value instead of `default_value`.
#' @param default_value Fill value outside the source domain.
#' @param on_conflict `"error"`, `"overwrite"` or `"suffix"`.
#' @param rtol Relative tolerance.
#' @param unordered Match points by spatial proximity instead of by index.
#' @param by `"type"`, `"component"`, or `"region"`/`"tag"`.
#' @param tag_name For `by = "region"`, the integer cell-data array name.
#' @param convert_mode `"linearize"`, `"simplexify"` or `"elevate"`.
#' @param levels How many times to apply the refinement templates.
#' @param ratio Fraction of faces to KEEP, in (0, 1]; negative = unset.
#' @param target_faces Absolute face count to stop at; negative = unset.
#' @param max_error Quadric-error ceiling; negative = unset.
#' @param placement `"optimal"`, `"midpoint"` or `"endpoint"`.
#' @param nparts Number of partitions.
#' @param imbalance KaHIP only: allowed imbalance fraction in (0, 1).
#' @param partition_mode KaHIP only: `"fast"`, `"eco"` or `"strong"`.
#' @param seed KaHIP only: random seed.
#' @param ghost_layers Reserved; must be 0.
#' @param weights Name of a scalar numeric cell-data array of per-cell weights.
#' @return A `mio_mesh`, or a named list containing one plus the operation's
#'   maps and counters. `mio_quality_counts()`, `mio_stats()` and `mio_diff()`
#'   return named lists; `mio_meshes_equal()` a single logical;
#'   `mio_compute_bandwidth()` a single number; `mio_split()` and
#'   `mio_partition()` lists of pieces; `mio_partition_labels()` a numeric
#'   vector of part ids.
#' @examples
#' \dontrun{
#' m <- mio_read("bracket.msh")
#' surf <- mio_extract_surface(m)
#' mio_write(surf, "bracket_surface.vtu")
#' }
#' @export
mio_extract_surface <- function(mesh, record_parent_ids = FALSE) {
  .Call(R_mio_extract_surface, mesh, isTRUE(record_parent_ids))
}

#' @rdname mio_extract_surface
#' @export
mio_extract_skin <- function(mesh, linearize = FALSE) {
  .Call(R_mio_extract_skin, mesh, isTRUE(linearize))
}

#' @rdname mio_extract_surface
#' @export
mio_attach_quality <- function(mesh) .Call(R_mio_attach_quality, mesh)

#' @rdname mio_extract_surface
#' @export
mio_quality_counts <- function(mesh) .Call(R_mio_quality_counts, mesh)

#' @rdname mio_extract_surface
#' @export
mio_transform <- function(mesh, matrix, rotate_vector_data = FALSE) {
  .Call(R_mio_transform, mesh, as.numeric(matrix), isTRUE(rotate_vector_data))
}

#' @rdname mio_extract_surface
#' @export
mio_clean <- function(mesh, weld = FALSE, atol = 1e-12, remove_orphans = TRUE,
                      drop_degenerate = TRUE, drop_duplicate_cells = TRUE) {
  .Call(
    R_mio_clean, mesh, isTRUE(weld), as.numeric(atol), isTRUE(remove_orphans),
    isTRUE(drop_degenerate), isTRUE(drop_duplicate_cells)
  )
}

#' @rdname mio_extract_surface
#' @export
mio_smooth <- function(mesh, method = "taubin", iterations = 10L, lambda = -1,
                       mu = -0.34, fix_boundary = TRUE, preserve_features = TRUE,
                       feature_angle = 30, guard_inversion = TRUE) {
  .Call(
    R_mio_smooth, mesh, as.character(method), as.integer(iterations),
    as.numeric(lambda), as.numeric(mu), isTRUE(fix_boundary),
    isTRUE(preserve_features), as.numeric(feature_angle), isTRUE(guard_inversion)
  )
}

#' @rdname mio_extract_surface
#' @export
mio_crop_bbox <- function(mesh, lo, hi, mode = "all", record_ids = FALSE) {
  .Call(
    R_mio_crop_bbox, mesh, as.numeric(lo), as.numeric(hi), .mio_mode(mode),
    isTRUE(record_ids)
  )
}

#' @rdname mio_extract_surface
#' @export
mio_crop_plane <- function(mesh, point, normal, mode = "all", record_ids = FALSE) {
  .Call(
    R_mio_crop_plane, mesh, as.numeric(point), as.numeric(normal), .mio_mode(mode),
    isTRUE(record_ids)
  )
}

#' @rdname mio_extract_surface
#' @export
mio_slice <- function(mesh, origin, normal, record_parent_ids = FALSE) {
  .Call(
    R_mio_slice, mesh, as.numeric(origin), as.numeric(normal),
    isTRUE(record_parent_ids)
  )
}

#' @rdname mio_extract_surface
#' @export
mio_isosurface <- function(mesh, array, isovalues, component = -1L,
                           record_parent_ids = FALSE) {
  .Call(
    R_mio_isosurface, mesh, as.character(array), as.numeric(isovalues),
    as.integer(component), isTRUE(record_parent_ids)
  )
}

#' @rdname mio_extract_surface
#' @export
mio_merge <- function(meshes, weld = FALSE, atol = 1e-12, source_tag = TRUE,
                      data_policy = "intersection", drop_duplicate_cells = FALSE) {
  policy <- switch(data_policy,
    intersection = 0L, fill = 1L,
    stop("`data_policy` must be \"intersection\" or \"fill\"")
  )
  .Call(
    R_mio_merge, meshes, isTRUE(weld), as.numeric(atol), isTRUE(source_tag),
    policy, isTRUE(drop_duplicate_cells)
  )
}

#' @rdname mio_extract_surface
#' @export
mio_interpolate <- function(source, target, method = "nearest", arrays = NULL,
                            extrapolate = FALSE, default_value = 0,
                            on_conflict = "error") {
  .Call(
    R_mio_interpolate, source, target, as.character(method),
    if (is.null(arrays)) NULL else as.character(arrays), isTRUE(extrapolate),
    as.numeric(default_value), as.character(on_conflict)
  )
}

#' @rdname mio_extract_surface
#' @export
mio_meshes_equal <- function(a, b, atol = 0, rtol = 0, unordered = FALSE) {
  .Call(R_mio_meshes_equal, a, b, as.numeric(atol), as.numeric(rtol), isTRUE(unordered))
}

#' @rdname mio_extract_surface
#' @export
mio_diff <- function(a, b, atol = 0, rtol = 0, unordered = FALSE) {
  .Call(R_mio_diff, a, b, as.numeric(atol), as.numeric(rtol), isTRUE(unordered))
}

#' @rdname mio_extract_surface
#' @export
mio_stats <- function(mesh) .Call(R_mio_stats, mesh)

#' @rdname mio_extract_surface
#' @export
mio_compute_bandwidth <- function(mesh) .Call(R_mio_compute_bandwidth, mesh)

#' @rdname mio_extract_surface
#' @export
mio_reorder <- function(mesh, method = "rcm") {
  .Call(R_mio_reorder, mesh, as.character(method))
}

#' @rdname mio_extract_surface
#' @export
mio_split <- function(mesh, by = "type", tag_name = NULL) {
  .Call(R_mio_split, mesh, as.character(by), tag_name)
}

#' @rdname mio_extract_surface
#' @export
mio_convert_cells <- function(mesh, convert_mode, record_parent_ids = FALSE) {
  .Call(
    R_mio_convert_cells, mesh, as.character(convert_mode),
    isTRUE(record_parent_ids)
  )
}

#' @rdname mio_extract_surface
#' @export
mio_refine <- function(mesh, levels = 1L, record_parent_ids = FALSE) {
  .Call(R_mio_refine, mesh, as.integer(levels), isTRUE(record_parent_ids))
}

#' @rdname mio_extract_surface
#' @export
mio_decimate <- function(mesh, ratio = -1, target_faces = -1, max_error = -1,
                         placement = "optimal", preserve_boundary = TRUE,
                         preserve_features = TRUE, feature_angle = 30) {
  .Call(
    R_mio_decimate, mesh, as.numeric(ratio), as.numeric(target_faces),
    as.numeric(max_error), as.character(placement), isTRUE(preserve_boundary),
    isTRUE(preserve_features), as.numeric(feature_angle)
  )
}

#' @rdname mio_extract_surface
#' @export
mio_partition <- function(mesh, nparts, method = "auto", imbalance = 0.03,
                          partition_mode = "eco", seed = 0L, record_ids = FALSE,
                          ghost_layers = 0L, weights = NULL) {
  .Call(
    R_mio_partition, mesh, as.integer(nparts), as.character(method),
    as.numeric(imbalance), as.character(partition_mode), as.integer(seed),
    isTRUE(record_ids), as.integer(ghost_layers), weights
  )
}

#' @rdname mio_extract_surface
#' @export
mio_partition_labels <- function(mesh, nparts, method = "auto", imbalance = 0.03,
                                 partition_mode = "eco", seed = 0L, weights = NULL) {
  .Call(
    R_mio_partition_labels, mesh, as.integer(nparts), as.character(method),
    as.numeric(imbalance), as.character(partition_mode), as.integer(seed),
    weights, mio_num_cells(mesh)
  )
}

# --- data operations ----------------------------------------------------------

#' Data operations
#'
#' These act on a mesh's point, cell and field data arrays and **never touch
#' geometry**: points, connectivity and block order come through bit-identical.
#'
#' The combined `data_manage` of the C++ API is not exposed across the C ABI (a
#' documented gap); `mio_data_drop()`, `mio_data_keep()` and
#' `mio_data_rename()` compose to the same effect.
#'
#' @param mesh A `mio_mesh` object.
#' @param location One of `"point"`, `"cell"`, `"field"`.
#' @param names Array names. `NULL` or `character(0)` means "every array at that
#'   location" for the averaging and conditioning entry points, and "nothing"
#'   for `mio_data_drop()`.
#' @param ignore_missing Skip names that do not exist.
#' @param from,to Old and new array names.
#' @param suffix Write to `name + suffix` instead of in place.
#' @param weight `"uniform"` or `"measure"` (weight by cell measure).
#' @param expression An elementwise expression over the arrays at `location`.
#'   The grammar accepts `+ - * /`, unary minus, parentheses, numeric literals,
#'   array names and the functions `abs`, `sqrt`, `min`, `max` and `norm`.
#'   Nothing else is evaluated: there is no arbitrary-code path.
#' @param output_name Name of the array to create.
#' @param overwrite Allow replacing an existing array.
#' @param cond_mode `"clamp"`, `"normalize"` or `"standardize"`.
#' @param lo,hi Clamp bounds, or the normalize target range.
#' @param scope `"component"` (each trailing component independently) or
#'   `"magnitude"` (rescale whole rows, preserving direction).
#' @param nan_policy `"ignore"`, `"replace"` or `"fail"`. Non-finite values are
#'   always excluded from reductions regardless.
#' @param nan_replacement Used when `nan_policy` is `"replace"`.
#' @return A new `mio_mesh`, except `mio_data_info()`, which returns a list of
#'   per-array summaries.
#' @export
mio_data_drop <- function(mesh, location, names, ignore_missing = FALSE) {
  .Call(
    R_mio_data_drop, mesh, .mio_location(location),
    if (is.null(names)) NULL else as.character(names), isTRUE(ignore_missing)
  )
}

#' @rdname mio_data_drop
#' @export
mio_data_keep <- function(mesh, location, names, ignore_missing = FALSE) {
  .Call(
    R_mio_data_keep, mesh, .mio_location(location),
    if (is.null(names)) NULL else as.character(names), isTRUE(ignore_missing)
  )
}

#' @rdname mio_data_drop
#' @export
mio_data_rename <- function(mesh, location, from, to) {
  .Call(
    R_mio_data_rename, mesh, .mio_location(location), as.character(from),
    as.character(to)
  )
}

#' @rdname mio_data_drop
#' @export
mio_data_point_to_cell <- function(mesh, names = NULL, suffix = NULL) {
  .Call(
    R_mio_data_point_to_cell, mesh,
    if (is.null(names)) NULL else as.character(names), suffix
  )
}

#' @rdname mio_data_drop
#' @export
mio_data_cell_to_point <- function(mesh, names = NULL, weight = "uniform",
                                   suffix = NULL) {
  w <- switch(weight,
    uniform = 0L, measure = 1L,
    stop("`weight` must be \"uniform\" or \"measure\"")
  )
  .Call(
    R_mio_data_cell_to_point, mesh,
    if (is.null(names)) NULL else as.character(names), w, suffix
  )
}

#' @rdname mio_data_drop
#' @export
mio_data_calc <- function(mesh, expression, output_name, location = "point",
                          overwrite = FALSE) {
  .Call(
    R_mio_data_calc, mesh, as.character(expression), .mio_location(location),
    as.character(output_name), isTRUE(overwrite)
  )
}

#' @rdname mio_data_drop
#' @export
mio_data_condition <- function(mesh, location, names = NULL, cond_mode = "clamp",
                               lo = 0, hi = 1, scope = "component",
                               nan_policy = "ignore", nan_replacement = 0,
                               suffix = NULL) {
  md <- switch(cond_mode,
    clamp = 0L, normalize = 1L, standardize = 2L,
    stop("`cond_mode` must be \"clamp\", \"normalize\" or \"standardize\"")
  )
  sc <- switch(scope,
    component = 0L, magnitude = 1L,
    stop("`scope` must be \"component\" or \"magnitude\"")
  )
  np <- switch(nan_policy,
    ignore = 0L, replace = 1L, fail = 2L,
    stop("`nan_policy` must be \"ignore\", \"replace\" or \"fail\"")
  )
  .Call(
    R_mio_data_condition, mesh, .mio_location(location),
    if (is.null(names)) NULL else as.character(names), md, as.numeric(lo),
    as.numeric(hi), sc, np, as.numeric(nan_replacement), suffix
  )
}

#' @rdname mio_data_drop
#' @export
mio_data_info <- function(mesh) .Call(R_mio_data_info, mesh)
