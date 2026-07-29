test_that("build information is available", {
  expect_true(nzchar(mio_version()))
  expect_true(mio_mesh_backend() %in% c("meshio", "native", "kratos"))
  expect_true(mio_format_readable("vtu"))
  expect_true(mio_format_writable("vtu"))
  expect_false(mio_format_writable("openfoam")) # read-only format
  expect_false(mio_format_readable("nonexistent"))
  expect_equal(mio_cell_type_num_nodes("tetra10"), 10L)
  expect_equal(mio_cell_type_dimension("triangle"), 2L)
  expect_equal(mio_cell_type_num_nodes("polygon"), -1L) # variable-size
})

test_that("a mesh can be built and inspected", {
  m <- fixture()
  on.exit(mio_release(m))

  expect_equal(mio_num_points(m), 5)
  expect_equal(mio_point_dim(m), 3)
  expect_equal(mio_num_cell_blocks(m), 1)
  expect_equal(mio_num_cells(m), 2)
  expect_equal(mio_cell_block_type(m, 1), "tetra")
  expect_equal(mio_cell_block_types(m), "tetra")

  info <- mio_cell_block_info(m, 1)
  expect_equal(info$num_cells, 2)
  expect_equal(info$nodes_per_cell, 4)
  expect_false(info$is_ragged)

  # Names come back sorted, on every mesh backend.
  expect_equal(mio_point_data_names(m), c("displacement", "temperature"))
  expect_equal(mio_cell_data_names(m), "material")
  expect_equal(mio_cell_data_num_blocks(m, "material"), 1)
})

test_that("the column-major shape identity holds", {
  m <- fixture()
  on.exit(mio_release(m))

  # An R (dim x num_points) matrix is the SAME layout as the C API's row-major
  # (num_points, dim) -- never transposed. With this non-square 3x5 fixture a
  # transpose would show up immediately as a 5x3 result or scrambled values.
  p <- mio_points(m)
  expect_equal(dim(p), c(3L, 5L))
  expect_equal(p, POINTS)
  expect_equal(p[, 2], c(1.1, 0.2, 0.3)) # point 2's coordinates
  expect_equal(p[1, ], c(0.0, 1.1, 0.4, 0.6, 1.4)) # every point's x

  # A 3-component vector field is (components x num_points).
  d <- mio_point_data(m, "displacement")
  expect_equal(dim(d), c(3L, 5L))
  expect_equal(as.vector(d), as.vector(VEC))

  # A scalar field stays a plain vector.
  expect_equal(as.vector(mio_point_data(m, "temperature")), TEMPERATURE)
  expect_equal(as.vector(mio_cell_data(m, "material", 1)), c(7, 9))
})

test_that("mio_connectivity is 1-based and mio_connectivity_raw is 0-based", {
  m <- fixture()
  on.exit(mio_release(m))

  c1 <- mio_connectivity(m, 1)
  c0 <- mio_connectivity_raw(m, 1)
  expect_equal(dim(c1), c(4L, 2L))
  expect_equal(c1, CONN)
  # Exactly the documented relationship between the two accessors.
  expect_equal(c1, c0 + 1)
})

test_that("data arrays carry their stored dtype as an attribute", {
  m <- fixture()
  on.exit(mio_release(m))
  # R has no native int64, so values arrive as double; the dtype the core
  # actually stored is reported separately rather than silently lost.
  expect_equal(attr(mio_point_data(m, "temperature"), "dtype"), "float64")
})

test_that("meshes round-trip through files", {
  m <- fixture()
  on.exit(mio_release(m))
  dir <- tempfile()
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  for (ext in c("vtu", "vtk")) {
    path <- file.path(dir, paste0("mesh.", ext))
    mio_write(m, path)
    expect_true(file.exists(path))

    r <- mio_read(path)
    expect_equal(mio_num_points(r), 5)
    expect_equal(mio_point_dim(r), 3)
    expect_equal(mio_cell_block_type(r, 1), "tetra")
    expect_equal(mio_points(r), POINTS, tolerance = 1e-10)
    expect_equal(mio_connectivity(r, 1), CONN)
    expect_true(mio_meshes_equal(m, r, atol = 1e-10, rtol = 1e-10))
    mio_release(r)
  }

  src <- file.path(dir, "mesh.vtu")
  dst <- file.path(dir, "converted.vtk")
  mio_convert(src, dst)
  expect_true(file.exists(dst))

  # sniff_format reads the leading bytes, not the extension.
  expect_equal(mio_sniff_format(src), "vtu")
})

test_that("selective reads and metadata work", {
  m <- fixture()
  dir <- tempfile()
  dir.create(dir)
  on.exit({
    mio_release(m)
    unlink(dir, recursive = TRUE)
  })
  path <- file.path(dir, "mesh.vtu")
  mio_write(m, path)

  r <- mio_read(path, points_only = TRUE)
  expect_equal(mio_num_points(r), 5)
  expect_length(mio_point_data_names(r), 0)
  mio_release(r)

  # `arrays = character(0)` is NOT the same request as `arrays = NULL`:
  # the first means "no arrays", the second "every array".
  r2 <- mio_read(path, arrays = "temperature")
  expect_equal(mio_point_data_names(r2), "temperature")
  mio_release(r2)

  r3 <- mio_read(path, arrays = character(0))
  expect_length(mio_point_data_names(r3), 0)
  mio_release(r3)

  meta <- mio_read_metadata(path)
  expect_equal(meta$num_points, 5)
  expect_equal(meta$point_dim, 3)
  expect_equal(meta$num_cells, 2)
  expect_equal(meta$num_cell_blocks, 1)
  expect_equal(meta$cell_block_types, "tetra")
  expect_false(meta$cell_block_is_ragged[1])
  expect_true("temperature" %in% meta$point_data_names)
  expect_true(mio_reader_supports_options("vtu"))
})

test_that("regions round-trip in memory", {
  m <- fixture()
  on.exit(mio_release(m))

  mio_add_region(m, "inlet", "point", c(1, 3, 5))
  mio_add_region(m, "solid", "cell", c(1, 2), dim = 3L, tag = 17)
  # (cell, facet) pairs: column-major, so 2 rows x N columns.
  mio_add_region(m, "wall", "side", matrix(c(1, 0, 2, 2), nrow = 2))

  rs <- mio_regions(m)
  expect_length(rs, 3)
  by_name <- rs
  names(by_name) <- vapply(rs, function(r) r$name, character(1))

  inlet <- by_name[["inlet"]]
  expect_equal(inlet$kind, "point")
  expect_equal(as.vector(inlet$entries), c(1, 3, 5)) # 1-based, as given

  solid <- by_name[["solid"]]
  expect_equal(solid$kind, "cell")
  expect_equal(solid$dim, 3L)
  expect_equal(solid$tag, 17)
  expect_equal(as.vector(solid$entries), c(1, 2))

  wall <- by_name[["wall"]]
  expect_equal(wall$kind, "side")
  expect_equal(dim(wall$entries), c(2L, 2L))
  # Canonicalized by the core. The CELL row is 1-based here; the FACET row is
  # an ordinal within the cell type, not a mesh index, so it is NOT shifted --
  # the same rule the Fortran and Julia bindings follow.
  expect_equal(wall$entries[1, ], c(1, 2)) # cells, 1-based
  expect_equal(wall$entries[2, ], c(0, 2)) # facets, unshifted

  expect_error(mio_add_region(m, "bad", "point", c(0, 1)), "1-based")
  expect_error(mio_add_region(m, "bad", "side", c(1, 2)), "value")
})

test_that("regions survive a file round-trip", {
  dir <- tempfile()
  dir.create(dir)
  m <- mio_mesh()
  on.exit({
    mio_release(m)
    unlink(dir, recursive = TRUE)
  })
  mio_set_points(m, POINTS)
  mio_add_cell_block(m, "triangle", matrix(c(1, 2, 3, 2, 3, 4), nrow = 3))
  mio_add_region(m, "physical", "cell", c(1, 2), dim = 2L, tag = 5)

  # Abaqus is one of the two Phase-1 region formats. Its *ELSET has no notion
  # of a dimension or an integer tag, so those come back unspecified (-1); the
  # group itself is what survives.
  path <- file.path(dir, "mesh.inp")
  mio_write(m, path)
  r <- mio_read(path)
  names <- vapply(mio_regions(r), function(x) x$name, character(1))
  expect_true("physical" %in% names)
  got <- mio_regions(r)[[which(names == "physical")]]
  expect_equal(got$kind, "cell")
  expect_equal(as.vector(got$entries), c(1, 2))
  mio_release(r)
})

test_that("geometry operations work", {
  m <- fixture()
  on.exit(mio_release(m))

  surf <- mio_extract_surface(m)
  expect_gt(mio_num_points(surf), 0)
  expect_equal(mio_cell_block_type(surf, 1), "triangle")
  mio_release(surf)

  surf2 <- mio_extract_surface(m, record_parent_ids = TRUE)
  expect_true("surface:parent_cell" %in% mio_cell_data_names(surf2))
  mio_release(surf2)

  skin <- mio_extract_skin(m, linearize = TRUE)
  expect_gt(mio_num_cells(skin), 0)
  mio_release(skin)

  q <- mio_attach_quality(m)
  expect_true(any(startsWith(mio_cell_data_names(q), "quality:")))
  mio_release(q)

  expect_equal(mio_quality_counts(m)$num_cells, 2)

  st <- mio_stats(m)
  expect_equal(st$num_points, 5)
  expect_equal(st$num_cells, 2)
  expect_equal(st$bbox_min[1], 0.0)
  expect_equal(st$bbox_max[1], 1.4)
  expect_gte(mio_compute_bandwidth(m), 0)
})

test_that("transform moves points and leaves topology alone", {
  m <- fixture()
  on.exit(mio_release(m))
  M <- rbind(
    c(1, 0, 0, 2),
    c(0, 1, 0, 3),
    c(0, 0, 1, 4),
    c(0, 0, 0, 1)
  )
  t <- mio_transform(m, M)
  p <- mio_points(t)
  expect_equal(p[, 1], c(2, 3, 4)) # the origin moved by the translation
  expect_equal(p[, 2], POINTS[, 2] + c(2, 3, 4))
  expect_equal(mio_connectivity(t, 1), CONN) # topology untouched
  mio_release(t)
})

test_that("refine reports 1-based maps with a 0 sentinel", {
  m <- fixture()
  on.exit(mio_release(m))

  r <- mio_refine(m, levels = 1L)
  expect_equal(mio_num_cells(r$mesh), 16) # 2 tetra -> 8 children each
  # Refinement never prunes, so the point map is the identity -- 1-based here,
  # like every copied index array.
  expect_equal(r$point_map, as.numeric(1:5))
  expect_length(r$cell_maps, 1)
  expect_equal(r$cell_maps[[1]], c(1, 9)) # each parent's first child
  mio_release(r$mesh)

  # convert_cells CAN prune, and a pruned point is reported as 0 -- never a
  # valid 1-based index, so the sentinel stays unambiguous.
  cc <- mio_convert_cells(m, "simplexify")
  expect_equal(mio_num_cells(cc$mesh), 2) # a tetra is already a simplex
  expect_true(all(cc$point_map >= 0))
  mio_release(cc$mesh)
})

test_that("selective refine closes up conformingly", {
  m <- fixture()
  on.exit(mio_release(m))

  # `cells` is 1-based here, like every index array this binding takes, and is
  # shifted to the C API's 0-based numbering inside.
  sel <- mio_refine(m, cells = c(1), record_levels = TRUE)
  expect_gt(mio_num_cells(sel$mesh), mio_num_cells(m))
  expect_lt(mio_num_cells(sel$mesh), 16) # not the uniform 8-per-cell
  expect_true("refine:level" %in% mio_cell_data_names(sel$mesh))
  mio_release(sel$mesh)

  # Propagation is the always-works baseline: on a connected mesh it reaches
  # every cell, which is exactly the uniform refinement.
  prop <- mio_refine(m, cells = c(1), closure = "propagate")
  expect_equal(mio_num_cells(prop$mesh), 16)
  mio_release(prop$mesh)

  # At most one selector, and both enum names are validated.
  expect_error(mio_refine(m, cells = c(1), region = "anything"))
  expect_error(mio_refine(m, closure = "blue"))
  expect_error(mio_refine(m, where_array = "q", where_op = "~="))
})

test_that("split and partition pieces own their handles", {
  m <- fixture()
  on.exit(mio_release(m))

  pieces <- mio_split(m, by = "type")
  expect_length(pieces, 1)
  expect_equal(names(pieces), "tetra")
  # The piece must still work after the C result handle is gone -- which it
  # does, because every piece takes ownership rather than borrowing.
  expect_equal(mio_num_cells(pieces[[1]]), 2)
  expect_equal(mio_points(pieces[[1]]), POINTS)
  mio_release(pieces[[1]])

  parts <- mio_partition(m, 2L, method = "sfc")
  expect_length(parts, 2)
  expect_equal(vapply(parts, function(p) p$part_id, integer(1)), c(0L, 1L))
  total <- sum(vapply(parts, function(p) mio_num_cells(p$mesh), numeric(1)))
  expect_equal(total, 2) # partition of unity
  for (p in parts) {
    expect_length(p$cell_maps, 1) # block structure kept 1:1
    mio_release(p$mesh)
  }

  labels <- mio_partition_labels(m, 2L, method = "sfc")
  expect_length(labels, 2)
  # Part IDs, not indices: deliberately NOT shifted to 1-based.
  expect_equal(sort(labels), c(0, 1))
})

test_that("reorder returns 1-based permutations", {
  m <- fixture()
  on.exit(mio_release(m))
  r <- mio_reorder(m, "rcm")
  expect_equal(mio_num_points(r$mesh), 5)
  expect_equal(sort(r$node_perm), as.numeric(1:5))
  expect_length(r$cell_perms, 1)
  expect_equal(sort(r$cell_perms[[1]]), as.numeric(1:2))
  mio_release(r$mesh)
})

test_that("clean, smooth and crop work", {
  m <- fixture()
  on.exit(mio_release(m))

  c1 <- mio_clean(m)
  expect_equal(mio_num_points(c1$mesh), 5) # nothing to remove here
  expect_equal(c1$points_welded, 0)
  mio_release(c1$mesh)

  s <- mio_smooth(m, iterations = 2L)
  expect_equal(mio_num_points(s$mesh), 5) # a pure coordinate move
  expect_equal(mio_num_cells(s$mesh), 2)
  expect_gte(s$max_displacement, 0)
  mio_release(s$mesh)

  b <- mio_crop_bbox(m, c(-1, -1, -1), c(10, 10, 10))
  expect_equal(mio_num_cells(b), 2) # the box holds everything
  mio_release(b)

  b2 <- mio_crop_bbox(m, c(100, 100, 100), c(101, 101, 101))
  expect_equal(mio_num_cells(b2), 0)
  mio_release(b2)

  p <- mio_crop_plane(m, c(0, 0, 0), c(1, 0, 0), mode = "any")
  expect_gte(mio_num_cells(p), 0)
  mio_release(p)
})

test_that("slice and isosurface work", {
  m <- fixture()
  on.exit(mio_release(m))

  sl <- mio_slice(m, c(0.5, 0.5, 0.5), c(0, 0, 1))
  expect_gte(mio_num_points(sl), 0) # may be empty; never an error
  mio_release(sl)

  iso <- mio_isosurface(m, "temperature", 3)
  expect_gte(mio_num_points(iso), 0)
  if (mio_num_cells(iso) > 0) {
    expect_true("iso:value" %in% mio_cell_data_names(iso))
    expect_true("iso:index" %in% mio_cell_data_names(iso))
  }
  mio_release(iso)

  # A cell_data array has no level set: that must fail by name.
  expect_error(mio_isosurface(m, "material", 8))
})

test_that("merge and interpolate work", {
  a <- fixture()
  b <- fixture()
  on.exit({
    mio_release(a)
    mio_release(b)
  })

  merged <- mio_merge(list(a, b))
  expect_equal(mio_num_points(merged), 10)
  expect_equal(mio_num_cells(merged), 4)
  expect_true("source_mesh_id" %in% mio_cell_data_names(merged))
  mio_release(merged)

  merged2 <- mio_merge(list(a, b), weld = TRUE, atol = 1e-9, source_tag = FALSE)
  expect_equal(mio_num_points(merged2), 5) # the two copies weld together
  mio_release(merged2)

  # The target already carries "temperature", so the default on_conflict
  # ("error") must refuse rather than silently overwrite.
  target <- fixture()
  expect_error(mio_interpolate(a, target, arrays = "temperature"))
  it <- mio_interpolate(a, target, arrays = "temperature", on_conflict = "suffix")
  expect_true("temperature_interp" %in% mio_point_data_names(it))
  mio_release(it)
  mio_release(target)
})

test_that("diff reports structured results", {
  a <- fixture()
  b <- fixture()
  on.exit({
    mio_release(a)
    mio_release(b)
  })

  rep <- mio_diff(a, b)
  expect_equal(rep$verdict, "identical")
  expect_equal(rep$points$num_exceeding, 0)

  moved <- POINTS
  moved[1, 1] <- 5.0
  mio_set_points(b, moved)
  rep2 <- mio_diff(a, b)
  expect_equal(rep2$verdict, "different")
  expect_equal(rep2$points$max_abs_error, 5.0)
  expect_false(mio_meshes_equal(a, b))
  expect_true(mio_meshes_equal(a, b, atol = 10)) # inside a loose tolerance
})

test_that("data operations never touch geometry", {
  m <- fixture()
  on.exit(mio_release(m))
  before_points <- mio_points(m)
  before_conn <- mio_connectivity(m, 1)

  d <- mio_data_drop(m, "point", "temperature")
  expect_equal(mio_point_data_names(d), "displacement")
  expect_equal(mio_points(d), before_points) # geometry bit-identical
  expect_equal(mio_connectivity(d, 1), before_conn)
  mio_release(d)

  k <- mio_data_keep(m, "point", "temperature")
  expect_equal(mio_point_data_names(k), "temperature")
  mio_release(k)

  r <- mio_data_rename(m, "point", "temperature", "T")
  expect_true("T" %in% mio_point_data_names(r))
  expect_false("temperature" %in% mio_point_data_names(r))
  mio_release(r)

  p2c <- mio_data_point_to_cell(m, "temperature")
  expect_true("temperature" %in% mio_cell_data_names(p2c))
  mio_release(p2c)

  c2p <- mio_data_cell_to_point(m, "material")
  expect_true("material" %in% mio_point_data_names(c2p))
  mio_release(c2p)

  calc <- mio_data_calc(m, "temperature * 2.0", "doubled")
  expect_equal(as.vector(mio_point_data(calc, "doubled")), c(2, 4, 6, 8, 10))
  mio_release(calc)

  cond <- mio_data_condition(m, "point", "temperature", cond_mode = "normalize")
  t <- as.vector(mio_point_data(cond, "temperature"))
  expect_equal(min(t), 0)
  expect_equal(max(t), 1)
  mio_release(cond)

  info <- mio_data_info(m)
  expect_length(info, 3)
  names <- vapply(info, function(x) x$name, character(1))
  temp <- info[[which(names == "temperature")]]
  expect_equal(temp$location, "point")
  expect_equal(temp$num_entries, 5)
  expect_equal(temp$num_components, 1)
  expect_equal(temp$min, 1)
  expect_equal(temp$max, 5)
  expect_equal(temp$mean, 3)
  expect_equal(temp$num_nan, 0)
})

test_that("errors carry the C API's own message", {
  # Every failure becomes an R condition carrying mio_last_error(); a status
  # code never reaches R.
  err <- tryCatch(mio_read("/definitely/not/a/mesh/file.vtu"), error = function(e) e)
  expect_s3_class(err, "error")
  expect_true(nzchar(conditionMessage(err)))
  expect_equal(conditionMessage(err), mio_last_error())

  m <- fixture()
  expect_error(mio_point_data(m, "no_such_array"))
  expect_error(mio_cell_block_type(m, 99))
  expect_error(mio_data_rename(m, "point", "missing", "x"))

  # Using a released mesh is an error, not a crash.
  mio_release(m)
  expect_error(mio_num_points(m), "released")
  mio_release(m) # idempotent
  expect_false(mio_is_open(m))

  # A foreign external pointer is rejected rather than dereferenced.
  expect_error(mio_num_points(list()), "mio_mesh")
})

test_that("a transient XDMF series round-trips", {
  # The one writer mio_write() cannot express: the grid goes out once and each
  # step is appended. "XML" keeps everything in the single .xdmf, so this runs
  # against a library built without HDF5 too.
  dir <- tempfile("mio_series")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  path <- file.path(dir, "series.xdmf")

  m <- fixture()
  s <- mio_xdmf_series(path, data_format = "XML")
  expect_s3_class(s, "mio_xdmf_series")
  expect_true(mio_xdmf_series_is_open(s))
  expect_equal(mio_xdmf_series_num_steps(s), 0)

  mio_xdmf_series_write_points_cells(s, m)

  # temperature = t + node id, so no two steps share a value and a step
  # mix-up cannot pass.
  times <- c(0, 0.5, 1)
  step <- mio_mesh()
  mio_set_points(step, POINTS)
  mio_add_cell_block(step, "tetra", CONN)
  for (t in times) {
    mio_add_point_data(step, "temperature", t + seq_len(5))
    mio_xdmf_series_write_data(s, t, step)
  }
  expect_equal(mio_xdmf_series_num_steps(s), 3)
  expect_output(print(s), "3 step")

  # The .xdmf is buffered until finalize; nothing is readable before it.
  mio_xdmf_series_finalize(s)
  mio_xdmf_series_finalize(s) # idempotent
  mio_xdmf_series_release(s)
  expect_false(mio_xdmf_series_is_open(s))
  mio_xdmf_series_release(s) # idempotent
  expect_error(mio_xdmf_series_num_steps(s), "released")
  expect_output(print(s), "released")

  meta <- mio_read_metadata(path)
  expect_equal(meta$time_values, times)

  for (k in seq_along(times)) {
    back <- mio_read(path, time_step = k - 1L)
    expect_equal(mio_num_points(back), 5)
    expect_equal(as.vector(mio_point_data(back, "temperature")), times[k] + seq_len(5))
  }

  # An unknown data format is an error carrying the C API's own message.
  expect_error(mio_xdmf_series(file.path(dir, "bad.xdmf"), data_format = "NoSuchFormat"))

  # A mesh and a series are separately tagged: neither is accepted for the
  # other, and a foreign pointer is rejected rather than dereferenced.
  expect_error(mio_xdmf_series_num_steps(m), "mio_xdmf_series")
  other <- mio_xdmf_series(file.path(dir, "b.xdmf"), data_format = "XML")
  expect_error(mio_num_points(other), "mio_mesh")
  expect_error(mio_xdmf_series_num_steps(list()), "mio_xdmf_series")
  # Released before the tempdir goes away: a finalizer running afterwards
  # would try to write the .xdmf into a directory that no longer exists.
  mio_xdmf_series_release(other)
})
