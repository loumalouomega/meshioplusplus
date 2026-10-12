test_that("An MFEM NURBS mesh keeps its patches in the format side channel", {
  path <- tempfile(fileext = ".mesh")
  output <- tempfile(fileext = ".mesh")
  on.exit(unlink(c(path, output)))
  w <- "0.70710678118654757"
  writeLines(c(
    "MFEM NURBS mesh v1.0", "dimension", "2", "elements", "1", "5 3 0 1 2 3",
    "boundary", "0", "edges", "4", "0 0 1", "0 3 2", "1 0 3", "1 1 2",
    "vertices", "4", "knotvectors", "2", "2 3 0 0 0 1 1 1", "1 2 0 0 1 1",
    "weights", "1", "1", "1", "1", w, w,
    "FiniteElementSpace", "FiniteElementCollection: NURBS2", "VDim: 2", "Ordering: 1",
    "1 0", "0 1", "0 2", "2 0", "1 1", "2 2"
  ), path)
  r <- mio_read_with_info(path, format = "mfem")
  expect_equal(mio_format_info_format(r$info), "mfem")
  copied <- mio_nurbs_info(r$info)
  expect_length(copied$patches, 1)
  p <- copied$patches[[1]]
  expect_equal(p$degrees, c(2, 1))
  expect_equal(p$knots[[2]], c(0, 0, 1, 1))
  expect_equal(dim(p$control_points), c(2, 6))
  expect_equal(p$control_points[, 2], c(1, 1))
  expect_equal(p$weights[2], sqrt(0.5))
  expect_equal(copied$field_names, character())
  expect_equal(copied$topology$num_vertices, 4)
  expect_equal(copied$topology$patch_vertices[, 1], c(1, 2, 3, 4))
  expect_equal(copied$topology$edges[, 1], c(0, 1, 2))
  expect_equal(copied$cell_patch, 1)
  mio_write_with_info(r$mesh, r$info, output, format = "mfem")
  expect_equal(readLines(output, n = 1), "MFEM NURBS mesh v1.0")
  back <- mio_read_with_info(output, format = "mfem")
  expect_equal(mio_nurbs_info(back$info)$patches, copied$patches)
  mio_release(back$mesh)
  mio_format_info_release(back$info)
  expect_error(mio_gmsh_info(r$info))
  mio_release(r$mesh)
  mio_format_info_release(r$info)
})
