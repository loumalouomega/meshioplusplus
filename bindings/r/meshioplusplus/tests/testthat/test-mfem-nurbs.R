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
  # the isogeometric layer: refine, write back, evaluate
  fine <- mio_iga_refine(r$info, subdivisions = 3, elevate = 1)
  expect_equal(mio_format_info_format(fine$info), "mfem")
  expect_equal(mio_nurbs_info(fine$info)$patches[[1]]$degrees, c(3, 2))
  refined_path <- tempfile(fileext = ".mesh")
  on.exit(unlink(refined_path), add = TRUE)
  mio_write_with_info(fine$mesh, fine$info, refined_path, format = "mfem")
  expect_equal(readLines(refined_path, n = 1), "MFEM NURBS mesh v1.0")
  uv <- matrix(c(0.3, 0.7, 0.9, 0.2), nrow = 2)
  a <- mio_iga_evaluate(r$info, 1, uv, derivatives = 1)
  b <- mio_iga_evaluate(fine$info, 1, uv)
  expect_equal(dim(a$points), c(2, 2))
  expect_lt(max(abs(a$points - b$points)), 1e-13)
  expect_equal(dim(a$derivatives), c(2, 2, 2))
  expect_null(a$second_derivatives)
  c2 <- mio_iga_evaluate(r$info, 1, uv, derivatives = 2)
  expect_equal(dim(c2$second_derivatives), c(2, 2, 2, 2))
  expect_equal(c2$second_derivatives[, 1, 2, ], c2$second_derivatives[, 2, 1, ])
  expect_null(a$normals)
  expect_null(a$inside)
  expect_length(a$fields, 0)
  expect_error(mio_iga_evaluate(r$info, 1, matrix(c(2, 0.5), nrow = 2)), "knot range")
  # closest point: the annulus is a surface in 2-D here
  cp <- mio_iga_closest_point(r$info, matrix(c(3, 0, 0, 1.5), nrow = 2))
  expect_equal(dim(cp$params), c(2, 2))
  expect_equal(cp$points[, 1], c(2, 0), tolerance = 1e-12)
  expect_equal(cp$distances, c(1, 0), tolerance = 1e-10)
  expect_error(mio_iga_closest_point(r$info, matrix(c(3, 0), nrow = 2), patch = 4), "patch")
  # knot removal undoes the refinement's insertions: two knots per direction
  red <- mio_iga_remove_knots(fine$info, direction = 1)
  expect_equal(red$report$num_removed, 2)
  expect_equal(red$report$removed, c(1L, 1L))
  expect_equal(red$report$knots, c(1 / 3, 2 / 3))
  expect_equal(mio_nurbs_info(red$info)$patches[[1]]$degrees, c(3, 2))
  mio_release(red$mesh)
  mio_format_info_release(red$info)
  expect_error(mio_iga_remove_knots(fine$info, knots = 0.5), "knot")
  expect_error(mio_iga_remove_knots(fine$info, direction = 4), "direction")
  # Bezier extraction: the annulus has unequal degrees; a bi-quadratic square is
  # one Bezier cell of nine points
  expect_error(mio_iga_bezier_extract(r$info, 1), "unequal degrees")
  square <- tempfile(fileext = ".mesh")
  on.exit(unlink(square), add = TRUE)
  writeLines(c(
    "MFEM NURBS mesh v1.0", "dimension", "2", "elements", "1", "5 3 0 1 2 3",
    "boundary", "0", "edges", "4", "0 0 1", "0 3 2", "1 0 3", "1 1 2",
    "vertices", "4", "knotvectors", "2", "2 3 0 0 0 1 1 1", "2 3 0 0 0 1 1 1",
    "weights", rep("1", 9),
    "FiniteElementSpace", "FiniteElementCollection: NURBS2", "VDim: 2", "Ordering: 1",
    "0 0", "1 0", "1 1", "0 1", "0.5 0", "0.5 1", "0 0.5", "1 0.5", "0.5 0.5"
  ), square)
  sq <- mio_read_with_info(square, format = "mfem")
  bez <- mio_iga_bezier_extract(sq$info, 1)
  expect_equal(mio_num_points(bez), 9)
  expect_error(mio_iga_bezier_extract(sq$info, 2), "out of range")
  mio_release(bez)
  mio_release(sq$mesh)
  mio_format_info_release(sq$info)
  sampled <- mio_iga_sample(r$info)
  expect_equal(mio_format_info_format(sampled$info), "mfem")
  mio_release(sampled$mesh)
  mio_format_info_release(sampled$info)
  mio_release(fine$mesh)
  mio_format_info_release(fine$info)
  expect_error(mio_gmsh_info(r$info))
  mio_release(r$mesh)
  mio_format_info_release(r$info)
})
