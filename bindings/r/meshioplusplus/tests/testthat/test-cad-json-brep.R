test_that("A Kratos CAD JSON file keeps its B-rep in the format side channel", {
  path <- tempfile(fileext = ".cad.json")
  output <- tempfile(fileext = ".cad.json")
  on.exit(unlink(c(path, output)))
  writeLines(c(
    '{"breps": [{"brep_id": 1, "faces": [{"brep_id": 2, "surface": {',
    '"is_trimmed": true, "degrees": [1, 1], "knot_vectors": [[0, 1], [0, 1]],',
    '"control_points": [[1, [0, 0, 0, 1]], [2, [2, 0, 0, 1]], [3, [0, 2, 0, 1]], [4, [2, 2, 0, 1]]]},',
    '"boundary_loops": [{"loop_type": "outer", "trimming_curves": [',
    '{"trim_index": 0, "curve_direction": true, "parameter_curve": {"degree": 1,',
    '"knot_vector": [0, 0, 1, 1], "active_range": [0, 1],',
    '"control_points": [[[0, 0, 0, 1]], [[1, 0, 0, 1]]]}}]}]}],',
    '"edges": [{"brep_id": 3, "topology": [{"brep_id": 2, "trim_index": 0}]}]}]}'
  ), path)
  r <- mio_read_with_info(path)
  expect_equal(mio_format_info_format(r$info), "cad_json")
  b <- mio_brep_info(r$info)
  expect_length(b$faces, 1)
  f <- b$faces[[1]]
  expect_equal(f$id, 2)
  expect_equal(f$patch, 1)
  expect_true(f$trimmed)
  expect_true(f$reduced_knots)
  expect_equal(f$control_point_ids, c(1, 2, 3, 4))
  expect_equal(b$patches[[1]]$knots[[1]], c(0, 0, 1, 1))
  expect_equal(dim(b$patches[[1]]$control_points), c(3, 4))
  expect_equal(f$loops[[1]]$type, "outer")
  expect_equal(f$loops[[1]]$trims[[1]]$trim_index, 0)
  expect_null(b$edges[[1]]$curve)
  expect_equal(b$edges[[1]]$topology[[1]]$relative_direction, -1)
  expect_equal(b$breps[[1]]$edges, 1)
  mio_write_with_info(r$mesh, r$info, output)
  back <- mio_read_with_info(output)
  expect_equal(mio_brep_info(back$info)$patches, b$patches)
  mio_release(back$mesh)
  mio_format_info_release(back$info)
  expect_error(mio_nurbs_info(r$info))
  mio_release(r$mesh)
  mio_format_info_release(r$info)
})

test_that("The isogeometric layer refines a trimmed face and reports the inside test", {
  path <- tempfile(fileext = ".cad.json")
  on.exit(unlink(path))
  writeLines(paste0(
    '{"breps": [{"brep_id": 1, "faces": [{"brep_id": 2, "surface": {"is_trimmed": true, ',
    '"degrees": [1, 1], "knot_vectors": [[0, 2], [0, 2]], "control_points": ',
    '[[1, [0, 0, 0, 1]], [2, [2, 0, 0, 1]], [3, [0, 2, 0, 1]], [4, [2, 2, 0, 1]]]}, ',
    '"boundary_loops": [{"loop_type": "outer", "trimming_curves": [',
    '{"trim_index": 0, "curve_direction": true, "parameter_curve": {"degree": 1, ',
    '"knot_vector": [0, 0, 1, 1], "active_range": [0, 1], "control_points": [[[0, 0, 0, 1]], [[1.5, 0, 0, 1]]]}}, ',
    '{"trim_index": 1, "curve_direction": true, "parameter_curve": {"degree": 1, ',
    '"knot_vector": [0, 0, 1, 1], "active_range": [0, 1], "control_points": [[[1.5, 0, 0, 1]], [[0, 1.5, 0, 1]]]}}, ',
    '{"trim_index": 2, "curve_direction": true, "parameter_curve": {"degree": 1, ',
    '"knot_vector": [0, 0, 1, 1], "active_range": [0, 1], "control_points": [[[0, 1.5, 0, 1]], [[0, 0, 0, 1]]]}}]}]}], ',
    '"edges": []}]}'), path)
  r <- mio_read_with_info(path)
  fine <- mio_iga_refine(r$info, subdivisions = c(2, 3), elevate = c(0, 1))
  expect_equal(mio_brep_info(fine$info)$patches[[1]]$degrees, c(1, 2))
  uv <- matrix(c(0.2, 0.2, 1.8, 1.8), nrow = 2)
  ev <- mio_iga_evaluate(fine$info, 1, uv, derivatives = 1)
  expect_equal(ev$inside, c(TRUE, FALSE))
  expect_equal(abs(ev$normals[3, 1]), 1, tolerance = 1e-13)
  expect_equal(dim(ev$derivatives), c(3, 2, 2))
  mio_release(fine$mesh)
  mio_format_info_release(fine$info)
  mio_release(r$mesh)
  mio_format_info_release(r$info)
})
