test_that("Gmsh periodic links use the format side channel", {
  path <- tempfile(fileext = ".msh")
  output <- tempfile(fileext = ".msh")
  on.exit(unlink(c(path, output)))
  writeLines(c(
    "$MeshFormat", "4.0 0 8", "$EndMeshFormat",
    "$Nodes", "1 3", "22 2 0 3", "30 1 0 0", "10 0 0 0", "20 0 1 0",
    "$EndNodes", "$Elements", "1 1", "22 2 2 1", "90 10 30 20", "$EndElements",
    "$Periodic", "1", "1 12 33", "3", "30 10", "20 10", "30 10", "$EndPeriodic"
  ), path)
  expect_error(mio_read(path, format = "gmsh"), "Periodic")
  r <- mio_read_with_info(path, format = "gmsh")
  copied <- mio_gmsh_info(r$info)
  expect_equal(copied$bounding_entities, list())
  expect_equal(copied$periodic[[1]]$entity, c(1, 12, 33))
  expect_equal(copied$periodic[[1]]$affine, numeric())
  expect_equal(copied$periodic[[1]]$node_pairs, matrix(c(1, 2, 3, 2, 1, 2), 2, 3))
  for (format in c("gmsh", "gmsh22")) {
    mio_write_with_info(r$mesh, r$info, output, format = format)
    back <- mio_read_with_info(output, format = "gmsh")
    expect_equal(mio_gmsh_info(back$info)$periodic, copied$periodic)
    mio_release(back$mesh)
    mio_format_info_release(back$info)
  }
  expect_error(mio_mdpa_info(r$info))
  mio_release(r$mesh)
  mio_format_info_release(r$info)
  expect_equal(copied$periodic[[1]]$node_pairs[2, 2], 2)
  expect_error(mio_gmsh_info(r$info), "released")
})
