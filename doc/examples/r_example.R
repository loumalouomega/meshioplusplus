# Minimal meshio++ R consumer: build a tet mesh, write it, read it back.
#
# Run with the C library on the loader path, e.g.
#   LD_LIBRARY_PATH=/opt/meshioplusplus/lib Rscript doc/examples/r_example.R

library(meshioplusplus)

cat(sprintf("meshio++ %s (backend: %s)\n", mio_version(), mio_mesh_backend()))

# (dim x num_points) and (nodes_per_cell x num_cells): column-major R is the
# same layout as the C API's row-major (num_points, dim) / (num_cells, npc), so
# nothing is transposed. Connectivity is 1-based here.
points <- matrix(
  c(
    0.0, 0.0, 0.0,
    1.1, 0.2, 0.3,
    0.4, 1.2, 0.5,
    0.6, 0.7, 1.3,
    1.4, 1.5, 1.6
  ),
  nrow = 3
)
conn <- matrix(c(1, 2, 3, 4, 2, 3, 4, 5), nrow = 4)

path <- file.path(tempdir(), "mio_example_r.vtu")

m <- mio_mesh()
mio_set_points(m, points)
mio_add_cell_block(m, "tetra", conn)
mio_write(m, path)
mio_release(m)

r <- mio_read(path)
stopifnot(
  mio_num_points(r) == 5,
  isTRUE(all.equal(mio_points(r), points)),
  mio_cell_block_type(r, 1) == "tetra",
  identical(mio_connectivity(r, 1), conn),
  # The _raw reader is 0-based -- the ABI's own numbering. It still copies:
  # R has no zero-copy borrow.
  identical(mio_connectivity_raw(r, 1), conn - 1)
)
mio_release(r)
unlink(path)

cat("r_example.R: OK\n")
