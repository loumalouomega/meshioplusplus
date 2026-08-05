# The fixture mirrors tests/fortran/test_fortran_api.f90 and the Julia
# binding's test suite deliberately: 5 points x 3 dims, 2 tetra x 4 nodes,
# 3-component vector data, and coordinates that are asymmetric in every axis.
# Nothing about it is square, so a transposed array mapping or a missed 1-based
# shift cannot cancel out and pass anyway.

POINTS <- matrix(
  c(
    0.0, 0.0, 0.0,
    1.1, 0.2, 0.3,
    0.4, 1.2, 0.5,
    0.6, 0.7, 1.3,
    1.4, 1.5, 1.6
  ),
  nrow = 3
)

# (nodes_per_cell x num_cells), 1-based.
CONN <- matrix(c(1, 2, 3, 4, 2, 3, 4, 5), nrow = 4)

VEC <- matrix(
  c(
    10, 11, 12,
    20, 21, 22,
    30, 31, 32,
    40, 41, 42,
    50, 51, 52
  ),
  nrow = 3
)

TEMPERATURE <- c(1, 2, 3, 4, 5)

fixture <- function() {
  m <- mio_mesh()
  mio_set_points(m, POINTS)
  mio_add_cell_block(m, "tetra", CONN)
  mio_add_point_data(m, "displacement", VEC)
  mio_add_point_data(m, "temperature", TEMPERATURE)
  mio_append_cell_data(m, "material", c(7, 9))
  m
}

# The unit cube as a closed, outward-wound triangle surface. The distance tests
# need a watertight surface, and the tetra fixture above is a volume.
cube_surface <- function() {
  m <- mio_mesh()
  mio_set_points(m, matrix(c(
    0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,
    0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1
  ), nrow = 3))
  mio_add_cell_block(m, "triangle", matrix(c(
    1, 3, 2, 1, 4, 3, 5, 6, 7, 5, 7, 8,
    1, 2, 6, 1, 6, 5, 2, 3, 7, 2, 7, 6,
    3, 4, 8, 3, 8, 7, 4, 1, 5, 4, 5, 8
  ), nrow = 3))
  m
}
