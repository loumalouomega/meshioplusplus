# Notebook-only display helpers for the meshio++ R examples.
#
# This binding rides on the flat C API (see doc/r.md), and CLAUDE.md documents
# that data-driven SVG colouring is a flat-ABI gap: registry.cpp's (path, mesh)
# writer lambdas can't carry per-call parameters, so mio_write() always emits
# the fixed default styling -- no color_by, no colorbar, no camera control
# (unlike the C++ notebooks, which call write_svg directly and bypass the
# registry). So where the C++ notebooks colour a render by a quality metric or
# a field, these notebooks show the same information as a small hand-rolled
# SVG chart instead -- honest about the gap rather than working around it.
#
# Nothing here is part of the package API; it exists only to give these
# notebooks something to display via IRdisplay.

library(meshioplusplus)

.mio_nb_counter <- local({
  n <- 0
  function() {
    n <<- n + 1
    n
  }
})

#' A fresh scratch file path with the given suffix, for one-shot round trips.
mio_nb_temp_path <- function(suffix) {
  file.path(tempdir(), sprintf("meshioplusplus_r_nb_%d%s", .mio_nb_counter(), suffix))
}

#' Write `mesh` through the C API's SVG writer (fixed default camera and
#' styling -- see the module note above) and display it inline.
render <- function(mesh) {
  path <- mio_nb_temp_path(".svg")
  mio_write(mesh, path, format = "svg")
  IRdisplay::display_svg(file = path)
}

#' An n x n x n hexahedron block on the lattice `[0, n*spacing]^3` -- the
#' small, easy-to-read fixture several operation demos build on
#' (convert_cells, refine, smooth, interpolate), mirroring the C++ and Julia
#' notebooks' `hex_block`.
hex_block <- function(n, spacing = 1.0) {
  pid <- function(i, j, k) i * (n + 1)^2 + j * (n + 1) + k + 1 # 1-based flat index
  npts <- (n + 1)^3
  pts <- matrix(0, nrow = 3, ncol = npts)
  for (i in 0:n) {
    for (j in 0:n) {
      for (k in 0:n) {
        pts[, pid(i, j, k)] <- c(i, j, k) * spacing
      }
    }
  }

  ncells <- n^3
  conn <- matrix(0, nrow = 8, ncol = ncells)
  c <- 0
  for (i in 0:(n - 1)) {
    for (j in 0:(n - 1)) {
      for (k in 0:(n - 1)) {
        c <- c + 1
        conn[, c] <- c(
          pid(i, j, k), pid(i + 1, j, k), pid(i + 1, j + 1, k), pid(i, j + 1, k),
          pid(i, j, k + 1), pid(i + 1, j, k + 1), pid(i + 1, j + 1, k + 1), pid(i, j + 1, k + 1)
        )
      }
    }
  }

  m <- mio_mesh()
  mio_set_points(m, pts)
  mio_add_cell_block(m, "hexahedron", conn)
  m
}

#' Horizontal bar chart -- e.g. file sizes per format, cell counts per type.
#' Uses base R graphics; IRkernel captures the device output as an inline PNG,
#' so this needs no extra package.
mio_bar_chart <- function(labels, values, title, color = "#3b82f6") {
  op <- par(mar = c(4, max(nchar(labels)) * 0.6 + 4, 3, 4))
  on.exit(par(op))
  barplot(rev(values),
    names.arg = rev(labels), horiz = TRUE, las = 1, col = color,
    main = title, cex.names = 0.9
  )
}

#' Histogram of a numeric vector -- base R graphics, no extra package.
mio_histogram <- function(values, title, xlabel, color = "#3b82f6") {
  hist(values,
    main = title, xlab = xlabel, col = color, border = "white",
    breaks = 12
  )
}
