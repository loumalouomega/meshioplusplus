\name{mio_error_size}
\alias{mio_error_size}
\title{A size field from a per-cell error estimate}
\usage{
mio_error_size(mesh, target_error, error_array = "error:zz", order = 1,
  min_size = 0, max_size = 0, output = "size", overwrite = FALSE)
}
\arguments{
\item{mesh}{A mesh carrying the indicator as cell data.}
\item{target_error}{The target global error, in the indicator's norm (positive).}
\item{error_array}{The one-component cell data indicator.}
\item{order}{The indicator's rate in the cell size (1 for linear elements).}
\item{min_size, max_size}{The sizes' bounds; 0 is none.}
\item{output}{The point data array written.}
\item{overwrite}{Replace an existing array of that name.}
}
\value{
A list: \code{mesh}, \code{num_cells} (evaluated), \code{num_skipped},
\code{num_clamped}, \code{min_size} and \code{max_size} (the least and greatest
size written).
}
\description{
A per-point size field from the per-cell error indicator \code{error_array}
(what \code{\link{mio_estimate_error}} writes) by the adaptive procedure of
Zienkiewicz and Zhu (1987; roadmap 6.10.2): with \code{N} cells evaluated, a
cell of diameter \code{h} whose indicator is \code{eta} asks
\code{h * (target_error / sqrt(N) / eta)^(1 / order)}, clamped to
\code{[min_size, max_size]} (0: no bound); a point takes the least of its
cells' (\code{NaN} where none is evaluated). The sizes are written as point
data \code{output}, what \code{\link{mio_mesh_cad_surface}}'s \code{background}
reads.
}
