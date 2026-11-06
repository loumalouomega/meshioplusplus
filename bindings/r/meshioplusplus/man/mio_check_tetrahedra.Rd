\name{mio_check_tetrahedra}
\alias{mio_check_tetrahedra}
\title{Whether the tetrahedra of a mesh are usable}
\usage{
mio_check_tetrahedra(mesh, zero_volume_factor = NULL)
}
\arguments{
\item{mesh}{A mesh.}
\item{zero_volume_factor}{The volume tolerance as a fraction of the longest
edge cubed (default \code{1e-6}); finite and not negative.}
}
\value{
A list: \code{num_tetrahedra}, \code{num_inverted}, \code{num_degenerate},
\code{num_overlaps}, \code{min_volume}, \code{max_volume}, \code{valid} (no
inverted or degenerate tetrahedron and no overlap) and \code{warnings}, a
character vector with one entry per warning.
}
\description{
Examines the \code{tetra} and \code{tetra10} cells (a quadratic tetrahedron by
its four corners) and counts the tetrahedra that are inverted (a signed volume
below minus the tolerance) or degenerate (a volume within it), and the pairs
that overlap across a shared face (no other tetrahedron uses it, and their
opposite nodes lie on the same side of it). The tolerance is
\code{zero_volume_factor} times the cube of the longest edge plus the effect of
six significant digits of written coordinates. The mesh is not changed.
}
