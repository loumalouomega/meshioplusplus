\name{mio_find_intersections}
\alias{mio_find_intersections}
\title{Find where a surface intersects itself, or two surfaces meet}
\usage{
mio_find_intersections(mesh, other = NULL, region = "", region_other = "",
  segments = TRUE, max_pairs = 0)
}
\arguments{
\item{mesh}{A surface mesh (run \code{mio_extract_surface} on a volume mesh).}
\item{other}{A second surface, or \code{NULL} to check \code{mesh} against itself.}
\item{region, region_other}{Named cell regions to restrict each surface to.}
\item{segments}{Build the line mesh of the intersection segments.}
\item{max_pairs}{Keep at most this many pairs (0: all).}
}
\value{
A list: \code{mesh} (with the \code{intersections:count} cell data),
\code{other} (likewise, or \code{NULL}), \code{segments} (a line mesh, or
\code{NULL}), \code{pairs} (an \code{n x 2} matrix of 1-based global cell ids),
\code{num_pairs}, \code{num_triangle_pairs}, \code{num_coplanar},
\code{num_segments}, \code{num_cells_a}, \code{num_cells_b}, \code{truncated},
\code{clean} and \code{warnings}.
}
\description{
Reduces the surfaces to triangles and decides every candidate pair with an
exact triangle-triangle test (roadmap 3.1.2), so every hit is exact. On one
surface, triangles sharing corners count only when they meet beyond them;
against \code{other}, any contact counts, touching included.
}
