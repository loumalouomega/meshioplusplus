\name{mio_mirror}
\alias{mio_mirror}
\title{Reflect a mesh across a plane}
\usage{
mio_mirror(mesh, origin = c(0, 0, 0), normal = c(1, 0, 0), keep_original = FALSE,
  weld = TRUE, tolerance = 1e-9, reflect_vector_data = FALSE)
}
\arguments{
\item{mesh}{A mesh.}
\item{origin}{A point of the plane, \code{c(x, y, z)}.}
\item{normal}{The plane's normal.}
\item{keep_original}{Return the input followed by its image.}
\item{weld}{With \code{keep_original}, share the points on the plane.}
\item{tolerance}{On the plane: a fraction of the bounding-box diagonal.}
\item{reflect_vector_data}{Reflect float vector and tensor data.}
}
\value{
A list: \code{mesh}, \code{num_flipped_cells}, \code{num_points_on_plane},
\code{num_dropped_cells}, \code{num_dropped_sides} and \code{warnings}.
}
\description{
Reflects \code{mesh} across the plane through \code{origin} with
\code{normal}, every cell right side out (roadmap 3.2.5): its nodes permuted
and \code{Side} region entries moved to their facets' new numbers. With
\code{keep_original} the input followed by its image, the points on the plane
shared when \code{weld}.
}
