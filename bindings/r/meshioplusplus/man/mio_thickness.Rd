\name{mio_thickness}
\alias{mio_thickness}
\title{Wall thickness or clearance per facet}
\usage{
mio_thickness(mesh, other = NULL, direction = "auto", max_distance = 0,
  threshold = 0, region = "", region_other = "")
}
\arguments{
\item{mesh}{A surface mesh.}
\item{other}{The surface to reach, or \code{NULL} for \code{mesh} itself.}
\item{direction}{\code{"auto"}, \code{"inward"} or \code{"outward"}.}
\item{max_distance}{The longest ray (0: the bounding-box diagonal).}
\item{threshold}{Count the hits below it (0: none).}
\item{region, region_other}{Named cell regions measured and reached.}
}
\value{
A list: \code{mesh} (with the \code{thickness:distance} and
\code{thickness:hit_cell} cell data), \code{direction}, \code{num_cells},
\code{num_hit}, \code{num_missed}, \code{num_below}, \code{min}, \code{max},
\code{mean}, \code{min_cell} (1-based, 0 when nothing was hit),
\code{min_point}, \code{max_distance} and \code{warnings}.
}
\description{
Casts a ray from each 2-D cell's centroid along its normal (roadmap 3.1.4):
\code{"inward"} stops at the far side of the wall, \code{"outward"} at the next
wall facing the ray, and \code{"auto"} is inward for one surface, outward
against \code{other}. The surface must be consistently oriented with outward
normals.
}
