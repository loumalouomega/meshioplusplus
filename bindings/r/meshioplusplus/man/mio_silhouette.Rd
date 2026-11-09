\name{mio_silhouette}
\alias{mio_silhouette}
\title{The outline of a part along a direction}
\usage{
mio_silhouette(mesh, normal = c(0, 0, 1), clearance = 0, round_angle = 180,
  arc_segments = 64, extrude = FALSE, region = NULL, merge_tolerance = 1e-9)
}
\arguments{
\item{mesh}{A surface or volume mesh.}

\item{normal}{The viewing direction.}

\item{clearance}{Grows the outline by this distance (negative shrinks it).}

\item{round_angle}{Corners turning by less than this many degrees are rounded.}

\item{arc_segments}{Polygon sides of a full turn of arc.}

\item{extrude}{Return the prism's closed triangle surface.}

\item{region}{Restrict to this cell region of a surface.}

\item{merge_tolerance}{Projected points this close, relative to the box diagonal, are one.}
}
\value{
A list: \code{mesh}, \code{num_loops}, \code{num_outer}, \code{num_holes}, \code{num_triangles}, \code{area}, \code{perimeter}, \code{min_height}, \code{max_height} and \code{warnings}.
}
\description{
The surface (or a volume mesh's skin) projected onto the plane normal to \code{normal} and unioned exactly, as closed chains with its holes, ready for \code{mio_fit_curves()}; \code{clearance} offsets it, and \code{extrude = TRUE} returns the closed surface of the prism over it instead (roadmap 7.6.3). See \code{doc/silhouette.md}.
}
