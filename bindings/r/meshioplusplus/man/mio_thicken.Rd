\name{mio_thicken}
\alias{mio_thicken}
\title{Thicken a surface into a closed shell}
\usage{
mio_thicken(mesh, thickness, side = "inward", smooth_iterations = 0)
}
\arguments{
\item{mesh}{A surface mesh.}

\item{thickness}{The wall's thickness.}

\item{side}{\code{"inward"}, \code{"outward"} or \code{"both"}.}

\item{smooth_iterations}{Smoothing passes over the directions.}
}
\value{
As \code{mio_offset_surface}.
}
\description{
The surface offset to the two sides of a wall \code{thickness} thick, the lower layer reversed and the boundary edges joined by quads; cell data \code{thicken:side} (roadmap 7.8.4). See \code{doc/offset_surface.md}.
}
