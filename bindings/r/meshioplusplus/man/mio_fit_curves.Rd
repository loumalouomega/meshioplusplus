\name{mio_fit_curves}
\alias{mio_fit_curves}
\title{Lines, arcs and circles from chained sections}
\usage{
mio_fit_curves(sections, tolerance = 0, snap_angle = 0, normal = c(0, 0, 0),
  min_arc_points = 4, min_arc_sweep = 5)
}
\arguments{
\item{sections}{A chained line mesh.}

\item{tolerance}{The largest residual; 0 takes 1e-3 of the diagonal.}

\item{snap_angle}{Snap lines within this many degrees of the axes or diagonals.}

\item{normal}{The sketch normal; zeros take it from the points.}

\item{min_arc_points, min_arc_sweep}{An arc spans at least this many points and degrees.}
}
\value{
A list: \code{mesh}, \code{pieces} (a data frame), \code{num_lines}, \code{num_arcs}, \code{num_circles}, \code{max_residual}, \code{tolerance} and \code{warnings}.
}
\description{
Each chain of a chained line mesh split into the runs a line or a circular arc explains within the tolerance, in its plane's 2-D frame (roadmap 7.6.2). See \code{doc/fit_curves.md}.
}
