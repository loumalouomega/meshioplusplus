\name{mio_constrain_primitives}
\alias{mio_constrain_primitives}
\title{Constrain fitted primitives under control}
\usage{
mio_constrain_primitives(
  mesh,
  face_tolerance = 0,
  angle_tolerance = -1,
  distance_tolerance = -1,
  min_fraction = -1,
  force = NULL,
  forbid = NULL
)
}
\arguments{
\item{mesh}{A fitted skin.}

\item{face_tolerance}{How far a face may lie from its patch (\code{0}: three times the largest patch tolerance).}

\item{angle_tolerance, distance_tolerance, min_fraction}{In place of the pass's own (negative keeps them).}

\item{force, forbid}{Relations to force or forbid: lists of \code{list(patch_a, patch_b, relation)}.}
}
\value{
A list: \code{mesh}, \code{num_constrained}, \code{num_groups}, \code{face_tolerance}, \code{relations} (a data frame of \code{patch_a}, \code{patch_b}, \code{relation}, \code{angle_deg}, \code{distance}, \code{forced}) and \code{unmet}.
}
\description{
Fit related primitives of a \code{mio_fit_primitives} skin together: coaxial groups, parallel and perpendicular planes, equal radii (roadmap 7.7.3). See \code{doc/constrain_primitives.md}.
}
