\name{mio_find_symmetry_plane}
\alias{mio_find_symmetry_plane}
\title{Symmetry plane of a mesh}
\usage{
mio_find_symmetry_plane(mesh, seed = NULL, tolerance = 0, min_share = 0.9,
  max_samples = 4000)
}
\arguments{
\item{mesh}{A \code{mio_mesh}.}

\item{seed}{\code{NULL}, or \code{list(origin, normal)}.}

\item{tolerance}{A sample within this of the surface counts as symmetric; 0 takes 0.005 of the diagonal.}

\item{min_share}{Refuse a plane below this share.}

\item{max_samples}{The samples of the fine fit.}
}
\value{
A list of \code{origin}, \code{normal}, \code{rms}, \code{share}, \code{tolerance} and \code{chosen_candidate}.
}
\description{
Each principal plane through the centroid (and an optional seed) is mirrored, registered back by ICP and read off the reflection-motion product (roadmap 7.3.6). See \code{doc/registration.md}.
}
