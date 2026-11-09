\name{mio_exaggerate}
\alias{mio_exaggerate}
\title{Exaggerate a deviation map for a figure}
\usage{
mio_exaggerate(mesh, field = "deviation:distance", scale = 10, clamp = 0,
  smooth_fraction = 0.01)
}
\arguments{
\item{mesh}{A surface carrying \code{field}.}
\item{field}{The signed deviation point data.}
\item{scale}{Displacement per unit deviation.}
\item{clamp}{Clamp the deviation first (0: automatic).}
\item{smooth_fraction}{Smoothing radius as a fraction of the diagonal.}
}
\value{
A list: \code{mesh}, \code{clamp}, \code{radius}, \code{num_moved} and \code{warnings}.
}
\description{
A copy of the surface \code{mesh} with each point displaced along its normal by
\code{scale} times its deviation, clamped to \code{+-clamp} (0: the suggested colour
range) and smoothed over \code{smooth_fraction} of the bounding-box diagonal
(roadmap 7.4.4).
}
