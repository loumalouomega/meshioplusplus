\name{mio_deviation}
\alias{mio_deviation}
\title{Signed deviation of a scan from a reference}
\usage{
mio_deviation(scan, reference = NULL, type = NULL, params = NULL, tolerance = 0,
  bins = 20, range = NULL, max_normal_angle = 90, reference_region = "",
  register = FALSE)
}
\arguments{
\item{scan}{A point cloud or surface mesh.}
\item{reference}{A surface mesh, or \code{NULL} for a primitive.}
\item{type}{\code{"plane"}, \code{"sphere"}, \code{"cylinder"}, \code{"cone"} or \code{"torus"}.}
\item{params}{The primitive's ten parameters.}
\item{tolerance}{Report the share within \code{+-tolerance} (0: none).}
\item{bins}{Histogram bins.}
\item{range}{\code{c(lo, hi)} of the histogram, or \code{NULL} for \code{+-max|d|}.}
\item{max_normal_angle}{The far-side rule's angle in degrees (180: off).}
\item{reference_region}{Named cell region of the reference.}

\item{register}{Register the scan to the reference first (roadmap 7.4.1): \code{mesh} is then the registered scan, and the list adds \code{registered}, \code{matrix} and \code{registration_rms}.}
}
\value{
A list: \code{mesh} (with the \code{deviation:distance} point data),
\code{num_points}, \code{num_measured}, \code{min}, \code{max}, \code{mean},
\code{rms}, \code{sigma}, \code{percentiles} (p1, p5, p50, p95, p99),
\code{histogram_edges}, \code{histogram_counts}, \code{num_outside_range},
\code{share_within}, \code{suggested_range}, \code{num_without_normal} and
\code{warnings}.
}
\description{
The signed distance of every point of \code{scan} to the surface
\code{reference} (positive outside it, excess material), or, with
\code{reference = NULL}, to a fitted primitive \code{type} with its ten
\code{fit:patch_params} entries (roadmap 7.4), and the statistics an
inspection reads from it.
}
