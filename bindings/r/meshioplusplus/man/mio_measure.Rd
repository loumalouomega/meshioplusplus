\name{mio_measure}
\alias{mio_measure}
\title{Elements, dimensions and GD&T tolerances}
\usage{
mio_measure(mesh, spec)
}
\arguments{
\item{mesh}{The scanned surface.}

\item{spec}{The spec as JSON text, or the path of a JSON file.}
}
\value{
A list: \code{mesh} (with its element table), \code{num_elements}, \code{num_dimensions}, \code{num_tolerances}, \code{num_fail}, \code{report_json} and \code{warnings}.
}
\description{
Measure a mesh by a JSON spec (roadmap 7.2): elements fitted from regions, seeds or patches, constructed or literal; dimensions and tolerances. See \code{doc/measure.md}.
}
