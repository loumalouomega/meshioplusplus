\name{mio_face_support}
\alias{mio_face_support}
\title{The scan's support of each face of a fitted solid}
\usage{
mio_face_support(
  solid,
  scan,
  tolerance = 0,
  measured_share = 0.9,
  inferred_share = 0.1,
  hole_data = "repair:hole"
)
}
\arguments{
\item{solid}{A \code{fit_cad} re-mesh.}

\item{scan}{The scan.}

\item{tolerance}{Within this of the scan a sample is measured (\code{0}: \code{1e-3} of the scan's diagonal).}

\item{measured_share, inferred_share}{The class thresholds.}

\item{hole_data}{The scan cell data marking fills.}
}
\value{
A list: \code{mesh}, \code{num_measured}, \code{num_partial}, \code{num_inferred}, \code{num_faces}, \code{share_within} and \code{tolerance}.
}
\description{
Class each face of a \code{mio_fit_cad} re-mesh as measured, partial or inferred (roadmap 7.7.1). See \code{doc/fit_cad_loop.md}.
}
