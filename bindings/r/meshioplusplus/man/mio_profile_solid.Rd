\name{mio_profile_solid}
\alias{mio_profile_solid}
\title{Extrude or revolve a sketch into a solid}
\usage{
mio_profile_solid(
  sketch,
  features,
  step_path = NULL,
  scan = NULL,
  provider = "",
  deflection = 0
)
}
\arguments{
\item{sketch}{A sketch.}

\item{features}{The features, in text form.}

\item{step_path}{Where the solid goes as STEP (\code{NULL}: nowhere).}

\item{scan}{A scan for an extrusion to its extent, or \code{NULL}.}

\item{provider}{The provider (\code{""}: the first).}

\item{deflection}{The tessellation's deflection (\code{0}: \code{1e-3} of the diagonal).}
}
\value{
A list: \code{mesh}, \code{volume}, \code{area}, \code{num_faces}, \code{valid} and \code{warnings}.
}
\description{
A \code{mio_fit_curves} sketch built into a solid by features applied in order, through a CAD provider (OpenCASCADE; roadmap 7.7.4): \code{extrude:D[:symmetric]}, \code{revolve:OX,OY,OZ:DX,DY,DZ[:ANGLE]}, \code{fillet:R}, \code{chamfer:D}, \code{shell:W}, separated by \code{;}. See \code{doc/profile_solid.md}.
}
