\name{mio_find_interface}
\alias{mio_find_interface}
\title{Find a Conforming or Proximity-Matched Interface}
\usage{
mio_find_interface(mesh, region_a, region_b, mesh_b = NULL,
  mode = "conforming", master = "a", gap_tolerance = 0,
  angle_tolerance = 30, overlap_tolerance = 0)
}
\arguments{
\item{mesh}{A mesh containing region A.}
\item{region_a, region_b}{Cell-region names.}
\item{mesh_b}{Optional second mesh containing region B.}
\item{mode}{\code{"conforming"} or \code{"proximity"}.}
\item{master}{\code{"a"} or \code{"b"}; selects the output facets and signed-gap normal.}
\item{gap_tolerance, angle_tolerance, overlap_tolerance}{Matching tolerances.}
}
\value{
A list containing the interface \code{mesh}, pair counts and measure, unmatched
boundary-facet counts, and \code{side_a}/\code{side_b} matrices. Cell ids in the
Side matrices are 1-based; local-facet ordinals remain 0-based.
}
\description{
Finds conforming or proximity-matched boundary facets between two Cell regions.
See \code{doc/region_adjacency.md}.
}
