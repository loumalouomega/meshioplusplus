\name{mio_region_adjacency}
\alias{mio_region_adjacency}
\title{Conforming Shared Facets Between Cell Regions}
\usage{
mio_region_adjacency(mesh, regions = NULL)
}
\arguments{
\item{mesh}{A mesh.}
\item{regions}{Optional character vector of Cell-region names.}
}
\value{
A mesh of conforming shared facets, with interface parent-cell/facet ids and
per-facet length or area in cell data.
}
\description{
When \code{regions} is \code{NULL} or empty, all Cell regions are selected. If
fewer than two exist, multiple cell blocks are used as groups. See
\code{doc/region_adjacency.md}.
}
