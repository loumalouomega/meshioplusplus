\name{mio_split_interface}
\alias{mio_split_interface}
\title{Split a Mesh Along a Side Region}
\usage{
mio_split_interface(mesh, side, add_cohesive = FALSE)
}
\arguments{
\item{mesh}{Input mesh.}
\item{side}{Name of a Side region on the input mesh.}
\item{add_cohesive}{Insert cohesive cells on selected two-owner facets.}
}
\value{
A list containing the split \code{mesh}, \code{num_duplicated_points}, and
\code{num_cohesive_cells}.
}
\description{
Duplicates incident-cell fans along a Side region while preserving data and
regions. Polyhedron inputs are rejected. See \code{doc/region_adjacency.md}.
}
