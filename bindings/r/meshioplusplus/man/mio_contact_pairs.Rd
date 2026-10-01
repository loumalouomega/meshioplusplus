\name{mio_contact_pairs}
\alias{mio_contact_pairs}
\title{Project Slave Points to a Master Boundary}
\usage{
mio_contact_pairs(slave_mesh, slave_points, master_cells,
  master_mesh = NULL, tolerance = 0, require_complete = FALSE)
}
\arguments{
\item{slave_mesh}{Mesh containing the slave Point region.}
\item{slave_points}{Point-region name.}
\item{master_cells}{Master Cell-region name.}
\item{master_mesh}{Optional second mesh; defaults to \code{slave_mesh}.}
\item{tolerance}{Maximum accepted projection distance; zero derives a scale.}
\item{require_complete}{Fail if any slave point is unmatched.}
}
\value{
A list of copied arrays: 1-based slave point and master cell ids, 0-based facet
and subfacet ordinals, local coordinates, closest points, signed gaps, normals,
and unmatched point ids.
}
\description{
Projects each Point-region entry to the closest queryable master facet. See
\code{doc/region_adjacency.md}.
}
