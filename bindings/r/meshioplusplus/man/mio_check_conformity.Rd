\name{mio_check_conformity}
\alias{mio_check_conformity}
\title{Whether a mesh is conforming}
\usage{
mio_check_conformity(mesh, tolerance = 1e-9, attach_flags = TRUE)
}
\arguments{
\item{mesh}{A mesh.}
\item{tolerance}{Distances at or below this fraction of the examined facets'
bounding-box diagonal count as zero.}
\item{attach_flags}{Attach the \code{conformity:hanging} point data.}
}
\value{
A list: \code{mesh}, \code{dimension}, \code{num_hanging_nodes},
\code{num_t_junctions}, \code{num_partial_overlaps},
\code{num_unwelded_facets}, \code{num_boundary_facets},
\code{num_non_manifold_facets}, \code{conforming}, \code{interface_mode}
(\code{"conforming"}, \code{"weld"} or \code{"proximity"}: the
\code{\link{mio_find_interface}} mode that suits the mesh) and
\code{warnings}.
}
\description{
Checks the mesh's volume cells (else its 2-D cells) for hanging nodes (a point
on a volume cell's face or edge that the cell does not use), T-junctions (the
same on a 2-D cell's edge), and partly overlapping or unwelded pairs of
unmatched facets (roadmap 3.2.1). With \code{attach_flags} the mesh carries the
\code{conformity:hanging} point data (1 hanging node, 2 T-junction).
}
