\name{mio_gmsh_info}
\alias{mio_gmsh_info}
\title{Copy a Gmsh format side channel}
\usage{mio_gmsh_info(info)}
\arguments{\item{info}{A format-info handle from \code{mio_read_with_info()}.}}
\description{Copy bounding-entity tags and periodic links into plain R values.}
\value{A list with \code{bounding_entities} (signed tags per cell block) and
\code{periodic}. Each periodic link has \code{entity} (dimension, slave tag,
master tag), \code{affine} (0 or 16 coefficients), and \code{node_pairs} (2,N)
with 1-based point rows. Entity tags remain file ids. Ordering and duplicates
survive I/O; mesh operations do not remap these pairs. The copy remains valid
after releasing the handle.}
