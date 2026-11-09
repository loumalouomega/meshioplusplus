\name{mio_sections}
\alias{mio_sections}
\title{Chained sections at parallel planes}
\usage{
mio_sections(mesh, direction, offsets, min_length = 0, record_parent_ids = FALSE)
}
\arguments{
\item{mesh}{A surface mesh.}

\item{direction}{The planes' normal.}

\item{offsets}{The planes' signed distances from the origin.}

\item{min_length}{Drop chains shorter than this.}

\item{record_parent_ids}{Attach \code{slice:parent_cell}.}
}
\value{
A \code{mio_mesh}: the merged line mesh with the \code{section:offset}, \code{section:chain} and \code{section:closed} cell data.
}
\description{
The sections of a surface at the planes \code{n . x = offsets[k]} (\code{n} the unit \code{direction}), each walked into polylines (roadmap 7.6.1). See \code{doc/slice.md}.
}
