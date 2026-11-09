\name{mio_remove_islands}
\alias{mio_remove_islands}
\title{Remove the small connected components of a mesh}
\usage{
mio_remove_islands(mesh, connectivity = "edge", min_diagonal = 0.01, min_area = 0,
  min_cells = 0, relative = TRUE, dry_run = FALSE)
}
\arguments{
\item{mesh}{A mesh.}

\item{connectivity}{\code{"edge"} or \code{"point"}.}

\item{min_diagonal, min_area, min_cells}{The thresholds (0 off).}

\item{relative}{Diagonal and area as fractions of the whole mesh's.}

\item{dry_run}{Remove nothing; label \code{islands:component} and \code{islands:removed}.}
}
\value{
A list: \code{mesh}, \code{num_components}, \code{num_removed}, \code{num_cells_removed} and \code{num_points_removed}.
}
\description{
Drop a scan's crumbs (roadmap 7.8.1): components (by shared edges or points) whose bounding-box diagonal, area or cell count is below a threshold. See \code{doc/remove_islands.md}.
}
