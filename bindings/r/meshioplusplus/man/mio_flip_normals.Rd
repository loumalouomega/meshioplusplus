\name{mio_flip_normals}
\alias{mio_flip_normals}
\title{Reverse the winding of cells}
\usage{
mio_flip_normals(mesh, region = NULL)
}
\arguments{
\item{mesh}{A mesh.}

\item{region}{A cell region's name, or \code{NULL} for every cell.}
}
\value{
A list: \code{mesh} and \code{num_flipped}.
}
\description{
Every cell, or the cells of a cell region, with no point moved; Side region entries move with their cells (roadmap 7.8.2). See \code{doc/flip_normals.md}.
}
