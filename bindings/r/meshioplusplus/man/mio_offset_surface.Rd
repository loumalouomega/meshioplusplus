\name{mio_offset_surface}
\alias{mio_offset_surface}
\title{Offset a surface}
\usage{
mio_offset_surface(mesh, distance, smooth_iterations = 0)
}
\arguments{
\item{mesh}{A surface mesh.}

\item{distance}{Along the faces' normals; negative against them.}

\item{smooth_iterations}{Smoothing passes over the directions.}
}
\value{
A list: \code{mesh}, \code{num_not_visible}, \code{num_clamped}, \code{num_folded}, \code{num_rim_quads}, \code{max_move} and \code{warnings}.
}
\description{
Move each point of a triangle and quad surface along its offset direction (the one \code{mio_boundary_layer} grows along) so every face's plane moves by \code{distance}; a cube offset by \code{d} is the cube grown by \code{d} (roadmap 7.8.4). See \code{doc/offset_surface.md}.
}
