\name{mio_hollow}
\alias{mio_hollow}
\title{Hollow a closed surface}
\usage{
mio_hollow(mesh, wall, smooth_iterations = 0)
}
\arguments{
\item{mesh}{A closed surface mesh.}

\item{wall}{The wall's thickness.}

\item{smooth_iterations}{Smoothing passes over the directions.}
}
\value{
As \code{mio_offset_surface}.
}
\description{
The surface with its offset by \code{-wall}, reversed, inside it; an open surface is refused (roadmap 7.8.4). See \code{doc/offset_surface.md}.
}
