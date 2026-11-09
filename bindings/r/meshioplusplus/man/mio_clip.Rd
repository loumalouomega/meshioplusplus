\name{mio_clip}
\alias{mio_clip}
\title{Cut a surface by a plane, keep a side and cap it}
\usage{
mio_clip(mesh, origin, normal, keep_below = TRUE, cap = TRUE)
}
\arguments{
\item{mesh}{A surface mesh.}

\item{origin, normal}{The plane.}

\item{keep_below}{Keep the side below (else above).}

\item{cap}{Triangulate the cut's closed loops.}
}
\value{
A list: \code{mesh}, \code{num_kept}, \code{num_split}, \code{num_cap_triangles}, \code{num_loops}, \code{num_open_chains} and \code{warnings}.
}
\description{
The plane through \code{origin} normal to \code{normal}; the side below it (or above) is kept, straddling triangles split, and the cut's closed loops capped (roadmap 7.8.3). See \code{doc/clip.md}.
}
