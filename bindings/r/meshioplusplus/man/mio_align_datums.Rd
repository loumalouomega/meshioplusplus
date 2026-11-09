\name{mio_align_datums}
\alias{mio_align_datums}
\title{The 3-2-1 datum alignment of a measured mesh}
\usage{
mio_align_datums(mesh, primary, secondary, origin = NULL)
}
\arguments{
\item{mesh}{A mesh carrying an element table (\code{mio_measure()}).}

\item{primary, secondary}{Element names.}

\item{origin}{An element name, or \code{NULL} for the primary's point.}
}
\value{
A list: \code{mesh} and the 4 x 4 \code{matrix}.
}
\description{
\code{primary} sets Z, \code{secondary} X, \code{origin} (a point-like element or a third plane) the origin (roadmap 7.3.4).
}
