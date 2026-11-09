\name{mio_measure_element_primitive}
\alias{mio_measure_element_primitive}
\title{A measured element as a fitted primitive}
\usage{
mio_measure_element_primitive(mesh, name)
}
\arguments{
\item{mesh}{A mesh carrying an element table.}

\item{name}{The element's name.}
}
\value{
A list: \code{type} and \code{params}.
}
\description{
The element \code{name} as a primitive type (1 plane, 2 sphere, 3 cylinder, 4 cone) and ten parameters, a deviation reference (roadmap 7.4.2).
}
