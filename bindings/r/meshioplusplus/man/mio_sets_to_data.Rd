\name{mio_sets_to_data}
\alias{mio_sets_to_data}
\alias{mio_data_to_sets}
\title{Convert sets and integer data}
\usage{
mio_sets_to_data(mesh, location, data_name = NULL, join_char = "-", order = NULL)

mio_data_to_sets(mesh, location, key)
}
\arguments{
\item{mesh}{A \code{mio_mesh} object.}
\item{location}{One of \code{"point"}, \code{"cell"}.}
\item{data_name}{Output field name, or \code{NULL} to join set names.}
\item{join_char}{Separator used when joining set names.}
\item{order}{Every set name exactly once, or \code{NULL} for native region-name order.}
\item{key}{Source scalar integer field to convert and remove.}
}
\value{A new \code{mio_mesh}.}
\description{
Convert point/cell region-backed sets into scalar integer labels, or one
scalar integer field into sets. Both return a new mesh. Labels are zero-based
(not indices): later overlapping sets win; uncovered entities have label -1.
Side regions, geometry, unrelated data and property sets are preserved.
}
