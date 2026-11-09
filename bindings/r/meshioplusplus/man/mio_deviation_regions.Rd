\name{mio_deviation_regions}
\alias{mio_deviation_regions}
\title{Excess and missing-material regions beyond a deviation tolerance}
\usage{
mio_deviation_regions(mesh, tolerance, field = "deviation:distance",
  prefix = "deviation:", min_area = 0)
}
\arguments{
\item{mesh}{A surface carrying \code{field} (see \code{mio_deviation()}).}
\item{tolerance}{The threshold, positive.}
\item{field}{The signed deviation point data.}
\item{prefix}{Region names are \code{<prefix>excess_<k>} and \code{<prefix>missing_<k>}.}
\item{min_area}{Drop regions smaller than this area.}
}
\value{
A list: \code{mesh}, per region (largest \code{|volume|} first) \code{name},
\code{sign}, \code{num_cells}, \code{area}, \code{extreme}, \code{mean},
\code{volume}, \code{centroid} (n x 3) and \code{bounds} (n x 6), then
\code{num_excess}, \code{num_missing}, \code{num_dropped} and \code{warnings}.
}
\description{
The connected regions (cells joined across shared edges) of a measured
surface whose deviation, the mean of each cell's corners' \code{field} values,
lies beyond \code{+tolerance} (excess) or \code{-tolerance} (missing material)
(roadmap 7.4.3), written as named cell regions with the \code{deviation:region}
cell data.
}
