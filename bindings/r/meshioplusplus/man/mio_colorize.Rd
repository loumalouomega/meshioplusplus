\name{mio_colorize}
\alias{mio_colorize}
\title{Bake a field into red, green and blue colours}
\usage{
mio_colorize(mesh, field, cmap = "viridis", vmin = NULL, vmax = NULL,
  symmetric = FALSE, steps = 0, nan_color = "#808080", component = NULL)
}
\arguments{
\item{mesh}{A mesh carrying \code{field} as point or cell data.}
\item{field}{The array to colour by.}
\item{cmap}{A colormap name.}
\item{vmin, vmax}{The range; \code{NULL} takes the field's min and max.}
\item{symmetric}{With the automatic range, \code{+-max|v|}.}
\item{steps}{Bands at fixed steps (0: continuous).}
\item{nan_color}{\code{"#rrggbb"} for a non-finite value.}
\item{component}{Colour by this 1-based component (\code{NULL}: the magnitude).}
}
\value{
A list: \code{mesh}, \code{vmin}, \code{vmax} and \code{location} (\code{"point"} or \code{"cell"}).
}
\description{
Maps \code{field} through a colormap into the UInt8 data \code{red}, \code{green} and
\code{blue} on the field's location (roadmap 7.4.4); \code{steps} of 2 or more bands
the map at fixed steps.
}
