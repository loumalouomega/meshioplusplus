\name{mio_extrude}
\alias{mio_extrude}
\alias{mio_revolve}
\title{Sweep a surface (or a curve) along a vector or around an axis}
\usage{
mio_extrude(mesh, direction = c(0, 0, 1), num_layers = 1, ratio = 1,
  offsets = numeric(0), record_layers = FALSE, max_cells = 20000000)

mio_revolve(mesh, axis = c(0, 1, 0), axis_point = c(0, 0, 0), angle = 360,
  num_layers = 16, tolerance = 1e-9, record_layers = FALSE,
  max_cells = 20000000)
}
\arguments{
\item{mesh}{A triangle and quad surface, or a curve of lines.}
\item{direction}{The whole displacement, \code{c(x, y, z)}.}
\item{num_layers}{Layers (steps of the turn).}
\item{ratio}{Each layer's thickness over the one below's.}
\item{offsets}{Explicit level fractions, or \code{numeric(0)}.}
\item{record_layers}{Attach \code{extrude:layer} (\code{revolve:layer}) cell data.}
\item{max_cells}{Refuse more profile cells times layers.}
\item{axis, axis_point}{The axis direction and a point of it.}
\item{angle}{Degrees, right-hand rule; 360 closes the turn.}
\item{tolerance}{On the axis: a fraction of the bounding-box diagonal.}
}
\value{
A list: \code{mesh}, \code{num_cells}, \code{num_hexahedra},
\code{num_wedges}, \code{num_pyramids}, \code{num_tetrahedra},
\code{num_quads}, \code{num_triangles}, \code{num_flipped},
\code{num_on_axis}, \code{min_scaled_jacobian}, \code{full} and
\code{warnings}.
}
\description{
\code{mio_extrude()} sweeps the triangles and quads of \code{mesh} into wedges
and hexahedra (or, with no 2-D cells, its lines into quads) along
\code{direction}, the whole displacement, in \code{num_layers} layers spaced
uniformly, by the geometric \code{ratio}, or at the explicit \code{offsets}
(fractions, increasing, ending at 1) (roadmap 3.1.2). \code{mio_revolve()}
turns them around the axis through \code{axis_point} along \code{axis} by
\code{angle} degrees in \code{num_layers} steps (3.1.3); points within
\code{tolerance} of the axis stay put and the cells around them collapse. The
result holds the swept cells only, with \code{Side} regions \code{bottom},
\code{top} (\code{start}, \code{end} for a partial turn) and \code{side}.
}
