\name{mio_nurbs_info}
\alias{mio_nurbs_info}
\title{Copy an MFEM NURBS side channel}
\usage{mio_nurbs_info(info)}
\arguments{\item{info}{A format-info handle from \code{mio_read_with_info()} of an MFEM mesh.}}
\description{Copy the patches an \code{MFEM NURBS mesh} was sampled from into plain R values (v16.40.0).}
\value{A list with \code{patches}, each a list of \code{degrees}, \code{knots}
(one vector per parametric direction), \code{control_points} (a
\code{(space_dimension, points)} matrix, direction 0 fastest), \code{weights}
and \code{fields} (one \code{(components, points)} matrix of control values per
field); \code{field_names}; \code{topology} (\code{dimension},
\code{num_vertices}, \code{patch_vertices}, \code{patch_attributes},
\code{boundary_vertices}, \code{boundary_attributes} and \code{edges} as
\code{(knot vector, vertex, vertex)} columns), \code{NULL} when the source
states none; and the cell map \code{cell_patch}, \code{cell_span}
(\code{(3, cells)}) and \code{boundary_cell_patch}. Vertices, patches and knot
spans are 1-based; attributes and knot-vector ids stay as the file spelled
them. The copy remains valid after releasing the handle;
\code{mio_write_with_info()} writes the patches back from the handle itself.}
