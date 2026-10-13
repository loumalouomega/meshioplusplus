\name{mio_brep_info}
\alias{mio_brep_info}
\title{Copy a B-rep side channel}
\usage{mio_brep_info(info)}
\arguments{\item{info}{A format-info handle from \code{mio_read_with_info()} of a CAD JSON (or STEP/IGES) file.}}
\description{Copy the trimmed NURBS faces, loops, edges and vertices of a Kratos CAD JSON handle (or a STEP/IGES one read through OpenCASCADE) into plain R values.}
\value{A list with \code{patches} (\code{degrees}, \code{knots} in the full
form, \code{control_points} as a \code{(3, points)} matrix, \code{weights});
\code{breps} (\code{id}, \code{name}, 1-based \code{faces}, \code{edges},
\code{vertices}, \code{extra}); \code{faces} (\code{id}, \code{name}, 1-based
\code{patch}, \code{swapped_normal}, \code{trimmed}, \code{rational},
\code{reduced_knots}, \code{control_point_ids}, \code{loops} of \code{type},
\code{trims} and \code{extra}, \code{embedded_edges}, \code{surface_extra},
\code{extra}), each trim \code{trim_index}, \code{direction}, \code{curve},
\code{extra}; \code{edges} (\code{id}, \code{name}, \code{curve} or
\code{NULL}, \code{topology}, \code{extra}); \code{vertices} (\code{id},
\code{name}, \code{topology}, \code{extra}), each topology entry \code{id},
\code{name}, \code{trim_index}, \code{relative_direction},
\code{local_coordinates}, \code{extra}; and \code{extra}, the file's other
members. A curve is \code{degree}, \code{knots}, \code{reduced_knots},
\code{rational}, \code{control_points}, \code{weights}, \code{ids},
\code{active_range}, \code{extra}. Ids are -1 where the file names the entity
instead; the \code{extra} members are compact JSON text. The copy remains
valid after releasing the handle; \code{mio_write_with_info()} writes the
B-rep back as CAD JSON from the handle itself.}
