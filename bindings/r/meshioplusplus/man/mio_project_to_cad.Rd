\name{mio_project_to_cad}
\alias{mio_project_to_cad}
\title{Move a mesh's points back onto its CAD model}
\usage{
mio_project_to_cad(mesh, info, face_array = "cad:face", edge_array = "cad:edge",
  max_distance = 0, trim_samples = 16)
}
\arguments{
\item{mesh}{A surface mesh carrying \code{face_array} (and \code{edge_array}) cell data.}
\item{info}{A format-info handle holding a B-rep.}
\item{face_array, edge_array}{The cell data naming each cell's face and edge.}
\item{max_distance}{A point farther than this stays (0: no limit).}
\item{trim_samples}{Samples per knot span when testing a face's trims.}
}
\value{
A list: \code{mesh}, \code{num_projected}, \code{num_on_edges}, \code{num_corners},
\code{num_ambiguous}, \code{num_too_far}, \code{num_failed}, \code{num_outside_trim},
\code{max_move} and \code{outside_trim} (1-based points projected outside their
face's trims, kept, never clamped).
}
\description{
Every point of \code{mesh} goes onto the CAD face or edge its cells name
(\code{cad:face} on triangles and quads, \code{cad:edge} on lines, as
\code{\link{mio_mesh_cad_surface}} writes them) of the B-rep held by \code{info}
(roadmap 6.10.1): a point on one edge's lines onto that edge, on two or more (a
corner) onto the lowest, a point whose cells name one face onto its surface
(Newton point inversion); a point on two faces and no edge stays. Data and
regions pass through.
}
