\name{mio_mesh_cad_surface}
\alias{mio_mesh_cad_surface}
\title{Triangulate a CAD model's faces}
\usage{
mio_mesh_cad_surface(info, target_size, max_chordal_deviation = 0,
  min_points_per_edge = 0, tolerance = 0, record_face_ids = TRUE, edges = TRUE,
  match_unreferenced_trims = TRUE, strict = FALSE, max_points = 50000000)
}
\arguments{
\item{info}{A format-info handle holding a B-rep.}
\item{target_size}{The target edge length (required, > 0).}
\item{max_chordal_deviation}{Sagitta cap of the boundary segments; 0 is off.}
\item{min_points_per_edge}{At least this many points inside every edge.}
\item{tolerance}{Geometric tolerance; 0 is 1e-7 of the model's size.}
\item{record_face_ids}{\code{cad:face} cell data and a \code{face:<k>} region per face.}
\item{edges}{A line block of the edges, with \code{cad:edge}.}
\item{match_unreferenced_trims}{Pair trims no edge names by their images.}
\item{strict}{Fail when any face cannot be meshed.}
\item{max_points}{The most points in the result.}
}
\description{Mesh generation from CAD (v17.13.0): a conforming triangle mesh of
the faces of the B-rep a format-info handle holds (Kratos CAD JSON, STEP or IGES,
from \code{mio_read_with_info()}) with edges of about \code{target_size}, with no
CAD kernel. Every edge is sampled once and shared by the faces that use it, so a
closed solid comes out watertight; each face is a constrained Delaunay
triangulation of its parameters lifted onto its surface. A face that cannot be
meshed is left out and reported, never meshed wrong; \code{strict = TRUE} makes
that an error.}
\value{A list: \code{mesh}, the counts (\code{num_points}, \code{num_triangles},
\code{num_edge_segments}, \code{num_faces}, \code{num_failed}, \code{num_edges},
\code{num_open_edges}, \code{num_non_manifold_edges}, \code{num_misoriented_edges},
\code{num_degenerate_dropped}), \code{tolerance}, the extreme angles and edge
lengths, \code{max_edge_gap}, \code{face_status} and \code{face_reason} (one entry
per face) and \code{warnings} (a character vector).}
