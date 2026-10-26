\name{mio_mesh_cad_surface}
\alias{mio_mesh_cad_surface}
\title{Triangulate a CAD model's faces}
\usage{
mio_mesh_cad_surface(info, target_size, max_chordal_deviation = 0,
  min_points_per_edge = 0, tolerance = 0, record_face_ids = TRUE, edges = TRUE,
  match_unreferenced_trims = TRUE, strict = FALSE, max_points = 50000000,
  min_angle = 0, smooth_iterations = 0, smooth_method = c("laplacian", "odt"),
  order = 1, grading = 0, curvature_deviation = 0, cells_across_gaps = 0,
  pole_rings = 0, face_sizes = NULL, edge_sizes = NULL, background = NULL,
  background_array = "size", metric = c("mean", "varying"),
  virtual_edges = FALSE, reparam_condition_number = 0,
  model_tolerances = FALSE, virtual_faces = FALSE, lop_corners = FALSE,
  metric_background = NULL, metric_array = "", metric_error = 0,
  untangle = FALSE, quads = c("none", "dominant", "all"),
  quad_pairs = c("matching", "greedy"))
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
\item{min_angle}{Refine until every triangle's smallest angle on the surface is
at least this (degrees, under 60); 0 keeps the unrefined mesh.}
\item{smooth_iterations}{Passes of guarded smoothing; 0 is none.}
\item{smooth_method}{\code{"laplacian"} or \code{"odt"}.}
\item{order}{1 linear triangles; 2 or 3 curved \code{triangle6} or
\code{triangle10} elements (and \code{line3} or \code{line4} edges) whose nodes
lie on the surfaces and edges; an element whose Jacobian is not positive
everywhere is kept straight-sided and counted in \code{num_curved_straightened}.}
\item{grading}{The size's largest growth, \code{|grad h| <= grading - 1}; 0 is
1.5. Acts only with \code{curvature_deviation} or \code{cells_across_gaps} on.}
\item{curvature_deviation}{Size each region so a chord sits at most this far
from the surface; 0 is off.}
\item{cells_across_gaps}{Size each region so at least this many elements span
its local feature size (a wall's thickness, a slot's width); 0 is off.}
\item{pole_rings}{Rings of points seeded around each pole, so it is no longer a
fan's apex; 0 keeps the fan, at most 8.}
\item{face_sizes}{Sizes given per face: a numeric vector named by 1-based face
indices, or a two-column matrix of (index, size); every point of the face is at
most that size (v17.19.0).}
\item{edge_sizes}{The same per edge (1-based edge indices).}
\item{background}{A mesh whose one-component point data \code{background_array}
is the size at its points, read at its nearest surface point; \code{NULL} for
none.}
\item{background_array}{The name of that point data array.}
\item{metric}{\code{"mean"} meshes each face in its mean metric;
\code{"varying"} follows the surface's metric as it varies over the face.}
\item{virtual_edges}{Sample chains of short edges meeting smoothly as one; the
chains are reported in \code{virtual_edges}.}
\item{reparam_condition_number}{Mesh a face whose metric's condition number over
its coarse mesh passes this in a flattening of that mesh instead of its
parameters (v17.20.0); 0 is off, otherwise above 1; order 1 only.}
\item{model_tolerances}{Accept the gaps the model itself declares (an edge's
\code{"tolerance"} and \code{"end_tolerances"} in its extra JSON).}
\item{virtual_faces}{Mesh a face narrower than \code{target_size} with the
neighbour it meets tangentially along the most edge; order 1 only.}
\item{lop_corners}{Lop acute corners (Pav and Walkington 2005) in place of
concentric shells; only with refinement (\code{min_angle} or sizes).}
\item{metric_background}{A mesh whose 9-component point data
\code{metric_array} is a field's Hessian: an anisotropic metric from it, with
the target interpolation error \code{metric_error}, sizes the triangles in every
direction (roadmap 6.10.2); \code{NULL} for none; order 1 only.}
\item{metric_array}{The name of that Hessian array (e.g. \code{"f:hessian"}).}
\item{metric_error}{The target interpolation error (positive).}
\item{untangle}{With \code{order} 2 or 3, move a face's own high-order nodes on
the surface to untangle a folded curved element before it is straightened
(roadmap 6.6.3; Toulorge et al. 2013).}
\item{quads}{\code{"none"}, \code{"dominant"} (pairs of triangles recombined
into quads where they make a good one; order 1 only) or \code{"all"} (the
quad-dominant mesh of twice the size split into quads only; \code{order = 2}
writes quad9) (roadmap 6.8).}
\item{quad_pairs}{How triangles pair into quads: \code{"matching"} (a maximum
matching of greatest total quality, Blossom-quad) or \code{"greedy"} (best pair
first) (roadmap 6.8.1).}
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
per face), \code{warnings} (a character vector), refinement's
\code{num_refinement_points}, \code{num_faces_below_bound}, \code{num_smoothed} and
\code{min_input_angle_deg}, \code{angle_histogram} (18 ten-degree bins),
\code{num_curved_straightened}, \code{face_size_limits} (per face, the 1-based
edges finer than its given size allows), \code{virtual_edges} (each chain's
1-based edges), \code{reparametrised} (a data frame: \code{face},
\code{condition_number}, \code{condition_number_disc}, \code{method} 1
conformal, 2 mean-value, 3 annulus, 0 declined, \code{reason}), \code{merged_faces}
(\code{sliver}, \code{host}, \code{width}, \code{merged}, \code{reason}),
\code{tolerance_uses} (\code{edge}, \code{declared}, \code{gap},
\code{where}), \code{lopped_corners} (\code{face}, \code{point},
\code{angle_deg}, \code{lfs0}, \code{lop_radius}, \code{turning_deg},
\code{stated_bound_deg}, \code{num_forced}, \code{min_forced_deg},
\code{within_two_radii}, \code{guaranteed}, \code{reason}), \code{num_untangled},
\code{num_quads}, \code{num_leftover_triangles}, \code{min_quad_scaled_jacobian}
and \code{corners}, a data frame of the corners under 60 degrees and poles (1-based
\code{face} and \code{point}, \code{u}, \code{v}, \code{x}, \code{y}, \code{z},
\code{angle_deg}, \code{min_forced_deg}, \code{reach}, \code{pole}).}
