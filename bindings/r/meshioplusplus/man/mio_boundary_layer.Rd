\name{mio_boundary_layer}
\alias{mio_boundary_layer}
\title{Grow boundary layers from a surface}
\usage{
mio_boundary_layer(mesh, first_height, num_layers = 3, growth = 1.2,
  reverse = FALSE, collision_fraction = 0.45, smooth_iterations = 5,
  max_shrink = 6, collapse = "points", normal_smooth_iterations = 0,
  front_check = TRUE)
}
\arguments{
\item{mesh}{A triangle and quad surface.}
\item{first_height}{The first layer's thickness.}
\item{num_layers, growth}{How many layers, and each over the one before.}
\item{reverse}{Grow against the faces' normals.}
\item{collision_fraction}{Of the distance a point's ray hits, in (0, 1].}
\item{smooth_iterations, max_shrink}{Smoothing passes; halvings before a point collapses.}
\item{collapse}{How a point whose cells fold gives way: \code{"points"} or \code{"layers"}.}
\item{normal_smooth_iterations}{Laplacian passes over the growth directions (0 to 1000).}
\item{front_check}{Test the last layer's faces against each other and the surface.}
}
\value{
A list: \code{mesh}, \code{num_cells}, \code{num_wedges}, \code{num_hexahedra},
\code{num_pyramids}, \code{num_tetrahedra}, \code{num_shrunk}, \code{num_collapsed},
\code{num_not_visible}, \code{min_scaled_jacobian}, \code{min_thickness},
\code{max_thickness}, \code{cap_closed}, \code{num_fewer_layers} and
\code{num_front_crossings}.
}
\description{
Prismatic layers grow from the triangle and quad surface \code{mesh} (roadmap
6.9): wedges from triangles and hexahedra from quads along each point's most
normal normal (Aubry and Löhner 2009), the first \code{first_height} thick and
each next \code{growth} times the one before; stopped at
\code{collision_fraction} of the distance to the surface ahead; the points of a
cell that would fold shrunk and then collapsed into pyramids and tetrahedra
(\code{collapse = "points"}) or their top layer dropped (\code{collapse =
"layers"}). The wall and the cap come out as faces, with
\code{boundary_layer:layer} cell data.
}
