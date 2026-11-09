\name{mio_auto_align}
\alias{mio_auto_align}
\title{Pose a scan from itself alone}
\usage{
mio_auto_align(mesh, angle_tolerance = 10, subdivisions = 3, region = "",
  rotate_vector_data = FALSE)
}
\arguments{
\item{mesh}{A \code{mio_mesh} with surface cells.}

\item{angle_tolerance}{Degrees within which a normal agrees, in (0, 45).}

\item{subdivisions}{The Gauss-sphere icosahedron's splits (0 to 5).}

\item{region}{A voting cell region; \code{""} for all.}

\item{rotate_vector_data}{Rotate vector and tensor point data too.}
}
\value{
A list of \code{mesh}, \code{matrix}, \code{mode} (\code{"face"} or \code{"axis"}), \code{primary}, \code{secondary}, \code{share}, \code{face_share}, \code{axis_share} and \code{opening_area}.
}
\description{
On its largest flat face, or a cylindrical part on the end of its axis, from an area-weighted histogram of its facet normals (roadmap 7.3.5). See \code{doc/registration.md}.
}
