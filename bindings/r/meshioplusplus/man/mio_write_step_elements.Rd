\name{mio_write_step_elements}
\alias{mio_write_step_elements}
\title{Fitted primitives as STEP elements}
\usage{
mio_write_step_elements(path, mesh, patches = numeric(0), segments = 64)
}
\arguments{
\item{path}{The \code{.step} file to create.}

\item{mesh}{A mesh carrying \code{segment:patch} and the \code{fit:patch_*} table.}

\item{patches}{Patch labels to write; empty writes every fitted patch.}

\item{segments}{Facets round a cylinder, cone or torus in the tessellation the B-rep writer is given.}
}
\value{
A list: \code{num_planes}, \code{num_spheres}, \code{num_cylinders}, \code{num_cones}, \code{num_tori}, \code{num_skipped} and \code{warnings}.
}
\description{
Write the fitted primitives of a mesh (\code{mio_segment()}, then \code{mio_fit_primitives()}) to a STEP file as exact geometry a CAD user can snap to: bounded planes, balls, lidded cylinders, frusta and tori (roadmap 7.6.5). See \code{doc/formats/step.md}.
}
