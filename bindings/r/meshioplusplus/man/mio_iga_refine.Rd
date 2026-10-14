\name{mio_iga_refine}
\alias{mio_iga_refine}
\alias{mio_iga_sample}
\alias{mio_iga_evaluate}
\title{Refine, sample and evaluate a NURBS model}
\usage{
mio_iga_refine(info, subdivisions = 1L, elevate = 0L, sample_subdivisions = 0L,
  trim_samples = 16L, edges = TRUE, trim = TRUE)

mio_iga_sample(info, sample_subdivisions = 0L, trim_samples = 16L, edges = TRUE,
  trim = TRUE)

mio_iga_evaluate(info, patch, params, derivatives = 0L, trim_samples = 16L)
}
\arguments{
\item{info}{A format-info handle holding a NURBS model.}
\item{subdivisions, elevate}{One integer, or one per parametric direction.}
\item{sample_subdivisions}{Cells per knot span of the sampling; 0 is automatic.}
\item{trim_samples}{B-rep: points per knot span a trimming curve is sampled at.}
\item{edges, trim}{B-rep: sample the edge curves; keep only the cells inside
the trimming loops.}
\item{patch}{Patch (B-rep: face) to evaluate, 1-based.}
\item{params}{A \code{(d, N)} matrix of parametric coordinates, one column per point.}
\item{derivatives}{0, or 1 for the first derivatives and the surface normals.}
}
\description{The isogeometric layer (v16.43.0): act on the NURBS model a reader
kept in a format-info handle (an MFEM NURBS mesh, Kratos CAD JSON, STEP or IGES,
from \code{mio_read_with_info()}) without going through its sampling.

\code{mio_iga_refine()} refines every patch (every face surface of a B-rep)
exactly: degree elevation by \code{elevate}, then knot insertion splitting every
nonzero knot span in \code{subdivisions} parts (an MFEM model refines every
direction alike). The geometry, the fields and a trimmed face's trimming loops
do not change. \code{mio_iga_sample()} samples the model into Lagrange cells.
\code{mio_iga_evaluate()} evaluates patch \code{patch} at the columns of
\code{params}.}
\value{\code{mio_iga_refine()} and \code{mio_iga_sample()}: \code{list(mesh = , info = )}, the
sampled mesh and the model (refined, or with the sampling's cell maps), for
\code{mio_write_with_info()}. \code{mio_iga_evaluate()}: a list with \code{points}
(\code{(space_dimension, N)}), \code{fields} (one \code{(components, N)} matrix per
NURBS field), \code{derivatives} (\code{(space_dimension, d, N)}) and \code{normals}
(\code{(3, N)}), both \code{NULL} unless \code{derivatives = 1}, and \code{inside}
(a logical vector for a B-rep face, else \code{NULL}).}
