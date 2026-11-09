\name{mio_fit_cad_iterate}
\alias{mio_fit_cad_iterate}
\title{Segment, fit and build a solid in rounds towards a target}
\usage{
mio_fit_cad_iterate(
  scan,
  step_path = NULL,
  part_class = "any",
  target_within = 0,
  target_explained = 0,
  segment_tolerance = 0,
  fit_tolerance = 0,
  max_rounds = 4,
  worst_faces = 3,
  refine_factor = 0.25,
  backend = "auto"
)
}
\arguments{
\item{scan}{A closed scan.}

\item{step_path}{Where the best solid goes as STEP (\code{NULL}: nowhere).}

\item{part_class}{\code{"any"}, \code{"analytic"}, \code{"cast"} or \code{"freeform"}.}

\item{target_within, target_explained}{The targets (\code{0}: the class's; off).}

\item{segment_tolerance, fit_tolerance}{Relative tolerances (\code{0}: the class's).}

\item{max_rounds, worst_faces, refine_factor}{The driver's budget.}

\item{backend}{\code{"auto"}, \code{"occt"} or \code{"native"}.}
}
\value{
A list: \code{mesh} (the best solid), \code{num_rounds}, \code{best_round}, \code{target_met}, \code{stop_reason}, \code{target_within}, \code{share_within}, \code{explained}, \code{hausdorff_max}, \code{num_faces}, \code{round_share_within} and \code{round_explained}.
}
\description{
Roadmap 7.7.2: each round segments, fits and builds; it stops at the target, at no improvement or after \code{max_rounds}, else re-segments the worst faces' patches finer. See \code{doc/fit_cad_loop.md}.
}
