\name{mio_register_landmarks}
\alias{mio_register_landmarks}
\title{Landmark registration}
\usage{
mio_register_landmarks(source, target, weights = NULL)
}
\arguments{
\item{source, target}{\code{n x 3} matrices of picked points, \code{n >= 3}, the sources not collinear.}

\item{weights}{\code{NULL}, or \code{n} non-negative weights.}
}
\value{
A list of \code{matrix} (4 x 4), \code{residuals} and \code{rms}.
}
\description{
The rigid motion laying picked points onto their targets: the weighted Horn solve in closed form (roadmap 7.3.2). See \code{doc/registration.md}.
}
