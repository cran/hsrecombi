% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{bestmapfun}
\alias{bestmapfun}
\title{Best fitting genetic-map function}
\usage{
bestmapfun(theta, dist_M)
}
\arguments{
\item{theta}{vector of recombination rates}

\item{dist_M}{vector of genetic positions}
}
\value{
list (LEN 2)
\describe{
  \item{mixing}{mixing parameter of system of genetic mapping functions}
  \item{mse}{minimum value of target function (theta - dist_M)^2}
}
}
\description{
Approximation of mixing parameter of system of map functions
}
\details{
The genetic mapping function that fits best to the genetic data
  (recombination rate and genetic distances) is obtained from Rao's system of
  genetic-map functions. The corresponding mixing parameter is estimated via
  1-dimensional constrained optimisation.
  See vignette for its application to estimated data.
}
\examples{
  theta <- seq(0, 0.5, 0.01)
  gendist <- -log(1 - 2 * theta) / 2
  bestmapfun(theta, gendist)
}
\references{
Rao, D.C., Morton, N.E., Lindsten, J., Hulten, M. & Yee, S (1977)
  A mapping function for man. Human Heredity 27: 99-104.
  \doi{10.1159/000152856}
}
