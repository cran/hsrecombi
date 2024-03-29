% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{rao inverse}
\alias{rao inverse}
\alias{rao.inv}
\title{Approximation to inverse of Rao's system of map functions}
\usage{
rao.inv(p, x)
}
\arguments{
\item{p}{mixing parameter (see details); \code{0 <= p <= 1}}

\item{x}{vector in Morgan units}
}
\value{
vector of recombination rates
}
\description{
Calculation of recombination rates from genetic distances given
  a mixing parameter
}
\details{
Mixing parameter \code{p=0} would match to Morgan, \code{p=0.25} to
  Carter, \code{p=0.5} to Kosambi and \code{p=1} to Haldane map function.
}
\examples{
  rao.inv(0.25, seq(0, 01, 0.1))
}
\references{
Rao, D.C., Morton, N.E., Lindsten, J., Hulten, M. & Yee, S (1977)
  A mapping function for man. Human Heredity 27: 99-104.
  \doi{10.1159/000152856}
}
