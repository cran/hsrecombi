% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{haldane}
\alias{haldane}
\title{Haldane's genetic map function}
\usage{
haldane(x, inverse = F)
}
\arguments{
\item{x}{vector of recombination rates}

\item{inverse}{logical, if FALSE recombination rate is mapped to Morgan unit,
if TRUE Morgan unit is mapped to recombination rate (default is FALSE)}
}
\value{
vector of genetic positions in Morgan units
}
\description{
Calculation of genetic distances from recombination rates
}
\examples{
  haldane(seq(0, 0.5, 0.01))
}
\references{
Haldane JBS (1919) The combination of linkage values, and the
  calculation of distances between the loci of linked factors. J Genet 8:
  299-309.
}
