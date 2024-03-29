% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{karlin}
\alias{karlin}
\title{Liberman and Karlin's genetic map function}
\usage{
karlin(N, x, inverse = F)
}
\arguments{
\item{N}{parameter (positive integer) required by the binomial model to
assess the count (of crossover) distribution; \code{N = 1} corresponds to
Morgan's map function}

\item{x}{vector of recombination rates}

\item{inverse}{logical, if FALSE recombination rate is mapped to Morgan unit,
if TRUE Morgan unit is mapped to recombination rate (default is FALSE)}
}
\value{
vector of genetic positions in Morgan units
}
\description{
Calculation of genetic distances from recombination rates given
  a parameter
}
\examples{
  karlin(2, seq(0, 0.5, 0.01))
}
\references{
Liberman, U. & Karlin, S. (1984) Theoretical models of genetic
  map functions. Theor Popul Biol 25:331-346.
}
