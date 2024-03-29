% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hsrecombi.R
\name{geneticPosition}
\alias{geneticPosition}
\title{Estimation of genetic position}
\usage{
geneticPosition(final, map1, exclude = NULL, threshold = 0.05)
}
\arguments{
\item{final}{table of results produced by \code{editraw} with pairwise
estimates of recombination rate between p SNPs within chromosome; minimum
required data frame with columns \code{SNP1}, \code{SNP2} and \code{theta}}

\item{map1}{data.frame containing information on physical map, at least:
\describe{
 \item{\code{SNP}}{SNP ID}
 \item{\code{locus_Mb}}{physical position in Mbp of SNP on chromosomes}
 \item{\code{Chr}}{chromosome of SNP}
}}

\item{exclude}{optional vector (LEN < p) of SNP IDs to be excluded (e.g.,
candidates of misplaced SNPs; default NULL)}

\item{threshold}{optional value; recombination rates <= threshold are
considered for smoothing approach assuming theta ~ Morgan (default 0.05)}
}
\value{
list (LEN 2)
\describe{
  \item{gen.cM}{vector (LEN p) of genetic positions of SNPs (in cM)}
  \item{gen.Mb}{vector (LEN p) of physical positions of SNPs (in Mbp)}
}
}
\description{
Estimation of genetic positions (in centi Morgan)
}
\details{
Smoothing of recombination rates (theta) <= 0.05 via quadratic
  optimization provides an approximation of genetic distances (in Morgan)
  between SNPs. The cumulative sum * 100 yields the genetic positions in cM.

  The minimization problem \code{(theta - D d)^2} is solved s.t. d > 0 where
  d is the vector of genetic distances between adjacent markers but theta is
  not restricted to adjacent markers. The incidence matrix D contains 1's for
  those intervals contributing to the total distance relevant for each theta.

  Estimates of theta = 1e-6 are neglected as these values coincide with start
  values and indicate that (because of a very flat likelihood surface) no
  meaningful estimate of recombination rate has been obtained.
}
\examples{
  ### test data
  data(targetregion)
  ### make list for paternal half-sib families
  hap <- makehaplist(daughterSire, hapSire)
  ### parameter estimates on a chromosome
  res <- hsrecombi(hap, genotype.chr)
  ### post-processing to achieve final and valid set of estimates
  final <- editraw(res, map.chr)
  ### approximation of genetic positions
  pos <- geneticPosition(final, map.chr)
}
\references{
Qanbari, S. & Wittenburg, D. (2020) Male recombination map of the autosomal
 genome in German Holstein. Genetics Selection Evolution 52:73.
 \doi{10.1186/s12711-020-00593-z}
}
