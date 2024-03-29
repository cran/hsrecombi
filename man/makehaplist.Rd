% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{makehaplist}
\alias{makehaplist}
\title{Make list of sire haplotypes}
\usage{
makehaplist(daughterSire, hapSire, nmin = 1)
}
\arguments{
\item{daughterSire}{vector (LEN n) of sire ID for each progeny}

\item{hapSire}{matrix (DIM 2N x p + 1) of sire haplotype at p SNPs; 2 lines
per sire, 1. columns contains sire ID}

\item{nmin}{scalar, minimum number of progeny required, default 1}
}
\value{
list (LEN 2) of lists. For each sire:
\describe{
  \item{\code{famID}}{list (LEN N) of vectors (LEN n.progeny) of progeny
   indices relating to lines in genotype matrix}
  \item{\code{sireHap}}{list (LEN N) of matrices (DIM 2 x p) of sire
   haplotypes (0, 1) on investigated chromosome}
}
}
\description{
List of sire haplotypes is set up in the format required for
  hsrecombi. Haplotypes (obtained by external software) are provided.
}
\examples{
  data(targetregion)
  hap <- makehaplist(daughterSire, hapSire)
}
