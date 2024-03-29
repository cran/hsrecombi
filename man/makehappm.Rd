% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{makehappm}
\alias{makehappm}
\title{Make list of imputed haplotypes and estimate recombination rate}
\usage{
makehappm(sireID, daughterSire, genotype.chr, nmin = 30, exclude = NULL)
}
\arguments{
\item{sireID}{vector (LEN N) of IDs of all sires}

\item{daughterSire}{vector (LEN n) of sire ID for each progeny}

\item{genotype.chr}{matrix (DIM n x p) of progeny genotypes (0, 1, 2) on a
single chromosome with p SNPs; 9 indicates missing genotype}

\item{nmin}{scalar, minimum required number of progeny for proper imputation,
default 30}

\item{exclude}{vector (LEN < p) of SNP IDs (for filtering column names of
\code{genotype.chr)} to be excluded from analysis}
}
\value{
list (LEN 2) of lists. For each sire:
\describe{
  \item{\code{famID}}{list (LEN N) of vectors (LEN n.progeny) of progeny
   indices relating to lines in genotype matrix}
  \item{\code{sireHap}}{list (LEN N) of matrices (DIM 2 x p) of sire
   haplotypes (0, 1) on investigated chromosome}
 \item{probRec}{vector (LEN p - 1) of proportion of recombinant progeny over
   all families between adjacent SNPs}
 \item{numberRec}{list (LEN N) of vectors (LEN n.progeny) of number of
   recombination events per animal}
 \item{gen}{vector (LEN p) of genetic positions of SNPs (in cM)}
}
}
\description{
List of sire haplotypes is set up in the format required for
  hsrecombi. Sire haplotypes are imputed from progeny genotypes using R
  package \code{hsphase}. Furthermore, recombination rate estimates between
  adjacent SNPs from hsphase are reported.
}
\examples{
  data(targetregion)
  hap <- makehappm(unique(daughterSire), daughterSire, genotype.chr, exclude = paste0('V', 301:310))
}
\references{
Ferdosi, M., Kinghorn, B., van der Werf, J., Lee, S. & Gondro, C. (2014)
  hsphase: an R package for pedigree reconstruction, detection of
  recombination events, phasing and imputation of half-sib family groups
  BMC Bioinformatics 15:172.
  \url{https://CRAN.R-project.org/package=hsphase}
}
