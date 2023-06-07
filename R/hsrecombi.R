#' @title Estimation of recombination rate and maternal LD
#' @name hsrecombi
#' @description Wrapper function for estimating recombination rate and maternal
#'   linkage disequilibrium between intra-chromosomal SNP pairs by calling EM
#'   algorithm
#' @details Paternal recombination rate and maternal linkage disequilibrium (LD)
#'   are estimated for pairs of biallelic markers (such as single nucleotide
#'   polymorphisms; SNPs) from progeny genotypes and sire haplotypes. At least
#'   one sire has to be double heterozygous at the investigated pairs of SNPs.
#'   All progeny are merged in two genomic families: (1) coupling phase family
#'   if sires are double heterozygous 0-0/1-1 and (2) repulsion phase family if
#'   sires are double heterozygous 0-1/1-0.
#'   So far it is recommended processing the chromosomes separately.
#'   If maternal half-sib families are used, the roles of sire/dam are swapped.
#'   Multiple families can be considered.
#' @param hap list (LEN 2) of lists
#' \describe{
#'   \item{famID}{list (LEN number of sires) of vectors (LEN n.progeny) of
#'    progeny indices relating to lines in genotype matrix}
#'   \item{sireHap}{list (LEN number of sires) of matrices (DIM 2 x p)
#'    of sire haplotypes (0, 1) on investigated chromosome}
#' }
#' @param genotype.chr matrix (DIM n x p) of all progeny genotypes (0, 1, 2) on
#'   a chromosome with p SNPs; 9 indicates missing genotype
#' @param exclude vector (LEN < p) of SNP IDs (for filtering column names of
#'   \code{genotype.chr)} to be excluded from analysis (default NULL)
#' @param only.adj logical; if \code{TRUE}, recombination rate is calculated
#'   only between neighbouring markers
#' @param prec scalar; precision of estimation
#' @return list (LEN p - 1) of data.frames; for each SNP, parameters are
#'   estimated with all following SNPs; two solutions (prefix sln1 and sln2) are
#'   obtained for two runs of the EM algorithm
#' \describe{
#'   \item{\code{SNP1}}{ID of 1. SNP}
#'   \item{\code{SNP2}}{ID of 2. SNP}
#'   \item{\code{D}}{maternal LD}
#'   \item{\code{fAA}}{frequency of maternal haplotype 1-1}
#'   \item{\code{fAB}}{frequency of maternal haplotype 1-0}
#'   \item{\code{fBA}}{frequency of maternal haplotype 0-1}
#'   \item{\code{fBB}}{frequency of maternal haplotype 0-0}
#'   \item{\code{p1}}{Maternal allele frequency (allele 1) at \code{SNP1}}
#'   \item{\code{p2}}{Maternal allele frequency (allele 1) at \code{SNP2}}
#'   \item{\code{nfam1}}{size of genomic family 1}
#'   \item{\code{nfam2}}{size of genomic family 2}
#'   \item{\code{error}}{0 if computations were without error; 1 if EM algorithm
#'     did not converge}
#'   \item{\code{iteration}}{number of EM iterations}
#'   \item{\code{theta}}{paternal recombination rate}
#'   \item{\code{r2}}{\eqn{r^2} of maternal LD}
#'   \item{\code{logL}}{value of log likelihood function}
#'   \item{\code{unimodal}}{1 if likelihood is unimodal; 0 if likelihood is
#'     bimodal}
#'   \item{\code{critical}}{0 if parameter estimates are unique; 1 if parameter
#'     estimates at both solutions are valid, then decision process follows in
#'     post-processing function "editraw"}
#' }
#'   Afterwards, solutions are compared and processed with function
#'   \code{editraw}, yielding the final estimates for each valid pair of SNPs.
#' @examples
#'   ### test data
#'   data(targetregion)
#'   ### make list for paternal half-sib families
#'   hap <- makehaplist(daughterSire, hapSire)
#'   ### parameter estimates on a chromosome
#'   res <- hsrecombi(hap, genotype.chr)
#'   ### post-processing to achieve final and valid set of estimates
#'   final <- editraw(res, map.chr)
#' @references
#'   Hampel, A., Teuscher, F., Gomez-Raya, L., Doschoris, M. & Wittenburg, D.
#'    (2018) Estimation of recombination rate and maternal linkage
#'    disequilibrium in half-sibs. Frontiers in Genetics 9:186.
#'    \doi{10.3389/fgene.2018.00186}
#'
#'   Gomez-Raya, L. (2012) Maximum likelihood estimation of linkage
#'    disequilibrium in half-sib families. Genetics 191:195-213.
#' @importFrom Rcpp evalCpp
#' @importFrom rlist list.rbind
#' @useDynLib hsrecombi
#' @export
hsrecombi <- function(hap, genotype.chr, exclude = NULL, only.adj = FALSE, prec = 1e-6){

  if(ncol(genotype.chr) != ncol(hap$sireHap[[1]])) stop('ERROR: inconsistency in number of SNPs')
  if(!all(genotype.chr %in% c(0, 1, 2, 9))) stop('ERROR: coding of genotypes must be 0, 1, 2, 9')

  coln <- colnames(genotype.chr)
  snp.chr <- which(!(coln %in% exclude))

  ls <- lapply(1:(length(snp.chr) - 1), function(j){

    if(j %% 100 == 1) message(paste('Processing SNP', j))
    q <- ifelse(!only.adj, length(snp.chr), j + 1)

    out <- lapply((j + 1):q, function(i){

      ## set-up genomic families for snp pairs
      GenFam1 <- GenFam2 <- c()
      for (l in 1:length(hap$sireHap)){

        if((hap$sireHap[[l]][1, i] != hap$sireHap[[l]][2, i]) && (hap$sireHap[[l]][1, j] != hap$sireHap[[l]][2, j])){
          if((hap$sireHap[[l]][1, j] == 0 & hap$sireHap[[l]][1, i] == 0) | (hap$sireHap[[l]][1, j] == 1 & hap$sireHap[[l]][1, i] == 1)){
            GenFam1 <- rbind(GenFam1, genotype.chr[hap$famID[[l]], c(j, i)])

          } else GenFam2 <- rbind(GenFam2, genotype.chr[hap$famID[[l]], c(j, i)])
        }
      }

      if(!is.null(GenFam1) | !is.null(GenFam2)){
        ## start values 1. run
        theta.start <- prec
        start <- startvalue(GenFam1, GenFam2, 0)
        ## 1.run
        sln1 <- LDHScpp(GenFam1, GenFam2, start$fAA.start, start$fAB.start, start$fBA.start, theta.start, F, prec)

        ## start values 2. run
        Ds <- sln1$D + (1 - 2 * start$p1) * (1 - 2 * start$p2) / 4
        theta.start <- (1 - 4 * Ds) / 2
        Dd <- (1 - 2 * sln1$theta) / 4 - (1 - 2 * start$p1) * (1 - 2 * start$p2) / 4
        fAA.start <- Dd + start$p1 * start$p2
        fAB.start <- -Dd + start$p1 * (1 - start$p2)
        fBA.start <- -Dd + (1 - start$p1) * start$p2

        ## 2. run if second mode exists
        if((!(theta.start > 1 | theta.start < 0)) & ((Dd >= start$L1) & (Dd <= start$L2))){
          sln2 <- LDHScpp(GenFam1, GenFam2, fAA.start, fAB.start, fBA.start, theta.start, F, prec)
          unimodal <- 0
        } else {
          sln2 <- sln1
          unimodal <- 1
        }

        c(SNP1 = snp.chr[j], SNP2 = snp.chr[i], sln1 = unlist(sln1), sln2 = unlist(sln2),
                            unimodal = unimodal, critical = start$critical)
      } else NULL

    })
    out <- data.frame(rlist::list.rbind(out))
    out$SNP1 <- coln[out$SNP1]
    out$SNP2 <- coln[out$SNP2]
    out
  })
  return(ls)
}



#' @title Estimation of genetic position
#' @name geneticPosition
#' @description Estimation of genetic positions (in centi Morgan)
#' @details Smoothing of recombination rates (theta) <= 0.05 via quadratic
#'   optimization provides an approximation of genetic distances (in Morgan)
#'   between SNPs. The cumulative sum * 100 yields the genetic positions in cM.
#'
#'   The minimization problem \code{(theta - D d)^2} is solved s.t. d > 0 where
#'   d is the vector of genetic distances between adjacent markers but theta is
#'   not restricted to adjacent markers. The incidence matrix D contains 1's for
#'   those intervals contributing to the total distance relevant for each theta.
#'
#'   Estimates of theta = 1e-6 are neglected as these values coincide with start
#'   values and indicate that (because of a very flat likelihood surface) no
#'   meaningful estimate of recombination rate has been obtained.
#' @param final table of results produced by \code{editraw} with pairwise
#'   estimates of recombination rate between p SNPs within chromosome; minimum
#'   required data frame with columns \code{SNP1}, \code{SNP2} and \code{theta}
#' @param map1 data.frame containing information on physical map, at least:
#' \describe{
#'  \item{\code{SNP}}{SNP ID}
#'  \item{\code{locus_Mb}}{physical position in Mbp of SNP on chromosomes}
#'  \item{\code{Chr}}{chromosome of SNP}
#' }
#' @param exclude optional vector (LEN < p) of SNP IDs to be excluded (e.g.,
#'   candidates of misplaced SNPs; default NULL)
#' @param threshold optional value; recombination rates <= threshold are
#'   considered for smoothing approach assuming theta ~ Morgan (default 0.05)
#' @return list (LEN 2)
#' \describe{
#'   \item{gen.cM}{vector (LEN p) of genetic positions of SNPs (in cM)}
#'   \item{gen.Mb}{vector (LEN p) of physical positions of SNPs (in Mbp)}
#' }
#' @references
#'  Qanbari, S. & Wittenburg, D. (2020) Male recombination map of the autosomal
#'  genome in German Holstein. Genetics Selection Evolution 52:73.
#'  \doi{10.1186/s12711-020-00593-z}
#'
#' @examples
#'   ### test data
#'   data(targetregion)
#'   ### make list for paternal half-sib families
#'   hap <- makehaplist(daughterSire, hapSire)
#'   ### parameter estimates on a chromosome
#'   res <- hsrecombi(hap, genotype.chr)
#'   ### post-processing to achieve final and valid set of estimates
#'   final <- editraw(res, map.chr)
#'   ### approximation of genetic positions
#'   pos <- geneticPosition(final, map.chr)
#' @import Matrix
#' @importFrom rlist list.select
#' @importFrom quadprog solve.QP
#' @export
geneticPosition <- function(final, map1, exclude = NULL, threshold = 0.05){
  part <- final[(final$theta <= threshold) & (final$theta >= 1e-5) &
                  !(final$SNP1 %in% exclude) & !(final$SNP2 %in% exclude), ]
  part$ID1 <- match(part$SNP1, map1$SNP)
  part$ID2 <- match(part$SNP2, map1$SNP)
  mini <- min(part$ID1)
  p <- max(part$ID2) - mini + 1

  ## quadratic optimization min(theta - D d)² s.t. d > 0
  # set up D matrix
  out <-  lapply(1:nrow(part), function(i){
    vec.j <- (part$ID1[i] - mini + 1):(part$ID2[i] - mini)
    vec.i <- rep(i, length(vec.j))
    list(vec.i = vec.i, vec.j = vec.j)
  })
  vec.i <- unlist(list.select(out, vec.i))
  vec.j <- unlist(list.select(out, vec.j))
  D <- sparseMatrix(i = vec.i, j = vec.j, x = 1)

  # components for quadprog, make D p.d.
  dvec <- crossprod(D, part$theta)
  Amat <- diag(nrow = p - 1, x = 1)
  Dmat <- crossprod(D) + diag(nrow = p - 1, x = 1e-6)

  sln <- solve.QP(Dmat, dvec, Amat)
  # due to numerics in optimization approach (although constraints have been defined)
  sln$solution[sln$solution < 0] <- 0
  gen.cM <- c(0, cumsum(sln$solution)) * 100

  gen <- rep(NA, nrow(map1))
  id <- seq(from = min(part$ID1), to = max(part$ID2))
  gen[id] <- gen.cM
  names(gen) <- map1$SNP
  id <- unique(c(part$SNP1, part$SNP2))
  gen[!(names(gen) %in% id)] <- NA

  return(list(pos.cM = gen, pos.Mb = map1$locus_Mb))
}

