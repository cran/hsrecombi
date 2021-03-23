#' @title Description of the targetregion data set
#' @name targetregion
#' @docType data
#' @description The data set contains sire haplotypes, assignment of progeny to
#'   sire, progeny genotypes and physical map information in a target region
#' \describe{
#'   The raw data can be downloaded at the source given below. Then,
#'   executing the following R code leads to the data provided in
#'   \code{targetregion.RData}.
#'   \item{\code{hapSire}}{matrix of sire haplotypes of each sire; 2 lines per
#'    sire; 1. column contains sireID}
#'   \item{\code{daughterSire}}{vector of sire ID for each progeny}
#'   \item{\code{genotype.chr}}{matrix of progeny genotypes}
#'   \item{\code{map.chr}}{SNP marker map in target region}
#' }
#' @source The data are available at RADAR
#'   \doi{10.22000/280}
#' @examples
#' \dontrun{
#' # download data from RADAR (requires about 1.4 GB)
#' url <- "https://www.radar-service.eu/radar-backend/archives/fqSPQoIvjtOGJlav/versions/1/content"
#' curl_download(url = url, 'tmp.tar')
#' untar('tmp.tar')
#' file.remove('tmp.tar')
#' path <- '10.22000-280/data/dataset'
#' ## list of haplotypes of sires for each chromosome
#' load(file.path(path, 'sire_haplotypes.RData'))
#' ## assign progeny to sire
#' daughterSire <- read.table(file.path(path, 'assign_to_family.txt'))[, 1]
#' ## progeny genotypes
#' X <- as.matrix(read.table(file.path(path, 'XFam-ARS.txt')))
#' ## physical and approximated genetic map
#' map <- read.table(file.path(path, 'map50K_ARS_reordered.txt'), header = T)
#' ## select target region
#' chr <- 1
#' window <- 301:600
#' ## map information of target region
#' map.chr <- map[map$Chr == chr, ][window, 1:5]
#' ## matrix of sire haplotypes in target region
#' hapSire <- rlist::list.rbind(haps[[chr]])
#' sireID <- 1:length(unique(daughterSire))
#' hapSire <- cbind(rep(sireID, each = 2), hapSire[, window])
#' ## matrix of progeny genotypes
#' genotype.chr <- X[, map.chr$SNP]
#' }
#' @importFrom rlist list.rbind
#' @importFrom curl curl_download
NULL
