## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  tidy.opts = list(width.cutoff = 100), tidy = TRUE
)

## ----setup--------------------------------------------------------------------
library(hsrecombi) 
library(AlphaSimR)
library(doParallel)
library(ggplot2)
library(rlist)
Sys.time()

## ----global-------------------------------------------------------------------
# number of chromosomes
nchr <- 2
# Number of simulated SNPs (e.g., p = 1500 resembles an average bovine chromosome)
p <- 150
# Number of simulated progeny in each half-sib family 
n <- 1000
# Directory for (simulated) data
path.dat <- 'data'
dir.create(path.dat, showWarnings = FALSE)
# Directory for output
path.res <- 'results'
dir.create(path.res, showWarnings = FALSE)
# Number of computing clusters allocated
nclust <- 2

## ----alphasim-----------------------------------------------------------------
founderPop <- runMacs2(nInd = 1000, nChr = nchr, segSites = p)
SP <- SimParam$new(founderPop)
SP$setSexes("yes_sys")
# Enable tracing location of recombination events
SP$setTrackRec(TRUE)  

pop <- newPop(founderPop)
N <- 10; ntotal <- N * n
my_pop <- selectCross(pop = pop, nFemale = 500, nMale = N, use = "rand", nCrosses = ntotal)

probRec <- list()
for(chr in 1:nchr){
  co.pat <- matrix(0, ncol = p, nrow = ntotal)
  for(i in 1:ntotal){
    if(nrow(SP$recHist[[1000 + i]][[chr]][[2]]) > 1){
      # 1. line contains 1 1 by default
      loci <- SP$recHist[[1000 + i]][[chr]][[2]][-1, 2] 
      co.pat[i, loci] <- 1 
    }
  }
  probRec[[chr]] <- colMeans(co.pat)
}

save(list = c('SP', 'founderPop', 'pop', 'my_pop', 'ntotal', 'probRec'), file = 'data/pop.Rdata')

## ----genetic-data-------------------------------------------------------------
PAT <- my_pop@father
rown <- paste(rep(unique(PAT), each = 2), c(1, 2), sep = '_')
H.pat <- pullSegSiteHaplo(pop)[rown, ]
X <- pullSegSiteGeno(my_pop)
# Physical position of markers in Mbp
map.snp <- lapply(founderPop@genMap, function(z) z * 100)

## ----plink-format-------------------------------------------------------------
map <- data.frame(Chr = rep(1:length(map.snp), unlist(lapply(map.snp, length))), 
                  Name = paste0('SNP', 1:length(unlist(map.snp))),
                  locus_Mb = unlist(map.snp), 
                  locus_bp = unlist(map.snp) * 1e+6)

colnames(X) <- map$Name
FID <- 'FAM001'
IID <- my_pop@id
MAT <- my_pop@mother
SEX <- 2
PHENOTYPE <- -9

for(chr in 1:nchr){
  write.table(map[map$Chr == chr, ], file.path(path.dat, paste0('map_chr', chr, '.map')), 
              col.names = F, row.names = F, quote = F)
  write.table(cbind(FID, IID, PAT, MAT, SEX, PHENOTYPE, X[, map$Chr == chr]), 
              file.path(path.dat, paste0('hsphase_input_chr', chr, '.raw')), col.names = T, row.names = F, quote = F) 
}

## ----parallel-computing-------------------------------------------------------
cl <- makeCluster(nclust)
registerDoParallel(cl)

## ----recombination-rate-------------------------------------------------------
out <- foreach(chr = 1:nchr, .packages = 'hsrecombi') %dopar% {
  
  # 1: Physical  map
  map <- read.table(file.path(path.dat, paste0('map_chr', chr, '.map')), col.names = c('Chr', 'SNP', 'locus_Mb', 'locus_bp'))
 
  # 2: Genotype matrix
  genomatrix <- data.table::fread(file.path(path.dat, paste0('hsphase_input_chr', chr, '.raw')))
  X <- as.matrix(genomatrix[, -c(1:6)])
  X[is.na(X)] <- 9 # required for hsphase
  
  # 3: Assign daughters to sire IDs
  daughterSire <- genomatrix$PAT
  
  # 4: Estimate sire haplotypes and format data 
  hap <- makehappm(unique(daughterSire), daughterSire, X)
  save('hap', file = file.path(path.res, paste0('hsphase_output_chr', chr, '.Rdata')))
  
  # Check order and dimension
  io <- sapply(1:nrow(map), function(z){grepl(x = colnames(X)[z], pattern = map$SNP[z])})
  if(sum(io) != nrow(map)) stop("ERROR in dimension")
  
  # 5: Estimate recombination rates
  res <- hsrecombi(hap, X)
  final <- editraw(res, map)
  save(list = c('final', 'map'), file = file.path(path.res, paste0("Results_chr", chr, ".Rdata")))
  
  ifelse(nrow(final) > 0, 'OK', 'no result')
}
print(which(unlist(out) == 'OK'))

## ----misplaced----------------------------------------------------------------
# 6a: Filter SNPs with unusually large recombination rate to neighbouring (30) SNPs
excl <- foreach(chr = 1:nchr, .packages = 'hsrecombi') %dopar% {
  load(file.path(path.res, paste0("Results_chr", chr, ".Rdata")))
  checkCandidates(final, map)
}

# 6b: Heatmap plot of recombination rates for visual verification, e.g.:
chr <- 2
load(file.path(path.res, paste0("Results_chr", chr, ".Rdata")))
cand <- excl[[chr]][1]
win <- match(cand, map$SNP) + (-100:100)
win <- win[(win >= 1) & (win <= nrow(map))]

target <- final[(final$SNP1 %in% map$SNP[win]) & (final$SNP2 %in% map$SNP[win]), ]
target$SNP1 <- match(target$SNP1, map$SNP)
target$SNP2 <- match(target$SNP2, map$SNP)

ggplot(data = target, aes(SNP2, SNP1, fill = theta)) + 
  geom_tile() +
  xlab("Locus 2") +
  ylab("Locus 1") +
  coord_equal() + 
  scale_y_continuous(trans = "reverse") +
  theme(panel.background = element_blank(), 
        panel.grid.major = element_line(colour = "grey", linewidth = 0.1),
        panel.grid.minor = element_line(colour = "grey")) +
  theme(text = element_text(size = 18)) +
  scale_fill_gradientn(colours = c('yellow', 'red'), limits = c(0, 1+1e-10), na.value = 'white')

# -> nothing conspicious
excl[[1]] <- excl[[2]] <- NA

## ----genetic-position---------------------------------------------------------
pos <- foreach(chr = 1:nchr, .packages = 'hsrecombi') %dopar% {
  load(file.path(path.res, paste0("Results_chr", chr, ".Rdata")))
  geneticPosition(final, map, exclude = excl[[chr]])
}

## ----stop-parallel------------------------------------------------------------
stopCluster(cl)

## ----plot-1-------------------------------------------------------------------
data <- data.frame(Chr = rep(1:length(pos), times = unlist(list.select(pos, length(pos.Mb)))), 
                   Mbp = unlist(list.select(pos, pos.Mb)), cM = unlist(list.select(pos, pos.cM)))
ggplot(data, aes(x = Mbp, y = cM)) + geom_point(na.rm = T) + facet_grid(Chr ~ .)

## ----plot-2-------------------------------------------------------------------
for (chr in 1:nchr) {
  load(file.path(path.res, paste0("Results_chr", chr, ".Rdata")))
  final$theta[(final$SNP1 %in% excl[[chr]]) | (final$SNP2 %in% excl[[chr]])] <- NA
  final$dist_M <- (pos[[chr]]$pos.cM[final$SNP2] - pos[[chr]]$pos.cM[final$SNP1]) / 100
  out <- bestmapfun(theta = final$theta, dist_M = final$dist_M)
  plot(final$dist_M, final$theta, xlab = 'genetic distance (Morgan)', ylab = 'recombination rate', col = 8)
  x <- seq(0, 0.5, 0.01); y <- rao(out$mixing, x); points(y, x, type = 'l', lwd = 2)
  legend('bottomright', legend = paste0('map function (p=', round(out$mixing, 3), ')'), lwd = 2, lty = 1, bty = 'n')
}

## ----selection----------------------------------------------------------------
chr <- 2

## ----check-simulated----------------------------------------------------------
load(file.path(path.dat, 'pop.Rdata'))
sim.cM <- cumsum(probRec[[chr]]) * 100

## ----check-deterministic------------------------------------------------------
load(file.path(path.res, paste0('hsphase_output_chr', chr, '.Rdata')))
hsphase.cM <- hap$gen

## ----plot-3-------------------------------------------------------------------
par(mar=c(5.1, 4.1, 4.1, 9.1), xpd = TRUE)
plot(pos[[chr]]$pos.Mb, pos[[chr]]$pos.cM, xlab = 'physical position (Mbp)', ylab = 'genetic position (cM)', 
     ylim = range(c(sim.cM, hsphase.cM, pos[[chr]]$pos.cM), na.rm = T), pch = 20)
points(pos[[chr]]$pos.Mb, sim.cM, pch = 20, col = 4)
points(pos[[chr]]$pos.Mb, hsphase.cM, pch = 20, col = 8)
legend('topleft', inset=c(1.01,0), legend = c('simulated', 'likelihood-based', 'deterministic'), pch = 20, col = c(4, 1, 8), bty = 'n')

## ----cleanup------------------------------------------------------------------
unlink(path.dat, recursive = TRUE)
unlink(path.res, recursive = TRUE)
Sys.time()

