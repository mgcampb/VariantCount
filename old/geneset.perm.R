t <- Sys.time()

#### Inputs
require(parallel)
tab <- read.delim("Hg19AllGenesRareVariantSummaryCounts.txt") # full table (huge!)
load("kegg.lengths.RData") # len = vector of unique lengths of kegg pathways
len <- unique(round(len,-1)) # save time by rounding to nearest 10
NP <- 10000
cpus <- detectCores()

#### Important variables
genes <- unique(tab[,1])
ids <- unique(tab[,2])
dx <- data.frame(id=ids, dx=tab[match(ids,tab[,2]),3])
n <- length(genes)
m <- length(ids)

#### Reorganize table into 2, one for hom alt counts and one for all other counts
tab1 <- matrix(0, n, m) # hom alt counts
tab2 <- matrix(0, n, m) # NOT hom alt (hom ref or het alt) counts
rownames(tab1) <- genes
rownames(tab2) <- genes
colnames(tab1) <- ids
colnames(tab2) <- ids
for(i in 1:m) {
  cat("Extracting sample", i, "/", m, ": elapsed time = ", Sys.time()-t,"\n")
  sel <- which(tab[,2]==ids[i])
  tab1[match(tab[sel,1], genes),i] <- tab[sel,6] 
  tab2[match(tab[sel,1], genes),i] <- tab[sel,4] + tab[sel,5]
}
save(tab1, tab2, dx, file="Hg19HomRareAndNotHomRareCountTables.RData")
rm(tab)

#### Define permutation function
perm <- function(sel) {
  case.ct1 <- sum(apply(tab2[sel,case], 2, sum)) # hom alt case
  case.ct2 <- sum(apply(tab1[sel,case], 2, sum)) # NOT hom alt case
  ctrl.ct1 <- sum(apply(tab2[sel,!case], 2, sum)) # hom alt ctrl
  ctrl.ct2 <- sum(apply(tab1[sel,!case], 2, sum)) # NOT hom alt ctrl
  p <- prop.test(matrix(c(case.ct1,ctrl.ct1,case.ct2,ctrl.ct2),2,2))$p.value
  if(is.na(p)) return(1) else return(p)
}

#### Parallelize
system("rm logfile")
cl <- makeCluster(cpus, outfile="logfile")
clusterExport(cl, c("tab1", "tab2", "case", "perm"))

#### Perform permutations
perm.dist <- matrix(0, length(len), NP)
for(i in 1:length(len)) {
  cat("Computing distribution for geneset size", i, "/", length(len), "\n")
  cat("\tElapsed time =", Sys.time()-t,"\n")
  perm.mat <- sapply(1:NP, function(j) sample(n, len[i]))
  perm.dist[i,] <- parApply(cl, perm.mat, 2, perm)
}

#### IF NEEDED: Get true p-values. Otherwise load from file
cat("\nCalculating true p-values\n\tElapsed time =", Sys.time()-t,"\n")
load("kegg.pathways.msigdb.synapses.RData") # kegg.pathways = table of pathways
paths <- unique(kegg.pathways[,2])
genelist <- sapply(paths, function(x)
  match(kegg.pathways[kegg.pathways[,2] %in% x,1], genes))
true.pvals <- parSapply(cl, genelist, perm)

stopCluster(cl)

### Calculate permutation p-values
cat("\nCalculating permutation p-values\n\tElapsed time =", Sys.time()-t,"\n")
len2 <- sapply(genelist, length)
map <- match(round(len2,-1), len)
perm.pvals <- sapply(1:length(map), function(x) sum(perm.dist[map[x],]<true.pvals[x])/NP)

### Save
save(true.pvals, perm.pvals, perm.dist, file="geneset.perm.results.RData")
cat("\nTotal elapsed time =", Sys.time()-t,"\n")
