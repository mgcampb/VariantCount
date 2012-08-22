t <- Sys.time()

#### Inputs
require(parallel)
load("Hg19HomRareAndNotHomRareCountTables.RData") 
  # tab1 = hom rare counts (genes (rows) x samples (columns))
  # tab2 = NOT hom rare counts (genes (rows) x samples (columns))
  # dx = id to dx mapping (column 1 = id, column 2 = dx)
load("kegg.pathways.msigdb.synapses.RData") 
  # kegg.pathways = table of pathways
NP <- 10000 # number of permutations
cpus <- detectCores() # number of cores to use

#### Important variables
genes <- rownames(tab1)
case <- dx[,2]=="CASE"
n <- length(genes)
paths <- unique(kegg.pathways[,2])
genelist <- sapply(paths, function(x) {
    m <- match(kegg.pathways[kegg.pathways[,2] %in% x,1], genes)
    m[!is.na(m)]
  }
)

#### Define permutation function
perm <- function(sel) {
  case.ct1 <- sum(apply(tab1[sel,case], 2, sum)) # hom alt case
  case.ct2 <- sum(apply(tab2[sel,case], 2, sum)) # NOT hom alt case
  ctrl.ct1 <- sum(apply(tab1[sel,!case], 2, sum)) # hom alt ctrl
  ctrl.ct2 <- sum(apply(tab2[sel,!case], 2, sum)) # NOT hom alt ctrl
  p <- prop.test(matrix(c(case.ct1,ctrl.ct1,case.ct2,ctrl.ct2),2,2))$p.value
  if(is.na(p)) return(1) else return(p)
}

#### Parallelize
system("rm logfile")
cl <- makeCluster(cpus, outfile="logfile")
clusterExport(cl, c("tab1", "tab2", "case", "perm"))

#### Perform permutations
len <- sort(unique(sapply(genelist,length)))
perm.dist <- matrix(0, length(len), NP) # distribution of chi sq test p-values
rownames(perm.dist) <- as.character(len)
for(i in 1:length(len)) {
  cat("Computing distribution for geneset size", i, "/", length(len), "\n")
  cat("\tElapsed time =", Sys.time()-t,"\n")
  perm.mat <- sapply(1:NP, function(j) sample(n, len[i]))
  perm.dist[i,] <- parApply(cl, perm.mat, 2, perm)
}

#### Get true p-values for KEGG pathways
cat("\nCalculating true p-values\n\tElapsed time =", Sys.time()-t,"\n")
true.pvals <- parSapply(cl, genelist, perm)
names(true.pvals) <- paths

stopCluster(cl)

### Calculate permutation p-values
cat("\nCalculating permutation p-values\n\tElapsed time =", Sys.time()-t,"\n")
map <- match(sapply(genelist, length), len)
perm.pvals <- sapply(1:length(map), function(x) sum(perm.dist[map[x],]<true.pvals[x])/NP)
names(perm.pvals) <- paths

### Save
save(true.pvals, perm.pvals, file="geneset.perm.pvalues.RData")
save(perm.dist, file="geneset.perm.dist.RData")
cat("\nTotal elapsed time =", Sys.time()-t,"\n")
