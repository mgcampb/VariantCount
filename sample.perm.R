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
perm <- function(sel,shuf=T) {
  if(shuf) case <- sample(case)
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
clusterExport(cl, c("tab1", "tab2", "genelist", "case", "perm"))

#### Perform permutations
M <- length(paths)
perm.dist <- matrix(0, M, NP) # distribution of chi sq test p-values
rownames(perm.dist) <- paths
for(i in 1:M) {
  cat("Computing distribution for pathway", i, "/", M, "\n")
  cat("\tElapsed time =", Sys.time()-t,"\n")
  perm.dist[i,] <- parSapply(cl, 1:NP, function(i) perm(genelist[[i]]))
}

#### Get true p-values for KEGG pathways
cat("\nCalculating true p-values\n\tElapsed time =", Sys.time()-t,"\n")
true.pvals <- parSapply(cl, genelist, function(x) perm(x, shuf=F))
names(true.pvals) <- paths

stopCluster(cl)

### Calculate permutation p-values
cat("\nCalculating permutation p-values\n\tElapsed time =", Sys.time()-t,"\n")
perm.pvals <- sapply(1:M, function(i) sum(perm.dist[i,]<true.pvals[i])/NP)
names(perm.pvals) <- paths

### Save
save(true.pvals, perm.pvals, file="sample.perm.pvalues.RData")
save(perm.dist, file="sample.perm.dist.RData")
cat("\nTotal elapsed time =", Sys.time()-t,"\n")
