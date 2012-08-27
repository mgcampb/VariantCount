geneset.chisq.with.permutation <- function(cl, geneset_names, chisq.pvals, single.perm.chisq.test) {
  
  #### Perform permutations
  len <- sort(unique(sapply(genelist,length)))
  geneset.perm.dist <- matrix(0, length(len), NP) # distribution of chi sq test p-values
  rownames(geneset.perm.dist) <- as.character(len)
  for(i in 1:length(len)) {
    cat("Computing distribution for geneset size", i, "/", length(len), "\n")
    cat("\tElapsed time =", Sys.time()-t,"\n")
    geneset.perm.mat <- sapply(1:NP, function(j) sample(n, len[i]))
    geneset.perm.dist[i,] <- parApply(cl, geneset.perm.mat, 2, single.perm.chisq.test)
  }
  
  ### Calculate permutation p-values
  cat("\nCalculating permutation p-values\n\tElapsed time =", Sys.time()-t,"\n")
  map <- match(sapply(genelist, length), len)
  geneset.perm.pvals <- sapply(1:length(map), function(x) sum(geneset.perm.dist[map[x],]<chisq.pvals[x])/NP)
  names(geneset.perm.pvals) <- geneset_names
  return(geneset.perm.pvals)
  
}