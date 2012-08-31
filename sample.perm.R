sample.fisher.with.permutation <- function(cl, genelist.rows, geneset_names, fisher.pvals, single.perm.fisher.test, t) {

  #### Perform sample.permutations
  M <- length(geneset_names)
  sample.perm.dist <- matrix(0, M, NP) # distribution of chi sq test p-values
  rownames(sample.perm.dist) <- geneset_names
  for(i in 1:M) {
    cat("Computing distribution for geneset", i, "/", M, "\n")
    cat("\tElapsed time =", Sys.time()-t,"\n")
    sample.perm.dist[i,] <- parSapply(cl, 1:NP, function(j) single.perm.fisher.test(genelist.rows[[i]]))
    #sample.perm.dist[i,] <- sapply(1:NP, function(j) single.perm.fisher.test(genelist.rows[[i]]))
  }
  
  ### Calculate sample.permutation p-values
  cat("\nCalculating sample.permutation p-values\n\tElapsed time =", Sys.time()-t,"\n")
  sample.perm.pvals <- sapply(1:M, function(i) sum(sample.perm.dist[i,] <= fisher.pvals[i])/NP)
  names(sample.perm.pvals) <- geneset_names
  return(sample.perm.pvals)

}