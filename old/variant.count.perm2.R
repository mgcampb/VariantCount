t <- Sys.time()

require(parallel)
cpus <- detectCores()
NP <- 1000

files <- list.files(pattern="Kegg.*txt") 
perm.pvals <- numeric(length(files))
names(perm.pvals) <- gsub("RareVariantSummaryCounts.txt","",files)

# this function performs one permutation
perm <- function(i) {
  case2 <- sample(rep(c(T,F),times=table(tab[,3])))
  case.ct2 <- apply(tab[case2,4:6], 2, sum)
  ctrl.ct2 <- apply(tab[!case2,4:6], 2, sum)
  case.ct2 <- c(case.ct2[3], sum(case.ct2)-case.ct2[3])
  ctrl.ct2 <- c(ctrl.ct2[3], sum(ctrl.ct2)-ctrl.ct2[3])
  prop.test(rbind(case.ct2, ctrl.ct2))$p.value
}

# start cluster
system("rm logfile")
cl <- makeCluster(cpus, outfile="logfile")
clusterExport(cl, c("perm", "t"))

for(i in 1:length(files)) {

  cat("Loading pathway",i,"/",length(files),":",names(perm.pvals)[i],
    "\n\tElapsed time =",Sys.time()-t,"\n")
  tab <- read.delim(files[i])
  clusterExport(cl, "tab")
  
  case <- tab[,3] == "CASE"
  case.ct <- apply(tab[case,4:6], 2, sum)
  ctrl.ct <- apply(tab[!case,4:6], 2, sum)
  case.ct <- c(case.ct[3], sum(case.ct)-case.ct[3])
  ctrl.ct <- c(ctrl.ct[3], sum(ctrl.ct)-ctrl.ct[3])
  p <- prop.test(rbind(case.ct, ctrl.ct))$p.value
  
  p.distr <- parSapply(cl, 1:NP, perm)
  perm.pvals[i] <- sum(p.distr < p) / NP
  
}

stopCluster(cl)

save(perm.pvals, file="perm.pvals.RData")

cat("Elapsed time:", Sys.time()-t,"\n")
