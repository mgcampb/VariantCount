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


system("rm logfile")
cl <- makeCluster(cpus, outfile="logfile")

case <- tab[,3] == "CASE"
case.ct <- apply(tab[case,4:6], 2, sum)
ctrl.ct <- apply(tab[!case,4:6], 2, sum)
case.ct <- c(case.ct[3], sum(case.ct)-case.ct[3])
ctrl.ct <- c(ctrl.ct[3], sum(ctrl.ct)-ctrl.ct[3])
ct <- rbind(case.ct, ctrl.ct)
pval <- prop.test(ct)$p.value





perm.distr <- parSapply(cl, 1:NP, perm)
stopCluster(cl)

save(p.case, p.ctrl, perm.distr, file="OUTPUT.RData")

cat("Elapsed time:", Sys.time()-t,"\n")
