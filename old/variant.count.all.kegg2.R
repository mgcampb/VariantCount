t <- Sys.time()

files <- list.files(pattern="Kegg.*txt") 
pvals <- numeric(length(files))
names(pvals) <- gsub("RareVariantSummaryCounts.txt","",files)
counts <- list()

for(i in 1:length(files)) {
  cat("Loading pathway",i,"/",length(files),"\n\tElapsed time =",Sys.time()-t,"\n")
  tab <- read.delim(files[i])
  case <- tab[,3] == "CASE"
  case.ct <- apply(tab[case,4:6], 2, sum)
  ctrl.ct <- apply(tab[!case,4:6], 2, sum)
  case.ct <- c(case.ct[3], sum(case.ct)-case.ct[3])
  ctrl.ct <- c(ctrl.ct[3], sum(ctrl.ct)-ctrl.ct[3])
  ct <- rbind(case.ct, ctrl.ct)
  colnames(ct) <- c("HomAlt","Not-HomAlt")
  rownames(ct) <- c("Case","Ctrl")
  pvals[i] <- prop.test(ct)$p.value
  counts[[i]] <- ct
}

names(counts) <- names(pvals)

save(pvals, counts, file="kegg.msigdb.results.7.12.12.RData")

cat("Total elapsed time:", Sys.time()-t,"\n")
