t <- Sys.time()

files <- list.files(pattern="Kegg.*txt") 
pvals <- numeric(length(files))
names(pvals) <- gsub("RareVariantSummaryCounts.txt","",files)
counts <- matrix(0, length(files), 6)
rownames(counts) <- names(pvals)
colnames(counts) <- c("HomRefCase", "HetCase", "HomAltCase", "HomRefCtrl", "HetCtrl", "HomAltCtrl")
counts2 <- list()

for(i in 1:length(files)) {
  cat(i, "/", length(files), ":\n\tLoading", files[i], "...\nElapsed time = ",
    Sys.time()-t,"\n")
  tab <- read.delim(files[i])
  case <- tab[,3] == "CASE"
  case.counts <- apply(tab[case,4:6], 2, sum)
  ctrl.counts <- apply(tab[!case,4:6], 2, sum)
  ct = rbind(rev(c(sum(case.counts[1],case.counts[2]),case.counts[3])),
    rev(c(sum(ctrl.counts[1],ctrl.counts[2]),ctrl.counts[3])))
  colnames(ct) = c("HomAlt","Not-HomAlt")
  rownames(ct) = c("Case","Ctrl")
  pvals[i] <- prop.test(ct)$p.value
  counts[i,] <- c(case.counts, ctrl.counts)
  counts2[[i]] <- ct
}

names(counts2) <- names(pvals)

x <- apply(counts, 2, mean)
means <- rbind(c(sum(x[1],x[2]),x[3]),c(sum(x[4],x[5]),x[6]))
means[1,] <- rev(means[1,])
means[2,] <- rev(means[2,])
colnames(means) <- c("HomAlt","Not-HomAlt")
rownames(means) <- c("Case","Ctrl")
prop.test(means)

x <- apply(counts, 2, sum)
sums <- rbind(c(sum(x[1],x[2]),x[3]),c(sum(x[4],x[5]),x[6]))
sums[1,] <- rev(sums[1,])
sums[2,] <- rev(sums[2,])
colnames(sums) <- c("HomAlt","Not-HomAlt")
rownames(sums) <- c("Case","Ctrl")
prop.test(sums)

cat("Total elapsed time:", Sys.time()-t,"\n")
