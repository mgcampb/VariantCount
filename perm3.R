t <- Sys.time()

#### Inputs
# libraries
require(parallel)
require(GSA)

# set up gene list
geneset_list <- GSA.read.gmt("todd.gmt")
genelist <- geneset_list[[1]]
names(genelist) <- geneset_list[[2]]
geneset_names <- names(genelist)

load("data/Hg19HomRareAndNotHomRareCountTables.RData") 
# tab.hom.counts = hom rare counts (genes (rows) x samples (columns))
# tab.non.hom.counts = NOT hom rare counts (genes (rows) x samples (columns))
# dx = id to dx mapping (column 1 = id, column 2 = dx)
tab.hom.counts <- tab1
tab.non.hom.counts <- tab2
rm(tab1, tab2)

#test_genes = unlist(genelist)
#test_cases = dx[dx[,2]=='CASE',1][1:10]
#test_controls = dx[dx[,2]=='CONTROL',1][1:10]
#test_columns = match(c(as.character(test_cases), as.character(test_controls)), colnames(tab.hom.counts))
#tab.hom.counts = tab.hom.counts[test_genes,test_columns]
#tab.non.hom.counts = tab.non.hom.counts[test_genes,test_columns]
#dx = dx[test_columns,]

NP <- 10000 # number of permutations
cpus <- detectCores() # number of cores to use

#### Important variables
genes <- rownames(tab.hom.counts)
n <- length(genes)
case <- dx[,2]=="CASE"
genelist.rows <- sapply(genelist, function(x) {
  m <- match(x, genes)
  l <- length(m[is.na(m)])
  if(l>0) cat(l,"genes are not in count data:", m[is.na(m)])
  return(m[!is.na(m)])
}
)

#### Define chi-squared p-value function with optional permutation
single.perm.chisq.test <- function(sel,shuf=T) {
  if(shuf) case <- sample(case)
  case.ct1 <- sum(apply(tab.hom.counts[sel,case], 2, sum)) # hom alt case
  case.ct2 <- sum(apply(tab.non.hom.counts[sel,case], 2, sum)) # NOT hom alt case
  ctrl.ct1 <- sum(apply(tab.hom.counts[sel,!case], 2, sum)) # hom alt ctrl
  ctrl.ct2 <- sum(apply(tab.non.hom.counts[sel,!case], 2, sum)) # NOT hom alt ctrl
  p <- prop.test(matrix(c(case.ct1,ctrl.ct1,case.ct2,ctrl.ct2),2,2))$p.value
  if(is.na(p)) return(1) else return(p)
}

#### Parallelize
system("rm logfile")
cl <- makeCluster(cpus, outfile="logfile")

clusterExport(cl, c("tab.hom.counts", "tab.non.hom.counts", "genelist.rows", "case", "single.perm.chisq.test"))

#### Get chi-squared p-values for genesets
cat("\nCalculating chi-squared p-values\n\tElapsed time =", Sys.time()-t,"\n")
chisq.pvals <- parSapply(cl, genelist.rows, function(x) single.perm.chisq.test(x, shuf=F))
names(chisq.pvals) <- geneset_names

#### Compute sample permutation p-values
source("sample.perm.R")
sample.perm.pvals <- sample.chisq.with.permutation(cl, genelist.rows, geneset_names, chisq.pvals, single.perm.chisq.test)
#### Compute geneset permutation p-values
source("geneset.perm2.R")
geneset.perm.pvals <- geneset.chisq.with.permutation(cl, geneset_names, chisq.pvals, single.perm.chisq.test)

stopCluster(cl)

#### Visualization
source("make.heatmap3.R")
require(reshape2)
tab.hom.counts.melt <- melt(tab.hom.counts)
names(tab.hom.counts.melt) <- c("gene", "id", "hom.alt")

sapply(1:length(geneset_names), function(i) {
  tab.hom.counts.geneset <- tab.hom.counts.melt[tab.hom.counts.melt[,1] %in% genelist[[i]],]
  case.rows <- which(dx[match(tab.hom.counts.geneset[,2],dx[,1]),2] == "CASE")
  control.rows <- which(dx[match(tab.hom.counts.geneset[,2],dx[,1]),2] == "CONTROL")
  tab.hom.counts.case.geneset <- tab.hom.counts.geneset[case.rows,]
  tab.hom.counts.control.geneset <- tab.hom.counts.geneset[control.rows,]
  tab.hom.counts.geneset <- tab.hom.counts.melt[,]
  make.heatmap(tab.hom.counts.case.geneset, tab.hom.counts.control.geneset, paste("figures/",geneset_names[i],".pdf",sep=""))
})

cat("\nTotal elapsed time =", Sys.time()-t,"\n")

