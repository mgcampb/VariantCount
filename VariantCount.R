t <- Sys.time()

#### Inputs
# libraries
require(parallel)
require(GSA)

# User-defined variable
NP <- 10000 # number of permutations
cpus <- detectCores() # number of cores to use
output_file <- "geneset_pvalues.csv"

# set up gene list
geneset_list <- GSA.read.gmt("data/todd.gmt")
genelist <- geneset_list[[1]]
names(genelist) <- geneset_list[[2]]
geneset_names <- names(genelist)

# get sample categories
dx <- read.csv("data/sample.category.csv")

# load variant count data
load("data/Hg19HomRareAndNotHomRareCountTables.RData") 
# tab.hom.counts = hom rare counts (genes (rows) x samples (columns))
# tab.non.hom.counts = NOT hom rare counts (genes (rows) x samples (columns))
tab.hom.counts <- tab1
tab.non.hom.counts <- tab2
rm(tab1, tab2)

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
single.perm.fisher.test <- function(sel,shuf=T) {
  if(shuf) case <- sample(case)
  case.ct1 <- sum(apply(tab.hom.counts[sel,case], 2, sum)) # hom alt case
  case.ct2 <- sum(apply(tab.non.hom.counts[sel,case], 2, sum)) # NOT hom alt case
  ctrl.ct1 <- sum(apply(tab.hom.counts[sel,!case], 2, sum)) # hom alt ctrl
  ctrl.ct2 <- sum(apply(tab.non.hom.counts[sel,!case], 2, sum)) # NOT hom alt ctrl
  m <- matrix(c(case.ct1,ctrl.ct1,case.ct2,ctrl.ct2),2,2)
  p <- fisher.test(m)$p.value
  if(is.na(p)) return(1) else return(p)
}

#### Parallelize
system("rm logfile")
cl <- makeCluster(cpus, outfile="logfile")

clusterExport(cl, c("tab.hom.counts", "tab.non.hom.counts", "genelist.rows", "case", "single.perm.fisher.test"))
clusterExport(cl,c("single.perm.fisher.test"))

#### Get chi-squared p-values for genesets
cat("\nCalculating chi-squared p-values\n\tElapsed time =", Sys.time()-t,"\n")
fisher.pvals <- parSapply(cl, genelist.rows, function(x) single.perm.fisher.test(x, shuf=F))
names(fisher.pvals) <- geneset_names

#### Compute sample permutation p-values
source("sample.perm.R")
sample.perm.pvals <- sample.fisher.with.permutation(cl, genelist.rows, geneset_names, fisher.pvals, single.perm.fisher.test)
#### Compute geneset permutation p-values
source("geneset.perm2.R")
geneset.perm.pvals <- geneset.fisher.with.permutation(cl, geneset_names, fisher.pvals, single.perm.fisher.test)

#stopCluster(cl)

#### Output table

output <- cbind(fisher.pvals, sample.perm.pvals, geneset.perm.pvals)
rownames(output) <- geneset_names
colnames(output) <- c("fisher_p", "sample_perm_p", "geneset_perm_p")
write.csv(output, file=output_file)

#### Visualization
source("make.heatmap4.R")
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
  fisher_p <- output[i,1]
  sample_perm_p <- output[i,2]
  geneset_perm_p <- output[i,3]
  make.heatmap(tab.hom.counts.case.geneset, tab.hom.counts.control.geneset, paste("figures/",geneset_names[i],".pdf",sep=""),fisher_p,sample_perm_p,geneset_perm_p)
})

cat("\nTotal elapsed time =", Sys.time()-t,"\n")

