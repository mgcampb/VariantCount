#### Inputs
# libraries
require(parallel)
require(GSA)

# set up gene list
geneset_list <- GSA.read.gmt("test.gmt")
genelist <- geneset_list[[1]]
names(genelist) <- geneset_list[[2]]

load("Hg19HomRareAndNotHomRareCountTables.RData") 
# tab1 = hom rare counts (genes (rows) x samples (columns))
# tab2 = NOT hom rare counts (genes (rows) x samples (columns))
# dx = id to dx mapping (column 1 = id, column 2 = dx)

NP <- 10000 # number of permutations
cpus <- detectCores() # number of cores to use

#### Important variables
genes <- rownames(tab1)
n <- length(genes)
case <- dx[,2]=="CASE"