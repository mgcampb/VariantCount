tab <- read.delim("KeggRibosomeRareVariantSummaryCounts.txt")
genes = unique(tab[,1])
ids = unique(tab[,2])
mat = matrix(0, length(genes), length(ids))
dx = tab[match(ids,tab[,2]),3]
ids2 = ids[order(dx)]
colnames(mat) = ids2
rownames(mat) = genes
for(j in 1:length(ids2)) {
    mat[,j] = tab[tab[,2]==ids2[j],6]
}

require(gplots)
col = c("gray90",heat.colors(max(mat))) # ,"orange","purple","cyan","green","brown")[1:(max(mat)+1)]
heatmap(mat, Rowv=NA, Colv=NA, scale="none",
        labCol=NA, col=col,
        ColSideColors=rep(c("red","blue"),times=table(dx)),margins=c(0,0))
legend(x=0,y=.7,as.character(0:max(mat)),
       fill=col,cex=0.6,bty="n")


