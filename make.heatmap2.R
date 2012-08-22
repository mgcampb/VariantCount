setwd("~/Desktop/Work/data/variant.count/")
files <- list.files("AllKeggPathways", "*.txt", full.names=T)

require(ggplot2)
require(grid)
require(gplots)

vplayout <- function(x, y)
viewport(layout.pos.row = x, layout.pos.col = y)

for(file in files) {
  tab <- read.delim(file)
  
  names(tab) <- c("gene", "id", "pheno", "hom.ref", "het.alt", "hom.alt")
  tab1 = tab[tab[,3]=="CASE",]
  tab2 = tab[tab[,3]=="CONTROL",]

  outfile = gsub("RareVariantSummaryCounts.txt","",gsub("AllKeggPathways/","",file))
  outfile = paste("~/Desktop/RareVariantKeggPathwayHeatmaps/",outfile,".pdf",sep="")
  pdf(outfile, height=10, width=20)
  
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1, 2)))
  
  base_size <- 9 
  col1 <- "steelblue"
  col2 <- "steelblue"
  col3 <- "white"
  p1 <- ggplot(tab1, aes(id, gene)) + 
       geom_tile(aes(fill=hom.alt),colour=col1) +      
       scale_fill_gradient(name="Homozygous\nalternate\ncounts", low=col2, high=col3, breaks=0:max(tab1[,6])) +
       theme_grey(base_size = base_size) + 
       labs(x = "", y = "") + 
       scale_x_discrete(expand = c(0, 0)) + 
       scale_y_discrete(expand = c(0, 0)) +     
       opts(title = "CASE", axis.ticks = theme_blank(), axis.text.x = theme_blank())  
  
  p2 <- ggplot(tab2, aes(id, gene)) + 
       geom_tile(aes(fill=hom.alt),colour=col1) +      
       scale_fill_gradient(name = "Homozygous\nalternate\ncounts", low=col2, high=col3, breaks=0:max(tab2[,6])) + 
       theme_grey(base_size = base_size) + 
       labs(x = "", y = "") + 
       scale_x_discrete(expand = c(0, 0)) + 
       scale_y_discrete(expand = c(0, 0)) +     
       opts(title = "CONTROL", axis.ticks = theme_blank(), axis.text.x = theme_blank())
  
  print(p1, vp=vplayout(1,1))
  print(p2, vp=vplayout(1,2))
  
  dev.off()
  
}