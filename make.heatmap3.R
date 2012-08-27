require(ggplot2)
require(grid)
require(gplots)

vplayout <- function(x, y)
viewport(layout.pos.row = x, layout.pos.col = y)

make.heatmap <- function(tab.hom.counts.case.geneset, tab.hom.counts.control.geneset, outfile) {

  pdf(outfile, height=10, width=20)
  
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1, 2)))
  
  base_size <- 9 
  col1 <- "steelblue"
  col2 <- "steelblue"
  col3 <- "white"
  p1 <- ggplot(tab.hom.counts.case.geneset, aes(id, gene)) + 
       geom_tile(aes(fill=hom.alt),colour=col1) +      
       scale_fill_gradient(name="Homozygous\nalternate\ncounts", low=col2, high=col3, breaks=0:max(tab.hom.counts.case.geneset[,3])) +
       theme_grey(base_size = base_size) + 
       labs(x = "", y = "") + 
       scale_x_discrete(expand = c(0, 0)) + 
       scale_y_discrete(expand = c(0, 0)) +     
       opts(title = "CASE", axis.ticks = theme_blank(), axis.text.x = theme_blank())  
  
  p2 <- ggplot(tab.hom.counts.control.geneset, aes(id, gene)) + 
       geom_tile(aes(fill=hom.alt),colour=col1) +      
       scale_fill_gradient(name = "Homozygous\nalternate\ncounts", low=col2, high=col3, breaks=0:max(tab.hom.counts.control.geneset[,3])) + 
       theme_grey(base_size = base_size) + 
       labs(x = "", y = "") + 
       scale_x_discrete(expand = c(0, 0)) + 
       scale_y_discrete(expand = c(0, 0)) +     
       opts(title = "CONTROL", axis.ticks = theme_blank(), axis.text.x = theme_blank())
  
  print(p1, vp=vplayout(1,1))
  print(p2, vp=vplayout(1,2))
  
  dev.off()
  
}