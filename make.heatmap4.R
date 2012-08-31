require(ggplot2)
require(grid)
require(gplots)

vplayout <- function(x, y)
viewport(layout.pos.row = x, layout.pos.col = y)

make.heatmap <- function(tab.hom.counts.case.geneset, tab.hom.counts.control.geneset, outfile, fisher_p,sample_perm_p,geneset_perm_p) {

  pdf(outfile, height=10, width=20)
  
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1, 2)))
  
  base_size <- 9 
  col1 <- "white"
  col2 <- "green"
  col3 <- "red"
  
  max_count = max(c(tab.hom.counts.case.geneset[,3],tab.hom.counts.control.geneset[,3]))
  tab.hom.counts.case.geneset[tab.hom.counts.case.geneset == 0] = -max_count
  tab.hom.counts.control.geneset[tab.hom.counts.control.geneset == 0] = -max_count
  p1 <- ggplot(tab.hom.counts.case.geneset, aes(id, gene)) + 
       geom_tile(aes(fill=hom.alt),colour=col1) +      
       scale_fill_gradient2(name="Homozygous\nalternate\ncounts", low=col1, mid=col2, high=col3, breaks=c(-max_count,1:max_count), labels=0:max_count) +
       theme_grey(base_size = base_size) + 
       labs(x = "", y = "") +
       scale_x_discrete(expand = c(0, 0)) + 
       scale_y_discrete(expand = c(0, 0)) +     
       opts(title = sprintf("CASE\n fisher_p: %.3f, sample_perm_p: %.3f, geneset_perm_p: %.3f",fisher_p,sample_perm_p,geneset_perm_p),
            axis.ticks = theme_blank(), axis.text.x = theme_blank())  
  
  p2 <- ggplot(tab.hom.counts.control.geneset, aes(id, gene)) + 
       geom_tile(aes(fill=hom.alt),colour=col1) +      
       scale_fill_gradient2(name = "Homozygous\nalternate\ncounts", low=col1, mid=col2, high=col3, breaks=c(-max_count,1:max_count), labels=0:max_count) + 
       theme_grey(base_size = base_size) + 
       labs(x = "", y = "") + 
       scale_x_discrete(expand = c(0, 0)) + 
       scale_y_discrete(expand = c(0, 0)) +     
       opts(title = sprintf("CONTROL\n fisher_p: %.3f, sample_perm_p: %.3f, geneset_perm_p: %.3f",fisher_p,sample_perm_p,geneset_perm_p),
            axis.ticks = theme_blank(), axis.text.x = theme_blank())  
  
  print(p1, vp=vplayout(1,1))
  print(p2, vp=vplayout(1,2))
  
  dev.off()
  
}