#create test dataset
#test_genes = unlist(genelist)
#test_cases = dx[dx[,2]=='CASE',1][1:10]
#test_controls = dx[dx[,2]=='CONTROL',1][1:10]
#test_columns = match(c(as.character(test_cases), as.character(test_controls)), colnames(tab.hom.counts))
#tab.hom.counts = tab.hom.counts[test_genes,test_columns]
#tab.non.hom.counts = tab.non.hom.counts[test_genes,test_columns]
#dx = dx[test_columns,]