load("/Users/malcolmcampbell/Desktop/Work/data/variant.count/KeggGlutamatergicRareVariantSummaryCounts.RData")

# "tab" is a table of rare variant counts for cases/controls
# in the kegg glutamatergic synapse.
#
# nrows = ngenes (125) x nsamples (790)
#
# For a variant to be included in this table, it must have 
# an alternate allele that is rare, and this allele must be
# present in at least one sample.
#
# Moreover every alternate allele in this table is rare (AF<1%).
#
# This script calculates the proportion of all genotypes that are 
# homozygous alternate, for cases and controls separately.
#
# It bootstraps to find a distribution for this number.
 
# calculate proportion of case/control variants that are 
# homozygous alternate
case = tab[,3] == "CASE"
cases = unique(tab[case,2])
ctrls = unique(tab[!case,2])
case.counts = apply(tab[case,4:6], 2, sum)
ctrl.counts = apply(tab[!case,4:6], 2, sum)
counts = rbind(rev(c(sum(case.counts[1],case.counts[2]),case.counts[3])),
  rev(c(sum(ctrl.counts[1],ctrl.counts[2]),ctrl.counts[3])))
colnames(counts) = c("HomoAlt","Not-HomoAlt")
rownames(counts) = c("Case","Ctrl")
p.case = counts[1,1]/sum(counts[1,])
p.ctrl = counts[2,1]/sum(counts[2,])


# define bootstrap function that samples with replacement
bootstrap <- function(i) {
  cases2 = sample(cases, length(cases), replace=T)
  ctrls2 = sample(ctrls, length(ctrls), replace=T)
  tab2 = tab[tab[,2] %in% union(cases2, ctrls2),]
  case2 = tab2[,3] == "CASE"
  case.counts2 = apply(tab2[case2,4:6], 2, sum)
  ctrl.counts2 = apply(tab2[!case2,4:6], 2, sum)
  p.case2 = case.counts2[3]/sum(case.counts2)
  p.ctrl2 = ctrl.counts2[3]/sum(ctrl.counts2)
  return(c(p.case2,p.ctrl2))
}

# distribution of proportion of homozygous rare variants out of all rare variants
# 1st col = case, 2nd col = ctrl
distr = t(sapply(1:1000, bootstrap))
colnames(distr) = paste(colnames(distr),c("Case","Control"),sep="")

# function to plot error bars
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

distr.pct <- 100*distr
p.bar <- apply(distr.pct,2,mean)
p.sd <- apply(distr.pct,2,sd)
pval <- 1-pbinom(p.bar[1]*sum(case.counts)/100-1,sum(case.counts),p.bar[2]/100)

title <- paste("Homozygous rare variants are more frequent in cases
in the KEGG Glutamatergic Synapse (p = ",signif(pval,2),")",sep="")
barx <- barplot(p.bar,col=c("red","blue"), ylim = c(0, max(p.bar)+max(p.sd)),
                ylab="% rare variants that are homozygous alternate", names = c("Case", "Control"), 
                main=title)
error.bar(barx,p.bar,p.sd)

# permutation p-value
perm <- function(i) {
  case2 = sample(rep(c(T,F),times=table(tab[,3])))
  case.counts2 = apply(tab[case2,4:6], 2, sum)
  ctrl.counts2 = apply(tab[!case2,4:6], 2, sum)
  p.case2 = case.counts2[3]/sum(case.counts2)
  p.ctrl2 = ctrl.counts2[3]/sum(ctrl.counts2)
  return(p.case2/p.ctrl2)
}

perm.distr = sapply(1:1000, perm)


