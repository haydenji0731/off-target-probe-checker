## do a negative control test
## if a gene doesn't have off targets, does our approach reject off targets as being likely based on expression

g <- 'MS4A1'
#offtargets <- c('ACTC1', 'POTEM', 'POTEJ', 'ACTB', 'POTEE', 'POTEI', 'POTEF', 'ACTA1', 'ACTG1', 'ACTBL2') ## just a test
offtargets <- c('CD19', 'HLA-A') ## biological correlates as control of negative control (doesn't work well surprisingly)

x = log10(xenium.summary.mm[g,]+1)
y = log10(singlecell.summary.mm[g,]+1)
A = log10(as.matrix(t(singlecell.summary.mm[offtargets,]))+1)
x
y
A ## all genes in scRNA-seq

## our general question is, solve for B x = y + B*A
par(mfrow=c(1,2))
plot(x,y)
plot(x,y + rowSums(A))

cor(x,y)
cor(x, y + rowSums(A))

## conclusion: can validate, can't discover off targets by this approach