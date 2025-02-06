## readin Caleb's data to do lasso modeling
caleb <- read.csv('~/Downloads/gene_summary_filtered.csv')
caleb$Matches

## some 'mismatches' due to gene annotation names being different
## can actually use these as examples of good probe sets
g  = 'ACTG2'
singlecell.summary.mm[g,]
xenium.summary.mm[g,]

par(mfrow=c(1,2), mar=rep(2,4)) 
cells.have <- intersect(rownames(emb.info), colnames(singlecell$counts))
MERINGUE::plotEmbedding(emb.info[cells.have,], col = log10(singlecell$counts[g, cells.have]+1), main=paste0('single cell:', g))
#points(emb.info[names(which(log10(singlecell$counts[g, cells.have]+1) > 0)),], col='red', pch=16)
cells.have2 <- intersect(rownames(emb.info), colnames(xenium$counts))
MERINGUE::plotEmbedding(emb.info[cells.have2,], col = log10(xenium$counts[g, cells.have2]+1), main=paste0('xenium:', g))
#points(emb.info[names(which(log10(xenium$counts[g, cells.have2]+1) > 0)),], col='red', pch=16)

## some genes have a lot of off target hits but some are not among the probe set of Visium
## sequencing should have these genes, though they may be non-expressed
offtargets <- c('ACTC1', 'POTEM', 'POTEJ', 'ACTB', 'POTEE', 'POTEI', 'POTEF', 'ACTA1', 'ACTG1', 'ACTBL2')

x = xenium.summary.mm[g,]
A = as.matrix(t(singlecell.summary.mm[c(g, offtargets),]))
x
A ## all genes in scRNA-seq

heatmap(cbind(xenium=xenium.summary.mm[g,], A), Rowv=NA, Colv=NA, scale='col')

library(glmnet)
lasso_model <- cv.glmnet(A, x, alpha=1, lower.limits=0)  # Enforce non-negative weights
best_lambda <- lasso_model$lambda.min
weights <- coef(lasso_model, s=best_lambda)[-1]  # Remove intercept
selected_columns <- which(weights > 0)
selected_weights <- weights[selected_columns]
print(selected_columns)  # Indices of chosen columns
print(selected_weights)  # Corresponding weights

colnames(A)[selected_columns]

par(mfrow=c(1,2), mar=rep(2,4)) 
cells.have <- intersect(rownames(emb.info), colnames(singlecell$counts))
new <- colSums(selected_weights * singlecell$counts[colnames(A)[selected_columns], cells.have])
MERINGUE::plotEmbedding(emb.info[cells.have,], col = log10(new+1), main=paste0('single cell: weighted composite'))

new <- colSums(singlecell$counts[colnames(A), cells.have])
MERINGUE::plotEmbedding(emb.info[cells.have,], col = log10(new+1), main=paste0('single cell: composite sum'))

