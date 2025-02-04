## use their existing celltype annotations:https://www.10xgenomics.com/products/xenium-in-situ/preview-dataset-human-breast
## make pseudobulks and show that each gene in xenium is a linear combination of genes in singlecell across celltypes
library(readxl)

xenium_annot <- read_excel('~/Downloads/Cell_Barcode_Type_Matrices.xlsx', sheet=4)
head(xenium_annot)

singlecell_annot <- read_excel('~/Downloads/Cell_Barcode_Type_Matrices.xlsx', sheet = 1)
head(singlecell_annot)

table(singlecell_annot$Barcode %in% colnames(singlecell$counts))
table(xenium_annot$Barcode %in% colnames(xenium$counts))



## make pseudobulks https://jef.works/blog/2020/04/06/quickly-creating-pseudobulks/
celltype <- singlecell_annot$Annotation
names(celltype) <- singlecell_annot$Barcode

cells.have <- intersect(names(celltype), rownames(emb.info)) ## from single_cell_exploration
celltype <- celltype[cells.have]
celltype <- as.factor(celltype)
par(mfrow=c(1,1))
MERINGUE::plotEmbedding(emb.info[cells.have,], groups=celltype, main='single cell cell-types', mark.clusters = TRUE, mark.cluster.cex = 1)

mm <- model.matrix(~ 0 + celltype)
colnames(mm) <- levels(celltype)

singlecell.summary.mm <- singlecell$counts[, cells.have] %*% mm
head(singlecell.summary.mm)





celltype <- xenium_annot$Cluster
names(celltype) <- xenium_annot$Barcode

cells.have <- intersect(names(celltype), rownames(emb.info)) ## from single_cell_exploration
celltype <- celltype[cells.have]
celltype <- as.factor(celltype)
par(mfrow=c(1,1))
MERINGUE::plotEmbedding(emb.info[cells.have,], groups=celltype, main='xenium cell-types', mark.clusters = TRUE, mark.cluster.cex = 1)

mm <- model.matrix(~ 0 + celltype)
colnames(mm) <- levels(celltype)

xenium.summary.mm <- xenium$counts[, cells.have] %*% mm
head(xenium.summary.mm)




## normalize after aggregating
singlecell.summary.mm <- MERINGUE::normalizeCounts(singlecell.summary.mm, log=FALSE)
xenium.summary.mm <- MERINGUE::normalizeCounts(xenium.summary.mm, log=FALSE)

######### compare
heatmap(t(as.matrix(singlecell.summary.mm[shared.genes,])), Rowv=NA, Colv=NA)
heatmap(t(as.matrix(xenium.summary.mm[shared.genes,])), Rowv=NA, Colv=NA)

ct = 13
colnames(singlecell.summary.mm)[ct]
colnames(xenium.summary.mm)[ct]
df <- data.frame(singlecell = singlecell.summary.mm[shared.genes,1],
                xenium = xenium.summary.mm[shared.genes,1],
                label = shared.genes)
ggplot(df, aes(x=singlecell, y=xenium, label=label)) + geom_point() + ggrepel::geom_label_repel() +
  scale_x_log10() + scale_y_log10() + ggtitle(colnames(singlecell.summary.mm)[ct])


########### focus on APOBEC3B
g <- 'APOBEC3B'
singlecell.summary.mm[g,]
xenium.summary.mm[g,]

## look across singlecell data...which gene is most correlated
results <- do.call(rbind, lapply(1:nrow(singlecell.summary.mm), function(i) {
  cpv <- cor.test(x=singlecell.summary.mm[i,], y=xenium.summary.mm[g,1:ncol(singlecell.summary.mm)])
  c(estimate=cpv$estimate, pv=cpv$p.value)
}))
rownames(results) <- rownames(singlecell.summary.mm)
head(results)

results['APOBEC3F',]
head(results[order(results[,1], decreasing=TRUE),]) ## surprisingly not highest correlated...so likely linear combo not just one gene...hummm

########### focus on CEACAM8
g <- 'CEACAM8'
singlecell.summary.mm[g,]
xenium.summary.mm[g,]

## look across singlecell data...which gene is most correlated
results <- do.call(rbind, lapply(1:nrow(singlecell.summary.mm), function(i) {
  cpv <- cor.test(x=singlecell.summary.mm[i,], y=xenium.summary.mm[g,1:ncol(singlecell.summary.mm)])
  c(estimate=cpv$estimate, pv=cpv$p.value)
}))
rownames(results) <- rownames(singlecell.summary.mm)
head(results)

results['CEACAM5',]
head(results[order(results[,1], decreasing=TRUE),])


