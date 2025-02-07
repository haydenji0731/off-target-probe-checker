## use their existing celltype annotations:https://www.10xgenomics.com/products/xenium-in-situ/preview-dataset-human-breast
## make pseudobulks and show that each gene in xenium is a linear combination of genes in singlecell across celltypes

## goal: evaluate feasibility of deconvolving xenium gene expression using single cell rna-seq data
dir <- '~/OneDrive - Johns Hopkins/Data_Public/xenium_data/'

############# read in xenium data
library(rhdf5)
h5file <- paste0(dir, "Xenium_FFPE_Human_Breast_Cancer_Rep1_cell_feature_matrix.h5")
h5ls(h5file)

h5 <- h5read(h5file, "matrix")
barcodes <- as.character(h5read(h5file, "matrix/barcodes"))

library(Matrix)
counts <- sparseMatrix(
  dims = h5$shape,
  i = as.numeric(h5$indices),
  p = as.numeric(h5$indptr),
  x = as.numeric(h5$data),
  index1 = FALSE
)
colnames(counts) <- barcodes
rownames(counts) <- as.data.frame(h5[["features"]])$name
head(counts)

csvfile <- paste0(dir, 'Xenium_FFPE_Human_Breast_Cancer_Rep1_cells.csv.gz')
pos.info <- read.csv(csvfile)
head(pos.info)

pos <- as.matrix(cbind(x=pos.info$x_centroid, y=pos.info$y_centroid))
rownames(pos) <- pos.info$cell_id
xenium <- list(pos=pos,
               counts=counts)

################# read in single-cell data
dir2 <- paste0(dir, 'Visium_xenium_section/single_cell/FRP/')
h5file <- paste0(dir2, "Chromium_FFPE_Human_Breast_Cancer_Chromium_FFPE_Human_Breast_Cancer_count_sample_filtered_feature_bc_matrix.h5")
h5ls(h5file)

h5 <- h5read(h5file, "matrix")
barcodes <- as.character(h5read(h5file, "matrix/barcodes"))

library(Matrix)
counts <- sparseMatrix(
  dims = h5$shape,
  i = as.numeric(h5$indices),
  p = as.numeric(h5$indptr),
  x = as.numeric(h5$data),
  index1 = FALSE
)
colnames(counts) <- barcodes
rownames(counts) <- as.data.frame(h5[["features"]])$name
head(counts)

singlecell <- list(counts=counts)

#################

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

#cells.have <- intersect(names(celltype), rownames(emb.info)) ## from single_cell_exploration
cells.have <- intersect(names(celltype), colnames(singlecell$counts)) ## from single_cell_exploration

celltype <- celltype[cells.have]
celltype <- as.factor(celltype)
#par(mfrow=c(1,1))
#MERINGUE::plotEmbedding(emb.info[cells.have,], groups=celltype, main='single cell cell-types', mark.clusters = TRUE, mark.cluster.cex = 1)

mm <- model.matrix(~ 0 + celltype)
colnames(mm) <- levels(celltype)

singlecell.summary.mm <- singlecell$counts[, cells.have] %*% mm
head(singlecell.summary.mm)





celltype <- xenium_annot$Cluster
names(celltype) <- xenium_annot$Barcode

#cells.have <- intersect(names(celltype), rownames(emb.info)) ## from single_cell_exploration
cells.have <- intersect(names(celltype), colnames(xenium$counts))
celltype <- celltype[cells.have]
celltype <- as.factor(celltype)
#par(mfrow=c(1,1))
#MERINGUE::plotEmbedding(emb.info[cells.have,], groups=celltype, main='xenium cell-types', mark.clusters = TRUE, mark.cluster.cex = 1)

mm <- model.matrix(~ 0 + celltype)
colnames(mm) <- levels(celltype)

xenium.summary.mm <- xenium$counts[, cells.have] %*% mm
head(xenium.summary.mm)




## normalize after aggregating
#singlecell.summary.mm <- MERINGUE::normalizeCounts(singlecell.summary.mm, log=FALSE)
#xenium.summary.mm <- MERINGUE::normalizeCounts(xenium.summary.mm, log=FALSE)

## limit to same genes as Xenium
shared.genes <- intersect(rownames(singlecell$counts), rownames(xenium$counts))
shared.genes


m1 <- singlecell.summary.mm[shared.genes,]
m1 <- scale(m1)
m1 <- scale(t(m1))
m1[is.na(m1)] <- 0
heatmap(t(m1), scale='col')

m2 <- xenium.summary.mm[shared.genes,]
m2 <- scale(m2)
m2 <- scale(t(m2))
m2[is.na(m2)] <- 0
heatmap(t(m2), scale='none')


######### compare
heatmap(t(as.matrix(singlecell.summary.mm[shared.genes,])), Rowv=NA, Colv=NA)
heatmap(t(as.matrix(xenium.summary.mm[shared.genes,])), Rowv=NA, Colv=NA)

ct = 2
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


