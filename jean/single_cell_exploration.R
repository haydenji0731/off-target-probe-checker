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

################## compare total expression
shared.genes <- intersect(rownames(xenium$counts), rownames(singlecell$counts))
length(shared.genes)

df <- data.frame(
  xenium_total = rowSums(xenium$counts[shared.genes,]),
  singlecell_total = rowSums(singlecell$counts[shared.genes,]),
  label=shared.genes
)
library(ggplot2)
ggplot(df, aes(y=xenium_total, x = singlecell_total, label=label)) + geom_point() + scale_x_log10() + scale_y_log10() +
  ggrepel::geom_label_repel()


################### do harmonized clustering
## https://portals.broadinstitute.org/harmony/mudan.html
cd <- cbind(xenium$counts[shared.genes,], singlecell$counts[shared.genes,])
meta <- c(rep('xenium', ncol(xenium$counts)), rep('singlecell', ncol(singlecell$counts)))
names(meta) <- c(colnames(xenium$counts), colnames(singlecell$counts))
meta <- factor(meta)
table(meta)

## filter out single cells without any counts
vi <- colSums(cd) > 10
cd <- cd[, vi]
meta <- meta[vi]
table(meta)

## normal clustering to sanity check
## CPM normalization
mat <- MERINGUE::normalizeCounts(cd, 
                              verbose=FALSE) 

## 30 PCs on overdispersed genes
pcs <- irlba::prcomp_irlba(t(mat), n=30)
names(pcs)
dim(pcs$x)

pc <- pcs$x
rownames(pc) <- colnames(cd)
MERINGUE::plotEmbedding(pc[,1:2], groups=meta, show.legend = TRUE)

## harmonize
library(harmony)
harmonized <- HarmonyMatrix(pc, meta, do_pca = FALSE, verbose = FALSE, theta=4)
dim(harmonized)

# UMAP embedding with regular PCs
emb.info <- uwot::umap(harmonized)
rownames(emb.info) <- rownames(pc)
head(emb.info)

# Plot
par(mfrow=c(1,1), mar=rep(2,4))
MERINGUE::plotEmbedding(emb.info, groups=meta, 
                     show.legend=TRUE, xlab=NA, ylab=NA, 
                     main='Harmonized tSNE Embedding',
                     verbose=FALSE)

par(mfrow=c(1,2), mar=rep(2,4)) 
g <- 'ACTG2'
cells.have <- intersect(rownames(emb.info), colnames(singlecell$counts))
MERINGUE::plotEmbedding(emb.info[cells.have,], col = log10(singlecell$counts[g, cells.have]+1), main=paste0('single cell:', g))
cells.have2 <- intersect(rownames(emb.info), colnames(xenium$counts))
MERINGUE::plotEmbedding(emb.info[cells.have2,], col = log10(xenium$counts[g, cells.have2]+1), main=paste0('xenium:', g))



####### summarize
cells.have <- intersect(rownames(emb.info), colnames(singlecell$counts))
mm <- model.matrix(~ 0 + com[cells.have])
colnames(mm) <- levels(com)

singlecell.summary.mm <- singlecell$counts[, cells.have] %*% mm
head(singlecell.summary.mm)

cells.have2 <- intersect(rownames(emb.info), colnames(xenium$counts))
mm <- model.matrix(~ 0 + com[cells.have2])
colnames(mm) <- levels(com[cells.have2])

xenium.summary.mm <- xenium$counts[, cells.have2] %*% mm
head(xenium.summary.mm)



m1 <- singlecell.summary.mm[shared.genes,]
m1 <- scale(m1)
m1 <- scale(t(m1))
m1[is.na(m1)] <- 0
heatmap(t(m1), scale='none')

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



