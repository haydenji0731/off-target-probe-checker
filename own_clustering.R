## do our own clustering
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
dir2 <- paste0(dir, 'Visium_xenium_section/single_cell/3prime/')
h5file <- paste0(dir2, "SC3pv3_GEX_Breast_Cancer_DTC_Aggr_count_filtered_feature_bc_matrix.h5")
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
harmonized <- HarmonyMatrix(pc, meta, do_pca = FALSE, verbose = FALSE, theta=8)
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
points(emb.info[names(which(log10(singlecell$counts[g, cells.have]+1) > 0)),], col='red', pch=16)

cells.have2 <- intersect(rownames(emb.info), colnames(xenium$counts))
MERINGUE::plotEmbedding(emb.info[cells.have2,], col = log10(xenium$counts[g, cells.have2]+1), main=paste0('xenium:', g))
#points(emb.info[names(which(log10(xenium$counts[g, cells.have2]+1) > 0)),], col='red', pch=16)


##################### do our own clustering
getApproxComMembership <- function (mat, k, nsubsample = nrow(mat) * 0.5, method = igraph::cluster_walktrap,
                                    seed = 0, vote = FALSE, verbose = TRUE)
{
  if (verbose) {
    print(paste0("Subsampling from ", nrow(mat), " cells to ",
                 nsubsample, " ... "))
  }
  set.seed(seed)
  subsample <- sample(rownames(mat), nsubsample)
  if (verbose) {
    print("Identifying cluster membership for subsample ... ")
  }
  pcs.sub <- mat[subsample, ]
  com.sub <- MERINGUE::getClusters(pcs.sub, k = k, method = method)
  if (verbose) {
    print("Imputing cluster membership for rest of cells ... ")
  }
  if (vote) {
    data <- mat[subsample, ]
    query <- mat[setdiff(rownames(mat), subsample), ]
    knn <- RANN::nn2(data, query, k = k)[[1]]
    rownames(knn) <- rownames(query)
    com.nonsub <- unlist(apply(knn, 1, function(x) {
      nn <- rownames(data)[x]
      nn.com <- com.sub[nn]
      return(names(sort(table(nn.com), decreasing = TRUE)[1]))
    }))
    com.all <- factor(c(com.sub, com.nonsub)[rownames(mat)])
  }
  else {
    df.sub <- data.frame(celltype = com.sub, pcs.sub)
    model <- MASS::lda(celltype ~ ., data = df.sub)
    df.all <- data.frame(mat)
    model.output <- stats::predict(model, df.all)
    com.all <- model.output$class
    names(com.all) <- rownames(df.all)
    if (verbose) {
      print("Model accuracy for subsample ...")
      print(table(com.all[names(com.sub)] == com.sub))
    }
  }
  return(com.all)
}
com <- getApproxComMembership(harmonized, k=10, method=igraph::cluster_louvain)

par(mfrow=c(1,2), mar=rep(2,4)) 
cells.have <- intersect(rownames(emb.info), colnames(singlecell$counts))
MERINGUE::plotEmbedding(emb.info[cells.have,], groups=com[cells.have], main='single cell clusters')
cells.have2 <- intersect(rownames(emb.info), colnames(xenium$counts))
MERINGUE::plotEmbedding(emb.info[cells.have2,], groups=com[cells.have2], main='xenium clusters')


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
#m1 <- scale(m1)
m1 <- scale(t(m1))
m1[is.na(m1)] <- 0
heatmap(t(m1), scale='none')

m2 <- xenium.summary.mm[shared.genes,]
#m2 <- scale(m2)
m2 <- scale(t(m2))
m2[is.na(m2)] <- 0
heatmap(t(m2), scale='none')

## use same clustering
hc <- hclust(dist(t(m1)))
plot(hc)
rc <- hclust(dist(m1))
plot(rc)  

heatmap(m1, Rowv=as.dendrogram(rc), Colv=as.dendrogram(hc), scale='none')
heatmap(m2, Rowv=as.dendrogram(rc), Colv=as.dendrogram(hc), scale='none')

heatmap(m1-m2, Rowv=as.dendrogram(rc), Colv=as.dendrogram(hc), scale='row')

par(mfrow=c(1,1))
barplot(sort(colSums(m1-m2)))


########### focus on a gene
g <- 'ACTG2'
singlecell.summary.mm[g,]
xenium.summary.mm[g,]

## until calebs predictions come...use regex orthologs
offtargets <- rownames(singlecell.summary.mm)[grepl('^ACT', rownames(singlecell.summary.mm))]
offtargets

## potentially off targets
#singlecell.summary.mm['ACTG1',]
#singlecell.summary.mm['ACTA1',]
#singlecell.summary.mm['ACTB',]
#singlecell.summary.mm['ACTC1',]
#singlecell.summary.mm['ACTA2',]

## see if the xenium can be modeled as a lienar cobmination of off targets
# model <- lm(xenium.summary.mm[g,] ~ 
#               singlecell.summary.mm[g,] + 
#               singlecell.summary.mm['ACTG1',] +
#               singlecell.summary.mm['ACTA1',] +
#               singlecell.summary.mm['ACTB',] +
#               singlecell.summary.mm['ACTC1',] +
#               singlecell.summary.mm['ACTA2',]
# )
# model

## must be non-negative though
#install.packages('nnls')
library(nnls)
B <- t(singlecell.summary.mm[offtargets,])
heatmap(cbind(xenium=xenium.summary.mm[g,], as.matrix(B)), scale='col', Colv=NA)

nnls_result <- nnls(B, xenium.summary.mm[g,])
names(nnls_result$x) <- offtargets
nnls_result$x

## not really working...seems like logic is off...
## inferred beta coefficients should either be 0 or 1...not scalars
offtargets[which(nnls_result$x==max(nnls_result$x))]


g1 <- g
g2 <- 'ACTB' 

## more correlated
plot(log10(singlecell.summary.mm[g,] + singlecell.summary.mm[g2,]), log10(xenium.summary.mm[g,]))

par(mfrow=c(1,2), mar=rep(2,4)) 
cells.have <- intersect(rownames(emb.info), colnames(singlecell$counts))
MERINGUE::plotEmbedding(emb.info[cells.have,], col = log10(singlecell$counts[g, cells.have]+1), main=paste0('single cell:', g))
cells.have2 <- intersect(rownames(emb.info), colnames(xenium$counts))
MERINGUE::plotEmbedding(emb.info[cells.have2,], col = log10(xenium$counts[g, cells.have2]+1), main=paste0('xenium:', g))

cells.have <- intersect(rownames(emb.info), colnames(singlecell$counts))
MERINGUE::plotEmbedding(emb.info[cells.have,], col = log10(singlecell$counts[g2, cells.have]+1), main=paste0('single cell:', g2))
cells.have2 <- intersect(rownames(emb.info), colnames(xenium$counts))
MERINGUE::plotEmbedding(emb.info[cells.have2,], col = log10(xenium$counts[g2, cells.have2]+1), main=paste0('xenium:', g2))



###################### visualize in visium

visium_df <- data.frame(visium$pos, gexp=
                          visium$counts[g1,] 
)
ggplot(visium_df, aes(x=x, y=y, col=gexp)) + geom_point(size=2, alpha=1) + coord_fixed() + 
  scale_color_gradient(low='lightgrey', high='red') +
  theme_void() + ggtitle(g1)

visium_df <- data.frame(visium$pos, gexp=
                          visium$counts[g2,] 
)
ggplot(visium_df, aes(x=x, y=y, col=gexp)) + geom_point(size=2, alpha=1) + coord_fixed() + 
  scale_color_gradient(low='lightgrey', high='red') +
  theme_void() + ggtitle(g2)

xenium_df <- data.frame(xenium$pos, gexp=xenium$counts[g,])
ggplot(xenium_df, aes(x=x, y=y, col=gexp)) + geom_point(size=0.001, alpha=0.5) + coord_fixed() + 
  scale_color_gradient(low='lightgrey', high='red') +
  theme_void() + ggtitle(g)

## seems like there's a lot more off target from visium though


