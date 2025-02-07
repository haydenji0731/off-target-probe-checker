## cleaned up analysis, export figures
## focus on coclustering with 3' sequencing since it's sequencing and more conventional
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
  ggrepel::geom_label_repel() + theme_classic()


################### do harmonized clustering
## https://portals.broadinstitute.org/harmony/mudan.html
cd <- cbind(xenium$counts[shared.genes,], singlecell$counts[shared.genes,])
meta <- c(rep('xenium', ncol(xenium$counts)), rep('singlecell', ncol(singlecell$counts)))
names(meta) <- c(colnames(xenium$counts), colnames(singlecell$counts))
meta <- factor(meta)
table(meta)

## filter out single cells without any counts
vi <- colSums(cd) > 1
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
colnames(emb.info) <- c('UMAP1', 'UMAP2')
head(emb.info)
save(harmonized, emb.info, file="embedding.RData")

# Plot
par(mfrow=c(1,1), mar=rep(2,4))
MERINGUE::plotEmbedding(emb.info, groups=meta, 
                        show.legend=TRUE, xlab=NA, ylab=NA, 
                        main='Harmonized UMAP Embedding',
                        verbose=FALSE)

## make ggplot version
df <- data.frame(emb.info, groups=meta)
df <- df[order(df$groups, decreasing=TRUE),]
ggplot(df, aes(x=UMAP1, y=UMAP2, col=groups)) + geom_point(size = 0.1, alpha=0.2) +
  theme_classic() + coord_fixed() + xlim(c(-10, 10)) + ylim(c(-10, 10))

par(mfrow=c(1,2), mar=rep(2,4)) 
g <- 'ACTG2'
cells.have <- intersect(rownames(emb.info), colnames(singlecell$counts))
MERINGUE::plotEmbedding(emb.info[cells.have,], col = log10(singlecell$counts[g, cells.have]+1), main=paste0('single cell:', g))
points(emb.info[names(which(log10(singlecell$counts[g, cells.have]+1) > 0)),], col='red', pch=16)

cells.have2 <- intersect(rownames(emb.info), colnames(xenium$counts))
MERINGUE::plotEmbedding(emb.info[cells.have2,], col = log10(xenium$counts[g, cells.have2]+1), main=paste0('xenium:', g))
#points(emb.info[names(which(log10(xenium$counts[g, cells.have2]+1) > 0)),], col='red', pch=16)

g <- 'ACTG2'
cells.have <- intersect(rownames(emb.info), colnames(singlecell$counts))
df <- data.frame(emb.info[cells.have,], gene = log10(singlecell$counts[g, cells.have]+1))
df <- df[order(df$gene),]
ggplot(df, aes(x=UMAP1, y=UMAP2, col=gene)) + geom_point(size = 0.1, alpha=0.2) +
  theme_classic() + coord_fixed() + scale_color_gradient(low = 'lightgrey', high='red') +
  ggtitle(g) + xlim(c(-10, 10)) + ylim(c(-10, 10))

cells.have2 <- intersect(rownames(emb.info), colnames(xenium$counts))
df <- data.frame(emb.info[cells.have2,], gene = log10(xenium$counts[g, cells.have2]+1))
df <- df[order(df$gene),]
ggplot(df, aes(x=UMAP1, y=UMAP2, col=gene)) + geom_point(size = 0.1, alpha=0.2) +
  theme_classic() + coord_fixed() + scale_color_gradient(low = 'lightgrey', high='red') +
  ggtitle(g) + xlim(c(-10, 10)) + ylim(c(-10, 10))


##################### do our own clustering
# getApproxComMembership <- function (mat, k, nsubsample = nrow(mat) * 0.5, method = igraph::cluster_walktrap,
#                                     seed = 0, vote = FALSE, verbose = TRUE)
# {
#   if (verbose) {
#     print(paste0("Subsampling from ", nrow(mat), " cells to ",
#                  nsubsample, " ... "))
#   }
#   set.seed(seed)
#   subsample <- sample(rownames(mat), nsubsample)
#   if (verbose) {
#     print("Identifying cluster membership for subsample ... ")
#   }
#   pcs.sub <- mat[subsample, ]
#   com.sub <- MERINGUE::getClusters(pcs.sub, k = k, method = method)
#   if (verbose) {
#     print("Imputing cluster membership for rest of cells ... ")
#   }
#   if (vote) {
#     data <- mat[subsample, ]
#     query <- mat[setdiff(rownames(mat), subsample), ]
#     knn <- RANN::nn2(data, query, k = k)[[1]]
#     rownames(knn) <- rownames(query)
#     com.nonsub <- unlist(apply(knn, 1, function(x) {
#       nn <- rownames(data)[x]
#       nn.com <- com.sub[nn]
#       return(names(sort(table(nn.com), decreasing = TRUE)[1]))
#     }))
#     com.all <- factor(c(com.sub, com.nonsub)[rownames(mat)])
#   }
#   else {
#     df.sub <- data.frame(celltype = com.sub, pcs.sub)
#     model <- MASS::lda(celltype ~ ., data = df.sub)
#     df.all <- data.frame(mat)
#     model.output <- stats::predict(model, df.all)
#     com.all <- model.output$class
#     names(com.all) <- rownames(df.all)
#     if (verbose) {
#       print("Model accuracy for subsample ...")
#       print(table(com.all[names(com.sub)] == com.sub))
#     }
#   }
#   return(com.all)
# }
# com <- getApproxComMembership(harmonized, k=5, method=igraph::cluster_infomap, nsubsample = nrow(harmonized)*0.75)

## honestly simple k-means clustering may work just as well here and be more reproducible
set.seed(0)
com <- kmeans(harmonized, centers = 50)$cluster
com <- as.factor(com)
head(com)

par(mfrow=c(1,2), mar=rep(2,4)) 
cells.have <- intersect(rownames(emb.info), colnames(singlecell$counts))
MERINGUE::plotEmbedding(emb.info[cells.have,], groups=com[cells.have], main='single cell clusters', mark.clusters = TRUE)
cells.have2 <- intersect(rownames(emb.info), colnames(xenium$counts))
MERINGUE::plotEmbedding(emb.info[cells.have2,], groups=com[cells.have2], main='xenium clusters')

cells.have <- intersect(rownames(emb.info), colnames(singlecell$counts))
df <- data.frame(emb.info[cells.have,], com=com[cells.have])
ggplot(df, aes(x=UMAP1, y=UMAP2, col=com)) + geom_point(size = 0.1, alpha=0.2) +
  theme_classic() + coord_fixed() + ggtitle('single cells') + xlim(c(-10, 10)) + ylim(c(-10, 10))

cells.have2 <- intersect(rownames(emb.info), colnames(xenium$counts))
df <- data.frame(emb.info[cells.have2,], com=com[cells.have2])
ggplot(df, aes(x=UMAP1, y=UMAP2, col=com)) + geom_point(size = 0.1, alpha=0.2) +
  theme_classic() + coord_fixed() + ggtitle('xenium') + xlim(c(-10, 10)) + ylim(c(-10, 10))

df <- data.frame(emb.info, com=com)
ggplot(df, aes(x=UMAP1, y=UMAP2, col=com)) + geom_point(size = 0.1, alpha=0.2) +
  theme_classic() + coord_fixed() + ggtitle('all') + xlim(c(-10, 10)) + ylim(c(-10, 10))


####### differential expression
cells.have2 <- intersect(rownames(emb.info), colnames(xenium$counts))
dg <- MERINGUE::getDifferentialGenes(xenium$counts[, cells.have2], com[cells.have2])
dg.sig <- dg$`34`
dg.sig <- dg.sig[dg.sig$p.adj < 0.05,]
dg.sig <- dg.sig[order(dg.sig$Z, decreasing=TRUE),]
head(dg.sig, n=30)

par(mfrow=c(1,2), mar=rep(2,4)) 
g <- 'TACSTD2'
cells.have <- intersect(rownames(emb.info), colnames(singlecell$counts))
MERINGUE::plotEmbedding(emb.info[cells.have,], col = log10(singlecell$counts[g, cells.have]+1), main=paste0('single cell:', g))
#points(emb.info[names(which(log10(singlecell$counts[g, cells.have]+1) > 0)),], col='red', pch=16)

cells.have2 <- intersect(rownames(emb.info), colnames(xenium$counts))
MERINGUE::plotEmbedding(emb.info[cells.have2,], col = log10(xenium$counts[g, cells.have2]+1), main=paste0('xenium:', g))
#points(emb.info[names(which(log10(xenium$counts[g, cells.have2]+1) > 0)),], col='red', pch=16)


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



m1 <- as.matrix(singlecell.summary.mm[shared.genes,])
#m1 <- scale(m1)
m1 <- scale(t(m1))
m1[is.na(m1)] <- 0
heatmap(t(m1), scale='none')

m2 <- as.matrix(xenium.summary.mm[shared.genes,])
#m2 <- scale(m2)
m2 <- scale(t(m2))
m2[is.na(m2)] <- 0
heatmap(t(m2), scale='none')

## use same clustering
#hc <- hclust(as.dist(1-cor(t(m2+m1))), method='ward.D2')
hc <- hclust(dist(m1))
plot(hc)
foo = 1-cor(m1,use = 'pairwise.complete.obs')
foo[is.na(foo)] <- 0
rc <- hclust(as.dist(foo), method='complete')
plot(rc)  

heatmap(m1[hc$labels[hc$order], rc$labels[rc$order]], Rowv=NA, Colv=NA, scale='none')
heatmap(m2[hc$labels[hc$order], rc$labels[rc$order]], Rowv=NA, Colv=NA, scale='none')

od <- order(rowSums(m1))
heatmap(m1[od, rc$labels[rc$order]], Rowv=NA, Colv=NA, scale='none')
heatmap(m2[od, rc$labels[rc$order]], Rowv=NA, Colv=NA, scale='none')

## maybe plit by genes with offtarget vs no offtarget?
sub <- c('ACTG2', 'CEACAM8', 'TUBB2B')
heatmap(m1[od, sub], Rowv=NA, Colv=NA, scale='none')
heatmap(m2[od, sub], Rowv=NA, Colv=NA, scale='none')



#heatmap(m1-m2, Rowv=as.dendrogram(rc), Colv=as.dendrogram(hc), scale='row')
#par(mfrow=c(1,1))
#barplot(sort(colSums(m1-m2)))


########### focus on a gene
g <- 'ACTG2'
singlecell.summary.mm[g,]
xenium.summary.mm[g,]

## until calebs predictions come...use regex orthologs
offtargets <- setdiff(rownames(singlecell.summary.mm)[grepl('^ACT', rownames(singlecell.summary.mm))], g)
#offtargets <- setdiff(rownames(singlecell.summary.mm), g)
#offtargets <- c('POTEE', 'POTEJ','ACTA2', 'ACTA2-AS1')
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
#library(nnls)
#B <- t(singlecell.summary.mm[offtargets,])
#heatmap(cbind(xenium=xenium.summary.mm[g,], as.matrix(B)), scale='none', Colv=NA)
#heatmap(cbind(xenium=xenium.summary.mm[g,], as.matrix(B)), scale='col', Colv=NA)

#nnls_result <- nnls(B, xenium.summary.mm[g,])
#names(nnls_result$x) <- offtargets
#nnls_result$x

## not really working...seems like logic is off...
## inferred beta coefficients should either be 0 or 1...not scalars
#offtargets[which(nnls_result$x==max(nnls_result$x))]

x = xenium.summary.mm[g,]
A = as.matrix(t(singlecell.summary.mm[c(g, offtargets),]))
x
A

## doesn't work if not perfect
# library(lpSolve)
# 
# # Set up binary constraint (each column used or not)
# f.obj <- rep(1, ncol(A))  # Minimize number of columns used
# f.con <- A                # Constraint: A * selection = x
# f.dir <- rep("=", nrow(A))
# f.rhs <- x
# # Solve the integer program
# solution <- lp("min", f.obj, f.con, f.dir, f.rhs, all.bin = TRUE)
# solution$solution

# Run LASSO (alpha=1 ensures L1 regularization)
library(glmnet)
lasso_model <- cv.glmnet(A, x, alpha=2, lower.limits=0)  # Enforce non-negative weights
best_lambda <- lasso_model$lambda.min
weights <- coef(lasso_model, s=best_lambda)[-1]  # Remove intercept
selected_columns <- which(weights > 0)
selected_weights <- weights[selected_columns]
print(selected_columns)  # Indices of chosen columns
print(selected_weights)  # Corresponding weights

colnames(A)[selected_columns]
 
##### double check
g1 <- g
g2 <- 'ACTL7B' 

## more correlated
plot(log10(singlecell.summary.mm[g,] + singlecell.summary.mm[g2,]), log10(xenium.summary.mm[g,]))

par(mfrow=c(1,2), mar=rep(2,4)) 
cells.have <- intersect(rownames(emb.info), colnames(singlecell$counts))
MERINGUE::plotEmbedding(emb.info[cells.have,], col = log10(singlecell$counts[g, cells.have]+1), main=paste0('single cell:', g))
cells.have2 <- intersect(rownames(emb.info), colnames(xenium$counts))
MERINGUE::plotEmbedding(emb.info[cells.have2,], col = log10(xenium$counts[g, cells.have2]+1), main=paste0('xenium:', g))

cells.have <- intersect(rownames(emb.info), colnames(singlecell$counts))
MERINGUE::plotEmbedding(emb.info[cells.have,], col = log10(singlecell$counts[g2, cells.have]+1), main=paste0('single cell:', g2))

cells.have <- intersect(rownames(emb.info), colnames(singlecell$counts))
MERINGUE::plotEmbedding(emb.info[cells.have,], col = log10(singlecell$counts[g2, cells.have]+
                                                             singlecell$counts[g, cells.have]+1), main=paste0('combined ', g, g2))





###################### visualize in visium
############### read in Visium data
library(Matrix)

## gene expression
gexp <- Matrix::readMM(paste0(dir, '/Visium_xenium_section/filtered_feature_bc_matrix/matrix.mtx.gz'))
barcodes <- read.csv(paste0(dir, 'Visium_xenium_section/filtered_feature_bc_matrix/barcodes.tsv.gz'), header=FALSE)
features <- read.table(paste0(dir, 'Visium_xenium_section/filtered_feature_bc_matrix/features.tsv.gz'), header=FALSE)
colnames(gexp) <- barcodes[,1]
rownames(gexp) <- features[,2]

## spatial positions of pixels
posinfo <- read.csv(paste0(dir, 'Visium_xenium_section/spatial/tissue_positions.csv'), header=FALSE)
pos <- posinfo[,5:6]
rownames(pos) <- posinfo[,1]
### restrict to same set
pos <- pos[colnames(gexp),]
colnames(pos) <- c('x', 'y')
pos[,1] <- as.numeric(pos[,1])
pos[,2] <- -as.numeric(pos[,2])  ## flip

visium <- list(pos=pos, counts=gexp)

##### visualize
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

visium_df <- data.frame(visium$pos, gexp=
                          visium$counts[g2,] + visium$counts[g1,] 
)
ggplot(visium_df, aes(x=x, y=y, col=gexp)) + geom_point(size=2, alpha=1) + coord_fixed() + 
  scale_color_gradient(low='lightgrey', high='red') +
  theme_void() + ggtitle('combined')


xenium_df2 <- data.frame(xenium$rastpos, gexp=xenium$pixelval[g,])
ggplot(xenium_df2, aes(x=x, y=y, col=gexp)) + geom_point(size=2, alpha=1) + coord_fixed() + 
  scale_color_gradient(low='lightgrey', high='red') +
  theme_void() + ggtitle(g)


## seems like there's a lot more off target from visium though


