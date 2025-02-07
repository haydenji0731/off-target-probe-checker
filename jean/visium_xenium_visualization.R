dir <- '~/OneDrive - Johns Hopkins/Data_Public/xenium_data/'

############## read xenium data
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

########### calculate Moran's I

library(SpatialExperiment)

pos <- as.matrix(cbind(x=pos.info$x_centroid, y=pos.info$y_centroid))
rownames(pos) <- pos.info$cell_id
se <- SpatialExperiment(
  assays = list(counts = counts),
  spatialCoords = pos[colnames(counts),]
)

## make same solution as visium
library(SEraster)
res <- 100
rastGexp <- SEraster::rasterizeGeneExpression(se,
                                              assay_name="counts",
                                              resolution = res,
                                              square=FALSE,
                                              n_threads = 1,
                                              fun = "sum")

## get SVGs
rastpos <- spatialCoords(rastGexp)
pixelval <- assay(rastGexp, 'pixelval')
# library(MERINGUE)
# w <- MERINGUE::getSpatialNeighbors(rastpos, filterDist = res)
# I_xenium <- getSpatialPatterns(pixelval, w)
# 
xenium <- list(pos=pos,
               counts=counts,
               rastpos = rastpos,
               pixelval = pixelval)

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

## limit to same genes as Xenium
shared.genes <- intersect(rownames(visium$counts), rownames(xenium$counts))
shared.genes

# library(MERINGUE)
# w <- MERINGUE::getSpatialNeighbors(pos, filterDist = res)
# I_visium <- getSpatialPatterns(gexp, w)
# 
# visium <- list(pos=pos,
#                counts=gexp,
#                I=I_visium)

######### visualize
library(ggplot2)

#g <- shared.genes[1]
#g <- 'CEACAM8'

sort(rowSums(xenium$counts))
g <- 'KRT7'
g <- 'APOBEC3B'
g <- 'C1QA'
g <- 'PLD4'
g <- 'CEACAM6'

xenium_df <- data.frame(xenium$pos, gexp=xenium$counts[g,])
ggplot(xenium_df, aes(x=x, y=y, col=gexp)) + geom_point(size=0.001, alpha=0.5) + coord_fixed() + 
  scale_color_gradient(low='lightgrey', high='red') +
  theme_void() + ggtitle(g)

xenium_df2 <- data.frame(xenium$rastpos, gexp=xenium$pixelval[g,])
ggplot(xenium_df2, aes(x=x, y=y, col=gexp)) + geom_point(size=1, alpha=1) + coord_fixed() + 
  scale_color_gradient(low='lightgrey', high='red') +
  theme_void() + ggtitle(g)

visium_df <- data.frame(visium$pos, gexp=visium$counts[g,])
ggplot(visium_df, aes(x=x, y=y, col=gexp)) + geom_point(size=1, alpha=1) + coord_fixed() + 
  scale_color_gradient(low='lightgrey', high='red') +
  theme_void() + ggtitle(g)


###################### CEACAM8 probes also bind CEACAM7 also CEACAM5
g <- g1 <- 'CEACAM8'
g2 <- 'CEACAM5' ## much more highly expressed
g3 <- 'CEACAM7'

visium_df <- data.frame(visium$pos, gexp=
                          visium$counts[g1,] 
)
ggplot(visium_df, aes(x=x, y=y, col=gexp)) + geom_point(size=1, alpha=1) + coord_fixed() + 
  scale_color_gradient(low='lightgrey', high='red') +
  theme_void() + ggtitle(g1)

visium_df <- data.frame(visium$pos, gexp=
                          visium$counts[g2,] 
)
ggplot(visium_df, aes(x=x, y=y, col=gexp)) + geom_point(size=1, alpha=1) + coord_fixed() + 
  scale_color_gradient(low='lightgrey', high='red') +
  theme_void() + ggtitle(g2)

visium_df <- data.frame(visium$pos, gexp=
                          visium$counts[g3,] 
)
ggplot(visium_df, aes(x=x, y=y, col=gexp)) + geom_point(size=1, alpha=1) + coord_fixed() + 
  scale_color_gradient(low='lightgrey', high='red') +
  theme_void() + ggtitle(g3)

visium_df <- data.frame(visium$pos, gexp=
                          visium$counts[g1,] +
                          visium$counts[g2,] +
                          visium$counts[g3,]
                          )
ggplot(visium_df, aes(x=x, y=y, col=gexp)) + geom_point(size=1, alpha=1) + coord_fixed() + 
  scale_color_gradient(low='lightgrey', high='red') +
  theme_void() + ggtitle('composite')

xenium_df2 <- data.frame(xenium$rastpos, gexp=xenium$pixelval[g,])
ggplot(xenium_df2, aes(x=x, y=y, col=gexp)) + geom_point(size=2, alpha=1) + coord_fixed() + 
  scale_color_gradient(low='lightgrey', high='red') +
  theme_void() + ggtitle(g)


###### compare magnitudes
visium_total <- rowSums(visium$counts[shared.genes,])
xenium_total <- rowSums(xenium$counts[shared.genes,])

setdiff(rownames(xenium$counts), shared.genes)[1:6] %in% rownames(visium$counts)

df <- data.frame(visium_total, xenium_total, gene = shared.genes)
ggplot(df, aes(x = visium_total, y = xenium_total)) + geom_point() 

ggplot(df, aes(x = visium_total, y = xenium_total)) + geom_point() + 
  scale_x_log10() + scale_y_log10()

ggplot(df, aes(x = visium_total, y = xenium_total, label=gene)) + geom_point() + 
  scale_x_log10() + scale_y_log10() + ggrepel::geom_label_repel()

visium_total[g]
xenium_total[g]

which(visium_total < 100 & xenium_total > 100)

