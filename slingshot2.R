library(Seurat)
library(tidyverse)
library(dplyr)
library(monocle3)
library(scran)
library(scater)
library(igraph)
library(cowplot)
library(slingshot)
library(tradeSeq)
library(uwot)
library(mclust)
library(RColorBrewer)

# Load and process the pre-treatment data
dataList <- Read10X(data.dir = "E:/RandomScience/1/SC3_v3_NextGem_DI_CellPlex_CRISPR_A549_30K_A549_Small_Pool_v2_No_Treatment_count_sample_feature_bc_matrix/sample_feature_bc_matrix/")
names(dataList)
data <- dataList$"Gene Expression"

data <- CreateSeuratObject(counts = data, project = 'Pre_Treatment2')
data <- NormalizeData(data)
data <- FindVariableFeatures(data)
data <- ScaleData(data)
data <- RunPCA(data)
data <- RunUMAP(data, dims = 1:20)
data <- FindNeighbors(data, dims = 1:20)
data <- FindClusters(data, resolution = 1)
data <- RunUMAP(data, n.neighbors = 10, dims = 1:50, spread = 2, min.dist = 0.3)


DimPlot(data, group.by = "RNA_snn_res.1")

dimred <- data@reductions$umap@cell.embeddings
clustering <- data$RNA_snn_res.1
counts <- as.matrix(data@assays$RNA@layers$counts)

lineages <- getLineages(data = dimred, clusterLabels = clustering, start.clus = "9")
lineages

pal <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"))
par(mfrow = c(1, 2))
plot(dimred[, 1:2], col = pal[clustering], cex = 0.5, pch = 16)
for (i in levels(clustering)) {
  text(mean(dimred[clustering == i, 1]), mean(dimred[clustering == i, 2]), labels = i, font = 2)
}
plot(dimred[, 1:2], col = pal[clustering], cex = 0.5, pch = 16)



for (lineage in lineages@metadata$lineages){
  last_element <- NULL
  
  for (element in lineage) {
    if (!is.null(last_element)) {
      x0 <- mean(dimred[clustering == element, 1])
      y0 <- mean(dimred[clustering == element, 2])
      
      x1 <- mean(dimred[clustering == last_element, 1])
      y1 <- mean(dimred[clustering == last_element, 2])
      
      segments(x0, y0, x1, y1, col = "black", lwd = 2)
    }
    print(element)
    
    last_element <- element
    
  }
}

for (i in levels(clustering)) {
  # Calculate the center coordinates for the white sphere
  x_center <- mean(dimred[clustering == i, 1])
  y_center <- mean(dimred[clustering == i, 2])
  
  # Draw a white circle (sphere) that is full and larger
  points(x_center, y_center, pch = 19, col = "red", cex = 3)  # Use pch = 19 and adjust cex as needed
  # Insert the text label
  text(x_center, y_center, labels = i, font = 2)
}



## NEW CODE


sce <- SingleCellExperiment(assays = List(counts = counts))
geneFilter <- apply(assays(sce)$counts,1,function(x){
  sum(x >= 3) >= 10
})
sce <- sce[geneFilter, ]


FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}
assays(sce)$norm <- FQnorm(assays(sce)$counts)

pca <- prcomp(t(log1p(assays(sce)$norm)), scale. = FALSE)
rd1 <- pca$x[,1:2]

plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)

rd2 <- uwot::umap(t(log1p(assays(sce)$norm)))
colnames(rd2) <- c('UMAP1', 'UMAP2')

plot(rd2, col = rgb(0,0,0,.5), pch=16, asp = 1)

reducedDims(sce) <- SimpleList(PCA = rd1, UMAP = rd2)



library(mclust, quietly = TRUE)


cl1 <- Mclust(rd1)$classification
colData(sce)$GMM <- cl1

library(RColorBrewer)
plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)


cl2 <- kmeans(rd1, centers = 4)$cluster
colData(sce)$kmeans <- cl2

plot(rd1, col = brewer.pal(9,"Set1")[cl2], pch=16, asp = 1)


sce <- slingshot(sce, clusterLabels = 'GMM', reducedDim = 'PCA')
summary(sce$slingPseudotime_1)

library(grDevices)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]

plot(reducedDims(sce)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')

plot(reducedDims(sce)$PCA, col = brewer.pal(9,'Set1')[sce$GMM], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')


library(tradeSeq)

# fit negative binomial GAM
sce <- fitGAM(sce)

# test for dynamic expression
ATres <- associationTest(sce)


topgenes <- rownames(ATres[order(ATres$pvalue), ])[1:250]
pst.ord <- order(sce$slingPseudotime_1, na.last = NA)
heatdata <- assays(sce)$counts[topgenes, pst.ord]
heatclus <- sce$GMM[pst.ord]

heatmap(log1p(heatdata), Colv = NA,
        ColSideColors = brewer.pal(9,"Set1")[heatclus])


