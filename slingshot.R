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

# Load and process the pre-treatment data
dataListPre <- Read10X(data.dir = "E:/RandomScience/1/SC3_v3_NextGem_DI_CellPlex_CRISPR_A549_30K_A549_Small_Pool_v2_No_Treatment_count_sample_feature_bc_matrix/sample_feature_bc_matrix/")
names(dataListPre)
dataPre <- dataListPre$"Gene Expression"

dataPre <- CreateSeuratObject(counts = dataPre, project = 'Pre_Treatment2')
dataPre <- NormalizeData(dataPre)
dataPre <- FindVariableFeatures(dataPre)
dataPre <- ScaleData(dataPre)
dataPre <- RunPCA(dataPre)
dataPre <- RunUMAP(dataPre, dims = 1:20)
dataPre <- FindNeighbors(dataPre, dims = 1:20)
dataPre <- FindClusters(dataPre, resolution = 1)

# Load and process the post-treatment data
dataListPost <- Read10X(data.dir = "E:/RandomScience/1/SC3_v3_NextGem_DI_CellPlex_CRISPR_A549_30K_A549_Small_Pool_v2_Resveratrol_count_sample_feature_bc_matrix/sample_feature_bc_matrix")
names(dataListPost)
dataPost <- dataListPost$"Gene Expression"

dataPost <- CreateSeuratObject(counts = dataPost, project = 'Post_Treatment2')
dataPost <- NormalizeData(dataPost)
dataPost <- FindVariableFeatures(dataPost)
dataPost <- ScaleData(dataPost)
dataPost <- RunPCA(dataPost)
dataPost <- RunUMAP(dataPost, dims = 1:20)
dataPost <- FindNeighbors(dataPost, dims = 1:20)
dataPost <- FindClusters(dataPost, resolution = 1)

# Integrate the datasets to align them
# There are multiple methods for integration; here we will use the standard integration workflow
anchors <- FindIntegrationAnchors(object.list = list(dataPre, dataPost))
dataCombined <- IntegrateData(anchorset = anchors)

# Perform the linear dimensional reduction and clustering on the integrated dataset
dataCombined <- ScaleData(dataCombined)
dataCombined <- FindVariableFeatures(dataCombined, nfeatures = 2000)
dataCombined <- RunPCA(dataCombined)
dataCombined <- FindNeighbors(dataCombined)
dataCombined <- FindClusters(dataCombined, resolution = 1)
dataCombined <- RunUMAP(dataCombined, n.neighbors = 10, dims = 1:50, spread = 2, min.dist = 0.3)
# dataCombined <- NormalizeData(dataCombined)

DimPlot(dataCombined, group.by = "integrated_snn_res.1")
dimred <- dataCombined@reductions$umap@cell.embeddings
clustering <- dataCombined$integrated_snn_res.1
counts.Pre_treatment2 <- GetAssayData(object = dataCombined, assay = "RNA", slot = "counts.Pre_Treatment2")
counts.Post_treatment2 <- GetAssayData(object = dataCombined, assay = "RNA", slot = "counts.Post_Treatment2")

pal <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"))

lineages <- getLineages(data = dimred, clusterLabels = clustering, start.clus = "9")



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


############# Diff Exp #####################################
curves <- getCurves(lineages, approx_points = 300, thresh = 0.01, stretch = 0.8, allow.breaks = FALSE, shrink = 0.99)
curves@metadata$curves
c
plot(dimred, col = pal[clustering], asp = 1, pch = 16)
lines(curves, lwd = 3, col = "black")

dim(counts.Post_treatment2)


sce <- fitGAM(counts = as.matrix(counts.Post_treatment2), sds = curves)

plotGeneCount(curves, filt_counts, clusters = clustering, models = sce)


plot_differential_expression <- function(feature_id) {
  feature_id <- pseudotime_association %>% filter(pvalue < 0.05) %>% top_n(1, -waldStat) %>% pull(feature_id)
  cowplot::plot_grid(plotGeneCount(curves, filt_counts, gene = feature_id[1], clusters = clustering, models = sce) + ggplot2::theme(legend.position = "none"), 
                     plotSmoothers(sce, as.matrix(counts), gene = feature_id[1]))
}


pseudotime_association <- associationTest(sce)
pseudotime_association$fdr <- p.adjust(pseudotime_association$pvalue, method = "fdr")
pseudotime_association <- pseudotime_association[order(pseudotime_association$pvalue), ]
pseudotime_association$feature_id <- rownames(pseudotime_association)






####################### STOP BEFORE THIS POINT ############################



# Visualize the combined UMAP plot with clusters
DimPlot(dataCombined, reduction = 'umap', group.by = 'seurat_clusters')

# To highlight differences pre- and post-treatment, use the "split.by" argument
DimPlot(dataCombined, reduction = 'umap', split.by = 'orig.ident')

# To visualize the expression of a marker gene pre- and post-treatment
FeaturePlot(dataCombined, features = 'ALX1', split.by = 'orig.ident')




# Markers
cluster_markers <- FindAllMarkers(dataCombined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster_markers$pct_diff <- cluster_markers$pct.1 - cluster_markers$pct.2
cluster_markers$pct.2.weighted <- -5 * cluster_markers$pct.2
cluster_markers$unique_score <- cluster_markers$pct.2.weighted + cluster_markers$pct.1


top_markers <- cluster_markers[order(-cluster_markers$pct_diff), ]
top_10_markers <- head(top_markers, 10)
print(top_10_markers)

top_unique_markers <- cluster_markers[order(-cluster_markers$unique_score), ]
top_10_unique_markers <- head(top_unique_markers, 10)
print(top_10_unique_markers)



DoHeatmap(dataCombined, features = head(cluster_markers[order(-cluster_markers$unique_score), ]$gene))

## Get the cell types
source("E:/RandomScience/1/cell_markers.R")

# Calculate the scores for each cell type
cell_type_scores <- sapply(names(cell_markers), function(cell_type) {
  markers <- cell_markers[[cell_type]]
  scores <- rowMeans(dataPre@assays$RNA@data[markers, ])
  return(scores)
}, simplify = "data.frame")

# Add the scores to the metadata of the Seurat object
for (cell_type in names(cell_markers)) {
  dataPre@meta.data[[paste0(cell_type, "_score")]] <- cell_type_scores[[cell_type]]
}

# Check the top cell types per cell based on the scores
cell_type_max_scores <- apply(cell_type_scores, 1, which.max)
dataPre@meta.data$max_cell_type <- names(cell_markers)[cell_type_max_scores]

# Let's view the cell types assigned based on the highest score
table(dataPre@meta.data$max_cell_type)

annotate_cluster <- function(cluster_id, cluster_markers, cell_markers) {
  # Extract the marker genes for the current cluster
  cluster_marker_genes <- cluster_markers$gene[cluster_markers$cluster == cluster_id]
  
  # Initialize a list to store the results
  matched_cell_types <- list()
  
  # Loop over the list of known cell type markers
  for (cell_type in names(cell_markers)) {
    # Find the intersection between cluster markers and cell type markers
    matching_markers <- intersect(cluster_marker_genes, cell_markers[[cell_type]])
    
    # If there are matching markers, add them to the result list
    if (length(matching_markers) > 0) {
      matched_cell_types[[cell_type]] <- matching_markers
    }
  }
  
  # Return the list of matched cell types and their markers
  return(matched_cell_types)
}

# Annotate each cluster
cluster_annotations <- lapply(unique(cluster_markers$cluster), function(cluster_id) {
  annotate_cluster(cluster_id, cluster_markers, cell_markers)
})

# If you want to print the annotations for each cluster
for (i in seq_along(cluster_annotations)) {
  cat("Cluster", names(cluster_annotations[i]), ":\n")
  print(cluster_annotations[[i]])
  cat("\n")
}




write.csv(cluster_markers, file = "cluster_markers.csv", row.names = FALSE)




# Find markers for a specific cluster assumed to be myeloid
myeloid.cluster.markers <- FindMarkers(dataPre, ident.1 = "3", min.pct = 0.25)

# View the top 10 markers for the cluster of interest
features_ <- rownames(head(myeloid.cluster.markers, 4))

FeaturePlot(dataCombined, features = features_, split.by = 'orig.ident')




dataCombined <- NormalizeData(dataCombined)

# Find variable features
dataCombined <- FindVariableFeatures(dataCombined, selection.method = "vst", nfeatures = 2000)

# Scale the data
dataCombined <- ScaleData(dataCombined, features = VariableFeatures(object = dataCombined))

# Run PCA
dataCombined <- RunPCA(dataCombined, features = VariableFeatures(object = dataCombined))

# Run UMAP using the first 10 PCs
dataCombined <- RunUMAP(dataCombined, reduction = "pca", dims = 1:10)

# Plot UMAP
DimPlot(dataCombined, reduction = "umap", group.by = "seurat_clusters")




# Monocle3
cdsPre <- new_cell_data_set(dataListPre$"Gene Expression")

# Create the CDS object for post-treatment data
cdsPost <- new_cell_data_set(dataListPost$"Gene Expression")

cdsCombined <- combine_cds(list(cdsPre, cdsPost))

cds_3d <- preprocess_cds(cdsCombined, method = "PCA", num_dim = 10)
cds_3d <- reduce_dimension(cds_3d, max_components = 3)
cds_3d <- cluster_cells(cds_3d)
cds_3d <- learn_graph(cds_3d)
cds_3d <- order_cells(cds_3d)

cds_3d_plot_obj <- plot_cells_3d(cds_3d)
cds_3d_plot_obj

cds_2d <- preprocess_cds(cdsCombined, method = "PCA", num_dim = 10)
cds_2d <- reduce_dimension(cds_2d, reduction_method = 'UMAP',
                           umap.metric = 'euclidean',
                           umap.n_neighbors = 30,
                           umap.min_dist = 0.3)
cds_2d <- cluster_cells(cds_2d)
cds_2d <- learn_graph(cds_2d)
plot_cells(cds_2d)

plot_cells(cds_2d, color_cells_by="cao_cell_type")


