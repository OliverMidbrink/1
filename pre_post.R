library(Seurat)
library(tidyverse)
library(dplyr)
library(monocle3)


# Load and process the pre-treatment data
dataPre <- Read10X(data.dir = "E:/RandomScience/1/aml035_pre_transplant_filtered_gene_bc_matrices/filtered_matrices_mex/hg19")
dataPre <- CreateSeuratObject(counts = dataPre, project = 'Pre_Treatment')
dataPre <- NormalizeData(dataPre)
dataPre <- FindVariableFeatures(dataPre)
dataPre <- ScaleData(dataPre)
dataPre <- RunPCA(dataPre)
dataPre <- RunUMAP(dataPre, dims = 1:20)
dataPre <- FindNeighbors(dataPre, dims = 1:20)
dataPre <- FindClusters(dataPre)

# Load and process the post-treatment data
dataPost <- Read10X(data.dir = "E:/RandomScience/1/aml035_post_transplant_filtered_gene_bc_matrices/filtered_matrices_mex/hg19")
dataPost <- CreateSeuratObject(counts = dataPost, project = 'Post_Treatment')
dataPost <- NormalizeData(dataPost)
dataPost <- FindVariableFeatures(dataPost)
dataPost <- ScaleData(dataPost)
dataPost <- RunPCA(dataPost)
dataPost <- RunUMAP(dataPost, dims = 1:20)
dataPost <- FindNeighbors(dataPost, dims = 1:20)
dataPost <- FindClusters(dataPost)

# Integrate the datasets to align them
# There are multiple methods for integration; here we will use the standard integration workflow
anchors <- FindIntegrationAnchors(object.list = list(dataPre, dataPost))
dataCombined <- IntegrateData(anchorset = anchors)

# Perform the linear dimensional reduction and clustering on the integrated dataset
dataCombined <- ScaleData(dataCombined)
dataCombined <- RunPCA(dataCombined)
dataCombined <- RunUMAP(dataCombined, dims = 1:30)
dataCombined <- FindNeighbors(dataCombined, dims = 1:30)
dataCombined <- FindClusters(dataCombined)

# Visualize the combined UMAP plot with clusters
DimPlot(dataCombined, reduction = 'umap', group.by = 'seurat_clusters')

# To highlight differences pre- and post-treatment, use the "split.by" argument
DimPlot(dataCombined, reduction = 'umap', split.by = 'orig.ident')

# To visualize the expression of a marker gene pre- and post-treatment
FeaturePlot(dataCombined, features = 'FCER1A', split.by = 'orig.ident')





myeloid.markers <- FindAllMarkers(dataPre, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Find markers for a specific cluster assumed to be myeloid
myeloid.cluster.markers <- FindMarkers(dataPre, ident.1 = "3", min.pct = 0.25)

# View the top 10 markers for the cluster of interest
features_ <- rownames(head(myeloid.cluster.markers, 4))

FeaturePlot(dataCombined, features = features_, split.by = 'orig.ident')


Idents(dataCombined) <- "RNA_snn_res.0.7"
DimPlot(dataCombined, reduction = "umap")


