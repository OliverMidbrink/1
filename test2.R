library(Seurat)
library(tidyverse)
library(dplyr)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(ggridges)


# Clustering of cells and naming of cell clusters

data <- Read10X(data.dir = "E:/RandomScience/1/SC3pv3_GEX_Human_PBMC_filtered_feature_bc_matrix")
data <- CreateSeuratObject(counts = data, project='E17.5', min.cells = 3, min.features = 200)

data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth()

data <- subset(data, subset = nFeature_RNA >200 & nFeature_RNA <2500 & percent.mt <5)

data <- NormalizeData(data)


data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(data), 10)
top10


plot1 <- VariableFeaturePlot(data)
LabelPoints(plot = plot1, points = top10, repel = T)


all.genes <- rownames(data)
data <- ScaleData(data, features = all.genes)

data <- RunPCA(data, features = VariableFeatures(object = data))

DimHeatmap(data, dims = 1:6, cells = 500, balanced = T)

ElbowPlot(data)

data <- FindNeighbors(data, dims = 1:20)

data <- FindClusters(data, resolution = c(0.3, 0.5, 0.7, 1))

head(data@meta.data)

p1 <- DimPlot(data, group.by = "RNA_snn_res.0.3", label = T)
p2 <- DimPlot(data, group.by = "RNA_snn_res.0.5", label = T)
p3 <- DimPlot(data, group.by = "RNA_snn_res.0.7", label = T)
p4 <- DimPlot(data, group.by = "RNA_snn_res.1", label = T)
p1 + p2 + p3 + p4

Idents(data) <- "RNA_snn_res.0.7"

data <- RunUMAP(data, dims = 1:20)

DimPlot(data, reduction = "umap")


# Analyze cluster gene expression
cluster_id <- 7  # Replace with the cluster ID you want to analyze

# Find marker genes for the specified cluster
marker_genes <- FindMarkers(data, ident.1 = cluster_id, test.use = "MAST", min.pct = 0.1, thresh.use = 0.05)

# View the top marker genes
head(marker_genes, n = 10)
top_marker_genes <- rownames(head(marker_genes, n = 10))


# Plot the results
VlnPlot(data, features = top_marker_genes, slot = "counts", log = TRUE)
FeaturePlot(data, features = top_marker_genes)


# Load the manual marker genes list from the separate file
source("E:/RandomScience/1/marker_genes.R")

# Plot other genes
VlnPlot(data, features = marker_genes, slot = "counts", log = TRUE)

FeaturePlot(data, features = c("CD11c",  "Eomes", "Aldh1l1",
                               "Tbr1",  "Olig2", "Sox2", "Cux2", "Neurog2"))



# Find the cell types of the clusters:
cluster_assignments <- Idents(data)

unique_clusters <- as.integer(levels(cluster_assignments))

# View the list of unique clusters
unique_clusters

source("E:/RandomScience/1/find_marker_genes.R")

all_marker_genes <- find_cluster_marker_genes(unique_clusters)
print(all_marker_genes)



