library(Seurat)
library(tidyverse)
library(dplyr)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(ggridges)

plotVariance <- function(data) {
  data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
  top10 <- head(VariableFeatures(data), 10)
  top10
  
  plot1 <- VariableFeaturePlot(data)
  LabelPoints(plot = plot1, points = top10, repel = T)
}


heatmap <- function(data){
  data <- ScaleData(data, features = all.genes)
  data <- RunPCA(data, features = VariableFeatures(object = data))
  DimHeatmap(data, dims = 1:20, cells = 500, balanced = T)
}

umap <- function(data){
  data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
  data <- subset(data, subset = nFeature_RNA >200 & nFeature_RNA <2500 & percent.mt <5)
  data <- NormalizeData(data)
  data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
  data <- ScaleData(data)
  data <- RunPCA(data, features = VariableFeatures(object = data))
  data <- RunUMAP(data, dims = 1:20)
  data <- FindNeighbors(data, dims = 1:20)
  data <- FindClusters(data, resolution = 0.3) # Adjust resolution as needed
  plot <- DimPlot(data, reduction = "umap", group.by = "seurat_clusters")
  plot
  return(plot)
}


dataPre <- Read10X(data.dir = "E:/RandomScience/1/aml035_pre_transplant_filtered_gene_bc_matrices/filtered_matrices_mex/hg19")
dataPre <- CreateSeuratObject(counts = dataPre, project='E17.5', min.cells = 3, min.features = 200)

dataPost <- Read10X(data.dir = "E:/RandomScience/1/aml035_post_transplant_filtered_gene_bc_matrices/filtered_matrices_mex/hg19")
dataPost <- CreateSeuratObject(counts = dataPost, project='E17.5', min.cells = 3, min.features = 200)

plot <- umap(dataPre)

ggsave("E:/RandomScience/1/umap_plot_pre.png", plot = plot)

