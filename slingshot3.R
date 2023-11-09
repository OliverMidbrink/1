library(Seurat)
library(SingleCellExperiment)
library(slingshot)
library(scran)

#setwd('/Users/oliver/Documents/Random/RandomScience/1')
setwd('E:/RandomScience/1')

# Step 1: Load data and create a Seurat object
dataList <- Read10X(data.dir = "./SC3_v3_NextGem_DI_CellPlex_CRISPR_A549_30K_A549_Small_Pool_v2_No_Treatment_count_sample_feature_bc_matrix/sample_feature_bc_matrix/")
names(dataList)
data <- dataList$"Gene Expression"

seurat_object <- CreateSeuratObject(counts = data, project="Pre_Treatment2")

# Step 2: Normalize and identify clusters using Seurat
seurat_object <- NormalizeData(seurat_object)
seurat_object <- FindVariableFeatures(seurat_object)
seurat_object <- ScaleData(seurat_object)
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
seurat_object <- FindNeighbors(seurat_object, dims = 1:20)
seurat_object <- FindClusters(seurat_object)

# Step 3: Convert to SingleCellExperiment object
sce_object <- as.SingleCellExperiment(seurat_object)

# Step 4: Run Slingshot
# Ensure you have a reducedDims entry in your sce_object, here we assume PCA has been stored
# You might need to set the reducedDims(sce_object) from the Seurat object PCA results
# Also, make sure `seurat_clusters` is the name of the cluster column in sce_object colData
sce_object <- slingshot(sce_object, clusterLabels = "seurat_clusters", reducedDim = 'PCA')

# The fitGAM function is part of the tradeSeq package and is used for fitting Generalized Additive Models
# It is not part of the Seurat or Slingshot workflow and requires additional context to be properly integrated
# Assuming you have already installed tradeSeq and have lineages identified, you could proceed
library(tradeSeq)
sce_object <- fitGAM(sce_object)

# Assuming 'gam_results' is the result of fitGAM:
saveRDS(sce_object, file = "./sce_object.rds")


