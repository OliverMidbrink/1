library(monocle3)

# The tutorial shown below and on subsequent pages uses two additional packages:
library(ggplot2)
library(dplyr)



# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, time_bin="130-170"){
  cell_ids <- which(colData(cds)[, "embryo.time.bin"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))


# Make the CDS object
cds <- load_mm_data(mat_path = "E:/RandomScience/1/data_filtered/matrix.mtx", 
                    feature_anno_path = "E:/RandomScience/1/data_filtered/features.tsv", 
                    cell_anno_path = "E:/RandomScience/1/data_filtered/barcodes.tsv")

## Step 1: Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 100)
cds <- align_cds(cds, alignment_group = "batch", residual_model_formula_str = "~ bg.300.loading + bg.400.loading + bg.500.1.loading + bg.500.2.loading + bg.r17.loading + bg.b01.loading + bg.b02.loading")

colnames(colData(cds))
# 2D plot
cds_2d <- reduce_dimension(cds,  reduction_method=c("UMAP"))
plot_cells(cds_2d, label_groups_by_cluster=FALSE,  color_cells_by = "n.umi")

cds_2d <- cluster_cells(cds_2d)
cds_2d <- learn_graph(cds)
plot_cells(cds_2d, label_groups_by_cluster=FALSE,  color_cells_by = "partition")


ciliated_genes <- c("ENSG00000142949",
                    "ATCCATTGTATACGGG-1",
                    "nhr-6",
                    "dmd-6",
                    "ceh-36",
                    "ham-1")

plot_cells(cds_2d,
           genes=ciliated_genes,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)




# Experimentation
gene_expression <- assays(cds_2d)$counts
colnames(assays(cds_2d)$counts)

genes_in_cds_2d <- colnames(assays(cds_2d)$counts)

# Check if each gene in 'ciliated_genes' is present in 'genes_in_cds_2d'
genes_not_found <- setdiff(ciliated_genes, genes_in_cds_2d)

if (length(genes_not_found) > 0) {
  cat("The following genes were not found in cds_2d:\n")
  cat(genes_not_found, sep = "\n")
} else {
  cat("All genes in ciliated_genes were found in cds_2d.\n")
}


# 3D plot
cds_3d <- reduce_dimension(cds, max_components = 3)
cds_3d <- cluster_cells(cds_3d, reduction_method=c("UMAP"))
cds_3d <- learn_graph(cds_3d)
cds_3d <- order_cells(cds_3d)

cds_3d_plot_obj <- plot_cells_3d(cds_3d, color_cells_by="pseudotime")

plot_cells_3d(cds_3d, color_cells_by="pseudotime")

