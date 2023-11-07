
## FIND the marker genes
find_cluster_marker_genes <- function(cluster_ids) {
  
  # List to store marker genes for each cluster
  all_marker_genes <- list()
  
  # Cluster IDs to analyze (replace with your cluster IDs)
  cluster_ids <- cluster_ids  # Add the cluster IDs you want to analyze
  
  # Iterate through cluster IDs
  for (cluster_id in cluster_ids) {
    # Find marker genes for the specified cluster
    marker_genes <- FindMarkers(data, ident.1 = cluster_id, test.use = "MAST", min.pct = 0.1, thresh.use = 0.05)
    
    # Store the top marker genes in the list
    top_marker_genes <- rownames(head(marker_genes, n = 10))
    all_marker_genes[[as.character(cluster_id)]] <- top_marker_genes
  }
  
  # View the list of marker genes for each cluster
  print(all_marker_genes)
  
  
  return(all_marker_genes)
}