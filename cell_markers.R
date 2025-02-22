cell_markers <- list(
  T_cells = c("CD3D", "CD3E", "CD3G", "CD2", "CD5", "CD7", "CD28", "CD45RA", "CD45RO", "CCR7"),
  B_cells = c("CD19", "CD20", "CD79A", "CD22", "CD24", "CD27", "CD38", "CD72", "CD40", "HLA-DR"),
  Natural_Killer_cells = c("NCAM1", "KIR2DL1", "KIR3DL1", "NKG2D", "NKG2A", "NKp30", "NKp46", "CD16", "CD56", "CD57"),
  Monocytes = c("CD14", "CD64", "CSF1R", "CD11b", "CD13", "CD33", "CD115", "CD116", "HLA-DR", "CCR2"),
  Macrophages = c("CD68", "CD80", "CD163", "CD14", "CD64", "CD32", "CD206", "CD209", "HLA-DRA", "CCR5"),
  Dendritic_cells = c("CD11c", "HLA-DRA", "CD86", "CD1c", "CD207", "CD209", "CD40", "CD80", "CD83", "CD123"),
  Neutrophils = c("CD66b", "CEACAM8", "FCGR3B", "CD15", "CD16", "CD11b", "CD62L", "CD63", "CD66a", "CD102"),
  Eosinophils = c("SIGLEC8", "IL5RA", "CCR3", "CD125", "CD9", "CD44", "CD69", "CD23", "CCR1", "CCR5"),
  Basophils = c("CD123", "FCER1A", "CPA3", "CD203c", "CD11b", "CD13", "CD44", "CD69", "CD123", "HLA-DR"),
  Megakaryocytes = c("CD41", "CD61", "VWF", "CD42b", "CD42d", "CD36", "GP1BA", "GP9", "GPIX", "CD105"),
  Erythrocytes = c("GYPA", "HBA1", "HBB", "HBA2", "HBG1", "HBE1", "CD235a", "CD233", "CD240D", "CD242"),
  Mesenchymal_stem_cells = c("NT5E", "THY1", "ENG", "CD105", "CD73", "CD90", "CD44", "CD29", "CD166", "CD106"),
  Hematopoietic_stem_cells = c("CD34", "CD90", "CD117", "CD133", "CD135", "CD49f", "CD38", "CD45RA", "CD123", "CD10"),
  Endothelial_cells = c("PECAM1", "VWF", "CDH5", "CD34", "CD31", "CD144", "CD102", "CD105", "CD146", "CD141"),
  Epithelial_cells = c("EPCAM", "CDH1", "KRT18", "KRT8", "KRT19", "CD24", "CD29", "CD44", "CD166", "CD326"),
  Fibroblasts = c("DCN", "FSP1", "COL1A1", "CD90", "PDGFRB", "FAP", "CD13", "CD44", "CD106", "CD29"),
  Adipocytes = c("ADIPOQ", "PPARG", "LEP", "FABP4", "PLIN1", "ADIRF", "CD36", "CD29", "CDH1", "GLUT4"),
  Myocytes = c("ACTA1", "MYH2", "MB", "ACTC1", "MYH7", "MYL2", "CKM", "TNNT1", "TNNI1", "MYOG"),
  Hepatocytes = c("ALB", "AAT", "AFP", "CYP3A4", "CYP2E1", "HNF4A", "ALAT", "ASAT", "GGT1", "FAH"),
  Neurons = c("MAP2", "NEFL", "GAP43", "NFH", "NFM", "SYN1", "RBFOX3", "SYP", "GAD1", "VACHT"),
  Astrocytes = c("GFAP", "SLC1A3", "ALDH1L1", "AQP4", "GLT1", "GLAST", "S100B", "SPARC", "EAAT1", "EAAT2"),
  Oligodendrocytes = c("MBP", "MOG", "CNP", "PLP1", "OLIG2", "SOX10", "NKX2-2", "MAG", "GALC", "CLDN11"),
  Microglia = c("TMEM119", "CX3CR1", "P2RY12", "CD68", "CD11b", "HLA-DRA", "ITGAM", "CSF1R", "TREM2", "IBA1"),
  Pancreatic_Beta_cells = c("INS", "GCG", "SST", "PDX1", "NKX6-1", "MAFB", "NEUROD1", "ISL1", "PCSK1", "GLP1R"),
  Myeloid_dendritic_cells = c("ITGAX", "CD1C", "CD11c", "CLEC9A", "CD40", "CD83", "HLA-DR", "CD86", "CD80", "CD207"),
  Plasmacytoid_dendritic_cells = c("IL3RA", "TCL1A", "LILRA4", "BDCA2", "BDCA4", "CLEC12A", "HLA-DR", "CD123", "TREML4", "CD2"),
  Plasma_cells = c("CD138", "MZB1", "XBP1", "CD319", "CD38", "CD27", "BLIMP1", "IRF4", "CD19", "BCMA"),
  Regulatory_T_cells = c("FOXP3", "CD25", "CTLA4", "CD4", "CD127", "GITR", "CD45RA", "CD62L", "CCR4", "CD39"),
  Th17_cells = c("RORC", "IL17A", "IL22", "CCR6", "CD161", "CD196", "IL23R", "IL26", "IL21", "CCL20"),
  Th1_cells = c("TBX21", "IFNG", "IL2", "CXCR3", "CCR5", "IL12RB2", "STAT4", "HLA-DR", "T-bet", "CD4"),
  Th2_cells = c("GATA3", "IL4", "IL5", "IL13", "CCR3", "CRTH2", "STAT6", "IL4R", "IL5RA", "CD294"),
  Memory_B_cells = c("CD27", "CD38", "IGHG3", "CD19", "CD24", "CD80", "CD95", "CD122", "CD148", "CD180"),
  Naive_B_cells = c("CD24", "CD38", "IGHD", "CD10", "CD20", "CD22", "CD23", "CD27", "CD45RA", "IGLL1"),
  CD4_positive_T_cells = c("CD4", "CD29", "LCK", "CD45RA", "CD45RO", "CD62L", "CCR7", "CXCR4", "CD154", "CD27"),
  CD8_positive_T_cells = c("CD8A", "CD8B", "PRF1", "CD27", "CD28", "CD45RA", "CD45RO", "CD57", "CD122", "CD244"),
  Germinal_center_B_cells = c("BCL6", "CD10", "CXCR5", "CD77", "CD95", "CD38", "FAS", "PECA", "CD179b", "CD179a"),
  Follicular_helper_T_cells = c("CXCR5", "PD1", "BCL6", "CD200", "SAP", "ICOS", "CD40L", "OX40", "CD84", "SLAMF6")
)
