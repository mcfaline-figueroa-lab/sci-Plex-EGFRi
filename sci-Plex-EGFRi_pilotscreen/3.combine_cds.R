library(tidyverse)
library(monocle3)

cds_3line <- readRDS("cds_QC_establishedlines.RDS")
cds_BT333 <- readRDS("cds_QC_BT333.RDS")

cds <- combine_cds(list(cds_3line, cds_BT333), cell_names_unique = TRUE)
rm(cds_3line, cds_BT333)

cds <- cds[,!colData(cds)$treatment %in% c("Water", "Puromycin")]

expressed_genes <- row.names(rowData(cds)[Matrix::rowSums(exprs(cds) > 0) > 
                                            dim(cds)[2]*0.05 ,])

cds <- preprocess_cds(cds, 
                      num_dim = 20, 
                      method = "PCA",
                      verbose = TRUE,
                      use_genes = expressed_genes)

cds <- monocle3::align_cds(cds,
                           alignment_group = "Cell.Line",
                           residual_model_formula_str = "~Cell.Line + log10.umi + batch")

cds <- reduce_dimension(cds, 
                        reduction_method = "UMAP", 
                        preprocess_method = "Aligned",
                        umap.min_dist = 0.05,
                        umap.n_neighbors = 15,
                        verbose = TRUE)

saveRDS(cds, "cds_pilot_screen.RDS")

