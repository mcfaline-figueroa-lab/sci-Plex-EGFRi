suppressPackageStartupMessages({
  library(ggplot2)
  library(plyr)
  library(dplyr)
  library(tidyr)
  library(devtools)
  library(monocle3)
})

# Download sci-Plex repo to load the following function and data frame
suppressPackageStartupMessages({
  source("~/Documents/github_repos/sci-plex/bin/cell_cycle.R")
  cc.genes = readRDS("~/Documents/github_repos/sci-plex/bin/cc.genes.RDS")
})

cds.list <- readRDS("cds.list.rds")

for(PDCL in names(cds.list)){
  cds.list[[PDCL]] <- estimate_size_factors(cds.list[[PDCL]])
  cds.list[[PDCL]] <- detect_genes(cds.list[[PDCL]])
  cds.list[[PDCL]] <- estimate_cell_cycle(cds.list[[PDCL]], g1s_markers = cc.genes$s.genes, g2m_markers = cc.genes$g2m.genes)
  colData(cds.list[[PDCL]])$log10_umi <- log10(colData(cds.list[[PDCL]])$n.umi)
  colData(cds.list[[PDCL]])$Cluster <- NULL
  }

expressed_genes.list <- list()

for(PDCL in names(cds.list)){
  expressed_genes.list[[PDCL]] <- row.names(rowData(cds.list[[PDCL]])[Matrix::rowSums(exprs(cds.list[[PDCL]]) > 0) > 
                                            dim(cds.list[[PDCL]])[2]*0.01 ,])
  print(length(expressed_genes.list[[PDCL]]))
}

for(PDCL in names(cds.list)){
  cds.list[[PDCL]] <- preprocess_cds(cds.list[[PDCL]],
                        method = "PCA",
                        num_dim = 50,
                        norm_method = "log",
                        use_genes = expressed_genes.list[[PDCL]])
}

dir.create("UMAPs")
plot_pc_variance_explained(cds.list[["BT112"]]) +
  theme(text = element_text(size = 6)) +
  geom_vline(xintercept = 20)
ggsave("UMAPs/PC_variance_explained_50_PCs_BT112.png", width = 2, height = 1.5, dpi = 600)

plot_pc_variance_explained(cds.list[["BT228"]]) +
  theme(text = element_text(size = 6)) +
  geom_vline(xintercept = 20)
ggsave("UMAPs/PC_variance_explained_50_PCs_BT228.png", width = 2, height = 1.5, dpi = 600)

plot_pc_variance_explained(cds.list[["BT333"]]) +
  theme(text = element_text(size = 6)) +
  geom_vline(xintercept = 20)
ggsave("UMAPs/PC_variance_explained_50_PCs_BT333.png", width = 2, height = 1.5, dpi = 600)

for(PDCL in names(cds.list)){
  cds.list[[PDCL]] <- preprocess_cds(cds.list[[PDCL]],
                        method = "PCA",
                        num_dim = 20,
                        norm_method = "log",
                        use_genes = expressed_genes.list[[PDCL]])
    
  cds.list[[PDCL]] <- reduce_dimension(cds.list[[PDCL]],
                          max_components = 2,
                          reduction_method = "UMAP",
                          umap.metric = "cosine",
                          umap.n_neighbors = 15,
                          umap.min_dist = 0.1,
                          umap.fast_sgd=FALSE,
                          cores=1,
                          verbose = T)
  
  colData(cds.list[[PDCL]])$UMAP1 <- reducedDims(cds.list[[PDCL]])[["UMAP"]][,1]
  colData(cds.list[[PDCL]])$UMAP2 <- reducedDims(cds.list[[PDCL]])[["UMAP"]][,2]
}

for(PDCL in names(cds.list)){
  
  cds.list[[PDCL]] <- reduce_dimension(cds.list[[PDCL]],
                                       max_components = 2,
                                       preprocess_method = "PCA",
                                       reduction_method = "UMAP",
                                       umap.metric = "cosine",
                                       umap.n_neighbors = 15,
                                       umap.min_dist = 0.1,
                                       umap.fast_sgd=FALSE,
                                       cores=1,
                                       verbose = T)
  
  colData(cds.list[[PDCL]])$UMAP1 <- reducedDims(cds.list[[PDCL]])[["UMAP"]][,1]
  colData(cds.list[[PDCL]])$UMAP2 <- reducedDims(cds.list[[PDCL]])[["UMAP"]][,2]
  
}

for(PDCL in names(cds.list)){
  cds.list[[PDCL]] <- cluster_cells(cds.list[[PDCL]],
                                    reduction_method = "PCA",
                                    resolution = 2e-3,
                                    num_iter = 1,
                                    random_seed = 2016L)
  colData(cds.list[[PDCL]])$Cluster <- clusters(cds.list[[PDCL]], 
                                                reduction_method = "PCA")

}

cluster_diff_test.list <- list()

for(PDCL in c("BT228","BT333")){
  
  cluster_diff_test.list[[PDCL]] <- list()
  
  for(cluster in unique(colData(cds.list[[PDCL]])$Cluster)){
    
    cds_subset <- cds.list[[PDCL]]
    colData(cds_subset)$Cluster_id <- sapply(colData(cds_subset)$Cluster, function(x){ifelse(x == cluster, paste0("_",cluster), "out")})
    colData(cds_subset)$Cluster_id <- factor(colData(cds_subset)$Cluster_id, levels = c("out", paste0("_",cluster)))
    
    cluster_diff_test.list[[PDCL]][[cluster]] <- fit_models(cds_subset[expressed_genes.list[[PDCL]]],
                                                            model_formula_str = "~Cluster_id",
                                                            cores = 1)
    
    cluster_diff_test.list[[PDCL]][[cluster]] <- coefficient_table(cluster_diff_test.list[[PDCL]][[cluster]]) %>% 
      dplyr::select(-model, -model_summary)
    
    cluster_diff_test.list[[PDCL]][[cluster]]$cell_line <- rep(PDCL, nrow(cluster_diff_test.list[[PDCL]][[cluster]]))
    cluster_diff_test.list[[PDCL]][[cluster]]$Cluster <- rep(cluster, nrow(cluster_diff_test.list[[PDCL]][[cluster]]))
    
    rm(cds_subset)
    
  }
  
}

for(PDCL in c("BT228","BT333")){
  
  cluster_diff_test.list[[PDCL]] <- do.call("rbind",cluster_diff_test.list[[PDCL]])
  cluster_diff_test.list[[PDCL]]$q_value <- p.adjust(cluster_diff_test.list[[PDCL]]$p_value, method = "BH")
  
}

cluster_diff_test_results <- do.call("rbind", cluster_diff_test.list)

# saveRDS(cluster_diff_test_results, "cluster_diff_test_results.rds")

plot_cells(cds.list[["BT112"]], color_cells_by = "Cluster", cell_size = 0.5, label_cell_groups = FALSE, cell_stroke = 0.05) +
  scale_color_brewer("PCA\ncluster",palette = "Set1") +
  theme_void() +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.2,"line"), 
        legend.key.height = unit(0.2,"line")) +
  guides(guides(fill = guide_legend(override.aes = list(size=1))))
ggsave("UMAPs/BT112_by_cluster.png", dpi = 900, width = 1, height = 1)

plot_cells(cds.list[["BT228"]], color_cells_by = "Cluster", cell_size = 0.5, label_cell_groups = FALSE, cell_stroke = 0.05) +
  scale_color_brewer("PCA\ncluster",palette = "Set1") +
  theme_void() +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.2,"line"), 
        legend.key.height = unit(0.2,"line")) +
  guides(guides(fill = guide_legend(override.aes = list(size=1))))
ggsave("UMAPs/BT228_by_cluster.png", dpi = 900, width = 1, height = 1)

plot_cells(cds.list[["BT333"]], color_cells_by = "Cluster", cell_size = 0.5, label_cell_groups = FALSE, cell_stroke = 0.05) +
  scale_color_brewer("PCA\ncluster",palette = "Set1") +
  theme_void() +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.2,"line"), 
        legend.key.height = unit(0.2,"line")) +
  guides(guides(fill = guide_legend(override.aes = list(size=1))))
ggsave("UMAPs/BT333_by_cluster.png", dpi = 900, width = 1, height = 1)

plot_cells(cds.list[["BT112"]], color_cells_by = "proliferation_index", cell_size = 0.5, label_cell_groups = FALSE, cell_stroke = 0.05) +
  viridis::scale_color_viridis("Proliferation\nindex") +
  theme_void() +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.3,"line"), 
        legend.key.height = unit(0.3,"line")) +
  guides(guides(fill = guide_legend(override.aes = list(size=1))))
ggsave("UMAPs/BT112_by_proliferation_index.png", dpi = 900, width = 1.2, height = 1)

plot_cells(cds.list[["BT228"]], color_cells_by = "proliferation_index", cell_size = 0.5, label_cell_groups = FALSE, cell_stroke = 0.05) +
  viridis::scale_color_viridis("Proliferation\nindex") +
  theme_void() +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.3,"line"), 
        legend.key.height = unit(0.3,"line")) +
  guides(guides(fill = guide_legend(override.aes = list(size=1))))
ggsave("UMAPs/BT228_by_proliferation_index.png", dpi = 900, width = 1.2, height = 1)

plot_cells(cds.list[["BT333"]], color_cells_by = "proliferation_index", cell_size = 0.5, label_cell_groups = FALSE, cell_stroke = 0.05) +
  viridis::scale_color_viridis("Proliferation\nindex") +
  theme_void() +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.3,"line"), 
        legend.key.height = unit(0.3,"line")) +
  guides(guides(fill = guide_legend(override.aes = list(size=1))))
ggsave("UMAPs/BT333_by_proliferation_index.png", dpi = 900, width = 1.2, height = 1)

plot_cells(cds.list[["BT112"]], color_cells_by = "log10_umi", cell_size = 0.5, label_cell_groups = FALSE, cell_stroke = 0.05) +
  viridis::scale_color_viridis("log10(UMIs)", option = "magma") +
  theme_void() +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.3,"line"), 
        legend.key.height = unit(0.3,"line")) +
  guides(guides(fill = guide_legend(override.aes = list(size=1))))
ggsave("UMAPs/BT112_by_log_umi.png", dpi = 900, width = 1.2, height = 1)

plot_cells(cds.list[["BT228"]], color_cells_by = "log10_umi", cell_size = 0.5, label_cell_groups = FALSE, cell_stroke = 0.05) +
  viridis::scale_color_viridis("log10(UMIs)", option = "magma") +
  theme_void() +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.3,"line"), 
        legend.key.height = unit(0.3,"line")) +
  guides(guides(fill = guide_legend(override.aes = list(size=1))))
ggsave("UMAPs/BT228_by_log_umi.png", dpi = 900, width = 1.2, height = 1)

plot_cells(cds.list[["BT333"]], color_cells_by = "log10_umi", cell_size = 0.5, label_cell_groups = FALSE, cell_stroke = 0.05) +
  viridis::scale_color_viridis("log10(UMIs)", option = "magma") +
  theme_void() +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.3,"line"), 
        legend.key.height = unit(0.3,"line")) +
  guides(guides(fill = guide_legend(override.aes = list(size=1))))
ggsave("UMAPs/BT333_by_log_umi.png", dpi = 900, width = 1.2, height = 1)

plot_cells(cds.list[["BT112"]], genes = c("FGFR1"), cell_size = 0.5, label_cell_groups = FALSE, cell_stroke = 0.05) +
  viridis::scale_color_viridis("Expression", option = "magma") +
  theme_void() +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.3,"line"), 
        legend.key.height = unit(0.3,"line")) +
  guides(guides(fill = guide_legend(override.aes = list(size=1))))
ggsave("UMAPs/BT112_by_FGFR1.png", dpi = 900, width = 1.2, height = 1)

plot_cells(cds.list[["BT228"]], genes = c("FGFR1"), cell_size = 0.5, label_cell_groups = FALSE, cell_stroke = 0.05) +
  viridis::scale_color_viridis("Expression", option = "magma") +
  theme_void() +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.3,"line"), 
        legend.key.height = unit(0.3,"line")) +
  guides(guides(fill = guide_legend(override.aes = list(size=1))))
ggsave("UMAPs/BT228_by_FGFR1.png", dpi = 900, width = 1.2, height = 1)

plot_cells(cds.list[["BT333"]], genes = c("FGFR1"), cell_size = 0.5, label_cell_groups = FALSE, cell_stroke = 0.05) +
  viridis::scale_color_viridis("Expression", option = "magma") +
  theme_void() +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.3,"line"), 
        legend.key.height = unit(0.3,"line")) +
  guides(guides(fill = guide_legend(override.aes = list(size=1))))
ggsave("UMAPs/BT333_by_FGFR1.png", dpi = 900, width = 1.2, height = 1)

plot_cells(cds.list[["BT112"]], genes = c("BRAF"), cell_size = 0.5, label_cell_groups = FALSE, cell_stroke = 0.05) +
  viridis::scale_color_viridis("Expression", option = "magma") +
  theme_void() +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.3,"line"), 
        legend.key.height = unit(0.3,"line")) +
  guides(guides(fill = guide_legend(override.aes = list(size=1))))
ggsave("UMAPs/BT112_by_BRAF.png", dpi = 900, width = 1.2, height = 1)

plot_cells(cds.list[["BT228"]], genes = c("BRAF"), cell_size = 0.5, label_cell_groups = FALSE, cell_stroke = 0.05) +
  viridis::scale_color_viridis("Expression", option = "magma") +
  theme_void() +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.3,"line"), 
        legend.key.height = unit(0.3,"line")) +
  guides(guides(fill = guide_legend(override.aes = list(size=1))))
ggsave("UMAPs/BT228_by_BRAF.png", dpi = 900, width = 1.2, height = 1)

plot_cells(cds.list[["BT333"]], genes = c("BRAF"), cell_size = 0.5, label_cell_groups = FALSE, cell_stroke = 0.05) +
  viridis::scale_color_viridis("Expression", option = "magma") +
  theme_void() +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.3,"line"), 
        legend.key.height = unit(0.3,"line")) +
  guides(guides(fill = guide_legend(override.aes = list(size=1))))
ggsave("UMAPs/BT333_by_BRAF.png", dpi = 900, width = 1.2, height = 1)

# Plot RTK pathway expression across all lines

cds <- combine_cds(cds.list)

plot_percent_cells_positive(cds[rowData(cds)$gene_short_name %in% c("EGFR","ERBB2","ERBB3","ERBB4","PDGFRA","PDGFRB","MET","FGFR1","FGFR2","FGFR3","FGFR4"),],
                            group_cells_by = "Cell.Line", ncol = 11,
                            panel_order = c("EGFR","ERBB2","ERBB3","ERBB4","PDGFRA","PDGFRB","MET","FGFR1","FGFR2","FGFR3","FGFR4")) +
  scale_fill_manual(values = c("BT112" = "navy","BT228" = "brown4","BT333" = "darkgreen")) +
  scale_color_manual(values = c("BT112" = "black","BT228" = "black","BT333" = "black")) +
  theme(text = element_text(size = 6),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  xlab("Patient-derived GBM model")
ggsave("Marker_plots/TK_levels_across_PDCLs.png", dpi = 900, width = 6.5, height = 1.5)

plot_percent_cells_positive(cds[rowData(cds)$gene_short_name %in% c("PTEN","KRAS","NRAS","HRAS","NF1","NF2","ARAF","BRAF","RAF1","PIK3CA","PIK3CB"),],
                            group_cells_by = "Cell.Line", ncol = 11,
                            panel_order = c("PTEN","KRAS","NRAS","HRAS","NF1","NF2","ARAF","BRAF","RAF1","PIK3CA","PIK3CB")) +
  scale_fill_manual(values = c("BT112" = "navy","BT228" = "brown4","BT333" = "darkgreen")) +
  scale_color_manual(values = c("BT112" = "black","BT228" = "black","BT333" = "black")) +
  theme(text = element_text(size = 6),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  xlab("Patient-derived GBM model")
ggsave("Marker_plots/RTK_pathway_levels_across_PDCLs.png", dpi = 900, width = 6.5, height = 1.5)

plot_percent_cells_positive(cds[rowData(cds)$gene_short_name %in% c("MDM2","MDM4","TP53","TP63","TP73","CDKN2A","CDKN2C","CDK4","CDK6","RB1","MGMT"),],
                            group_cells_by = "Cell.Line", ncol = 11,
                            panel_order = c("MDM2","MDM4","TP53","TP63","TP73","CDKN2A","CDKN2C","CDK4","CDK6","RB1","MGMT")) +
  scale_fill_manual(values = c("BT112" = "navy","BT228" = "brown4","BT333" = "darkgreen")) +
  scale_color_manual(values = c("BT112" = "black","BT228" = "black","BT333" = "black")) +
  theme(text = element_text(size = 6),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  xlab("Patient-derived GBM model")
ggsave("Marker_plots/p53_RB_pathway_levels_across_PDCLs.png", dpi = 900, width = 6.5, height = 1.5)



















