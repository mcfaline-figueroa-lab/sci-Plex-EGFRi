library(tidyverse)
library(monocle3)
this.dir <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(this.dir))
rm(this.dir)

source("cell_cycle.R")
cc.genes = readRDS("cc.genes.RDS")

cds <- readRDS("sci-Plex-EGFRi-GFfree-cds-combined.RDS")

# Filtering to just DMSO treated cells
cds <- cds[,colData(cds)$drug == "DMSO"]
colData(cds)$combined <- paste0(colData(cds)$cell_line,"_",colData(cds)$GF_condition)

# ==============================================================================================
# UMAPs to visualize differences
# ==============================================================================================

cds.list <- list()
for(PDCL in c("BT112", "BT333")){
  cds.list[[PDCL]] <- cds[,colData(cds)$cell_line == PDCL]
  cds.list[[PDCL]] <- estimate_size_factors(cds.list[[PDCL]])
  cds.list[[PDCL]] <- detect_genes(cds.list[[PDCL]])
  colData(cds.list[[PDCL]])$log10_umi <- log10(colData(cds.list[[PDCL]])$n_umi)
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

plot_pc_variance_explained(cds.list[["BT333"]]) +
  theme(text = element_text(size = 6)) +
  geom_vline(xintercept = 20)
ggsave("UMAPs/PC_variance_explained_50_PCs_BT333.png", width = 2, height = 1.5, dpi = 600)


for(PDCL in c("BT112", "BT333")){
  cds.list[[PDCL]] <- preprocess_cds(cds.list[[PDCL]],
                                     method = "PCA",
                                     num_dim = 15,
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

plot_cells(cds.list[["BT112"]], color_cells_by = "GF_condition", cell_size = 0.6, label_cell_groups = FALSE, cell_stroke = 0.05) +
  # scale_color_brewer(name = "GF Condition", palette = "Set1") +
  scale_color_manual(name = "GF Condition", 
                     values = RColorBrewer::brewer.pal(name = "Set1",n = 5)[3:5]) +
  theme_void() +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.2,"line"), 
        legend.key.height = unit(0.2,"line")) +
  guides(fill = guide_legend(title = "GF Condition",override.aes = list(size=1)))
ggsave("UMAPs/BT112_by_GF_condition.png", dpi = 900, width = 1.5, height = 1, bg = "white")

plot_cells(cds.list[["BT333"]], color_cells_by = "GF_condition", cell_size = 0.6, label_cell_groups = FALSE, cell_stroke = 0.05) +
  scale_color_manual(name = "GF Condition", 
                     values = RColorBrewer::brewer.pal(name = "Set1",n = 5)[3:5]) +
  theme_void() +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.2,"line"), 
        legend.key.height = unit(0.2,"line")) +
  guides(fill = guide_legend(title = "GF Condition",override.aes = list(size=1)))
ggsave("UMAPs/BT333_by_GF_condition.png", dpi = 900, width = 1.5, height = 1, bg = "white")

# ==============================================================================================
# Clustering DMSO-only GF conditions
# ==============================================================================================

for (PDCL in c("BT112", "BT333")) {
  colData(cds.list[[PDCL]])$cell <- colData(cds.list[[PDCL]])$cell_ID
  
  cds.list[[PDCL]] <- cluster_cells(cds.list[[PDCL]],
                                    reduction_method = "PCA",
                                    resolution = 5e-3)
  message("# Clusters: ", length(unique(clusters(cds.list[[PDCL]], reduction_method = "PCA"))))
  
  colData(cds.list[[PDCL]])$Cluster <- clusters(cds.list[[PDCL]], reduction_method = "PCA")
  
  col_clust <- case_when(
    PDCL == "BT112" ~ RColorBrewer::brewer.pal(name = "Set3",n = 5)[3:4],
    PDCL == "BT333" ~ RColorBrewer::brewer.pal(name = "Set3",n = 5)[c(2,5)]
  )
  plot_cells(cds.list[[PDCL]], color_cells_by = "Cluster", cell_size = 0.6, 
             label_cell_groups = FALSE, cell_stroke = 0.05) +
    # scale_color_brewer(name = "GF Condition", palette = "Set1") +
    scale_color_manual(name = "GF Condition", 
                       values = col_clust) +
    theme_void() +
    theme(text = element_text(size = 6),
          legend.key.width = unit(0.2,"line"), 
          legend.key.height = unit(0.2,"line")) +
    guides(fill = guide_legend(title = "Cluster",override.aes = list(size=1)))
  
  ggsave(paste0("UMAPs/", PDCL,"_by_cluster.png"), 
         dpi = 900, width = 1.35, height = 1, bg = "white")
  
  cluster_df_for_plot <- colData(cds.list[[PDCL]]) %>% as.data.frame() %>%
    group_by(GF_condition, Cluster) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(GF_condition) %>%
    mutate(prop = count / sum(count))
  
  ggplot(cluster_df_for_plot, aes(x = GF_condition, y = prop, fill = as.factor(Cluster))) +
    geom_bar(stat = "identity", position = "stack", width = 0.7) +
    theme(text = element_text(size = 6),
          legend.key.width = unit(0.3,"line"), 
          legend.key.height = unit(0.2,"line"),
          legend.position = "bottom",
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.ticks = element_line(linewidth = unit(0.35, "mm")),
          axis.ticks.length = unit(0.75, "mm"),
          legend.margin = margin(t = -5)
          ) +
    scale_fill_manual(name = "Cluster",
                       values = col_clust) +
    ylab("Proportion") +
    guides(fill = guide_legend(byrow = T, nrow = 1, override.aes = list(size = 1))) +
    monocle3:::monocle_theme_opts()
  ggsave(paste0("UMAPs/", PDCL,"_cluster_proportion_bar_plot.png"), 
         dpi = 900, width = 1.25, height = 1.5, bg = "white")
  
}

saveRDS(cds.list, "sci-Plex-EGFRi-GFfree-clustered-cds-list.RDS")

# ==============================================================================================
# Running DEG test to compare clusters
# ==============================================================================================

treatment_diff_test.list <- list()
expressed_genes.list <- list()
for (PDCL in c("BT112", "BT333")) {
  expressed_genes.list[[PDCL]] <- rowData(cds.list[[PDCL]]) %>% as.data.frame() %>%
    filter(num_cells_expressed >= nrow(colData(cds.list[[PDCL]])) * 0.05) %>%
    select(id) %>%
    pull()
  
  message("Running DEG analysis for ", PDCL, " using ", 
          length(expressed_genes.list[[PDCL]])," expressed genes")
  
  colData(cds.list[[PDCL]])$Cluster_factor <- factor(paste0("Cluster",
                                                            colData(cds.list[[PDCL]])$Cluster),
                                              levels = c("Cluster1","Cluster2"))
  
  treatment_diff_test.list[[PDCL]] <- fit_models(cds.list[[PDCL]][expressed_genes.list[[PDCL]],],
                                                                           model_formula_str = "~ Cluster_factor",
                                                                           cores = 1)
  treatment_diff_test.list[[PDCL]] <- coefficient_table(treatment_diff_test.list[[PDCL]]) %>% 
    dplyr::select(-model,-model_summary)
  
  treatment_diff_test.list[[PDCL]]$cell_line <- PDCL
  
  message("Finished ",PDCL)
  
}
saveRDS(treatment_diff_test.list, "UMAPs/cluster_treatment_diff_test.list.RDS")
treatment_diff_test.list <- readRDS("UMAPs/cluster_treatment_diff_test.list.RDS")

treatment_diff_test_results <- do.call("rbind",treatment_diff_test.list)

# For BT112: Positive effect represents genes upregulated for No GF
# For BT333: Negative effect represents genes upregulated for No GF
sig_genes_results <- treatment_diff_test_results %>%
  filter(grepl("Cluster", term),
         q_value < 0.01, 
         abs(normalized_effect) > 0.05)
  
# ==============================================================================================
# GSEA for each PDCL
# ==============================================================================================

library(fgsea)
# https://stephenturner.github.io/deseq-to-fgsea/
# https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html

pathways.oncogenic <- gmtPathways("c6.all.v2023.2.Hs.symbols.gmt")

# ranks must look like a named vector of values
##         A1BG     A1BG-AS1         A1CF          A2M      A2M-AS1 
##  0.679946437 -1.793291412 -0.126192849 -1.259539478  0.875346116 

library(gridExtra)
fgseaRes_list <- list()
for (PDCL in c("BT112", "BT333")) {
  ranks <- treatment_diff_test_results %>%
    filter(cell_line == PDCL) %>%
    filter(grepl(pattern = "Cluster", term)) %>%
    arrange(desc(normalized_effect)) %>%
    select(gene_short_name, normalized_effect) %>%
    distinct(gene_short_name, .keep_all = T)
  
  rank_vector <- ranks$normalized_effect
  names(rank_vector) <- ranks$gene_short_name
  
  fgseaRes <- fgsea(pathways=pathways.oncogenic, stats=rank_vector)
  fgseaRes_list[[PDCL]] <- fgseaRes %>% 
    mutate(cell_line = PDCL)
  
}

fgseaRes_final <- do.call("rbind",fgseaRes_list)

fgseaRes_final <- fgseaRes_final %>%
  mutate(q_value = padj) %>%
  mutate(signif = ifelse(q_value <= 0.05, "yes", "no"))

fgseaRes_plot <- fgseaRes_final %>%
  filter(q_value <= 0.01) %>%
  filter(pathway %in% c("RPS14_DN.V1_DN",
                        "RB_P130_DN.V1_UP",
                        "RB_P107_DN.V1_UP",
                        "VEGF_A_UP.V1_DN",
                        "E2F1_UP.V1_UP",
                        "MTOR_UP.V1_UP",
                        "MTOR_UP.V1_DN",
                        "ERBB2_UP.V1_DN",
                        "ERBB2_UP.V1_UP",
                        "VEGF_A_UP.V1_UP"))

ggplot(fgseaRes_final %>%
         mutate(NES = case_when(
           cell_line == "BT333" ~ -NES, # flipping so relative to No GF condition
           TRUE ~ NES
         )) %>%
         filter(pathway %in% fgseaRes_plot$pathway) %>%
         mutate(pathway = fct_reorder(.f = pathway, .x = NES, .fun = "mean")) %>%
         arrange(pathway) %>%
         mutate(plot_signif = case_when(
           padj <= 0.05 ~ "*",
           TRUE ~ NA
         )),
       aes(x = NES, y = pathway, fill = cell_line)) +
  geom_bar(stat = "identity", position = "dodge2", aes(fill = cell_line), width = 0.9,
           color = "black", linewidth = 0.1) +
  geom_text(aes(label = plot_signif, group = cell_line),
            position = position_dodge(width = 0.9),
            hjust = 0) +
  scale_fill_manual(name = "PDCL",
                    labels = c("BT112 Cluster 2",
                               "BT333 Cluster 1"),
                    values = c(RColorBrewer::brewer.pal(name = "Set3",n = 5)[4],
                               RColorBrewer::brewer.pal(name = "Set3",n = 5)[2])) +
  xlab("NES") +
  theme(axis.title.y = element_blank(),
        text = element_text(size = 6),
        legend.position = "bottom",
        legend.key.size = unit(3, "mm"),
        legend.margin = margin(t = -7, l = -100),
        axis.ticks = element_line(linewidth = 0.1),
        axis.ticks.length = unit(0.6, "mm")) +
  guides(fill = guide_legend(byrow = T, nrow = 1)) +
  monocle3:::monocle_theme_opts()

ggsave("UMAPs/cluster_GSEA_oncogenic_bar_plot.png",
       dpi = 600, height = 2.5, width = 2.75)

