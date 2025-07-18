library(tidyverse)
library(monocle3)
this.dir <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(this.dir))
rm(this.dir)

cds <- readRDS("sci-Plex-EGFRi-siRNA-cds.RDS")

# ==============================================================================================
# UMAPs to visualize differences
# ==============================================================================================

expressed_genes <- row.names(rowData(cds)[Matrix::rowSums(exprs(cds) > 0) > 
                                            dim(cds)[2]*0.05 ,])

cds <- preprocess_cds(cds,
                      method = "PCA",
                      num_dim = 50,
                      norm_method = "log",
                      use_genes = expressed_genes)

dir.create("UMAPs")
plot_pc_variance_explained(cds) +
  theme(text = element_text(size = 6)) +
  geom_vline(xintercept = 20)
ggsave("UMAPs/PC_variance_explained_50_PCs_BT112.png", width = 2, height = 1.5, dpi = 600)

cds <- preprocess_cds(cds,
                      method = "PCA",
                      num_dim = 20,
                      use_genes = expressed_genes)

cds <- reduce_dimension(cds,
                        max_components = 2,
                        reduction_method = "UMAP",
                        umap.metric = "cosine",
                        umap.n_neighbors = 20,
                        umap.min_dist = 0.05,
                        umap.fast_sgd=FALSE,
                        cores=1,
                        verbose = T)

colData(cds)$UMAP1 <- reducedDims(cds)[["UMAP"]][,1]
colData(cds)$UMAP2 <- reducedDims(cds)[["UMAP"]][,2]

colData(cds)$siRNA_condition <- factor(colData(cds)$siRNA_condition, levels = c("siNTC", "siEGFR1", "siEGFR2"))

plot_cells(cds, color_cells_by = "siRNA_condition", cell_size = 0.3, label_cell_groups = FALSE, cell_stroke = 0.05) +
  scale_color_manual(name = "siRNA Condition", 
                     values = RColorBrewer::brewer.pal(name = "Dark2",n = 3)) +
  theme_void() +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.2,"line"), 
        legend.key.height = unit(0.6,"line")) +
  guides(color = guide_legend(title = "siRNA Condition",
                              override.aes = list(size=2)))
ggsave("UMAPs/BT112_by_siRNA_condition.png", 
       dpi = 900, width = 1.65, height = 1, bg = "white")

plot_cells(cds, color_cells_by = "timepoint", cell_size = 0.3, label_cell_groups = FALSE, cell_stroke = 0.05) +
  scale_color_manual(name = "timepoint", 
                     values = viridis::viridis(n = 4)[2:4]) +
  theme_void() +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.2,"line"), 
        legend.key.height = unit(0.6,"line")) +
  guides(color = guide_legend(title = "Timepoint",
                              override.aes = list(size=2)))
ggsave("UMAPs/BT112_by_timepoint.png", 
       dpi = 900, width = 1.5, height = 1, bg = "white")

plot_cells(cds, genes = "EGFR", cell_size = 0.3, 
           label_cell_groups = FALSE, cell_stroke = 0.05, 
           scale_to_range = F) +
  viridis::scale_color_viridis(option = "A") +
  theme_void() +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.2,"line"), 
        legend.key.height = unit(0.5,"line"),
        strip.text = element_blank()) +
  guides(color = guide_colorbar(title = "EGFR Exprs",
                              override.aes = list(size=2)))
ggsave("UMAPs/BT112_by_EGFR_exprs.png", 
       dpi = 900, width = 1.5, height = 1, bg = "white")

plot_cells(cds, color_cells_by = "siRNA_condition", cell_size = 0.3, label_cell_groups = FALSE, cell_stroke = 0.05) +
  facet_wrap(~siRNA_condition, nrow = 1) +
  scale_color_manual(name = "siRNA Condition", 
                     values = RColorBrewer::brewer.pal(name = "Dark2",n = 3)) +
  theme_void() +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.2,"line"), 
        legend.key.height = unit(0.6,"line")) +
  guides(color = guide_legend(title = "siRNA Condition",
                              override.aes = list(size=2)))
ggsave("UMAPs/BT112_by_siRNA_condition_facet.png", 
       dpi = 900, width = 2.5, height = 1, bg = "white")

plot_cells(cds, color_cells_by = "timepoint", cell_size = 0.3, label_cell_groups = FALSE, cell_stroke = 0.05) +
  facet_wrap(~timepoint, nrow = 1) +
  scale_color_manual(name = "timepoint", 
                     values = viridis::viridis(n = 4)[2:4]) +
  theme_void() +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.2,"line"), 
        legend.key.height = unit(0.6,"line")) +
  guides(color = guide_legend(title = "Timepoint",
                              override.aes = list(size=2)))
ggsave("UMAPs/BT112_by_timepoint_facet.png", 
       dpi = 900, width = 2.35, height = 1, bg = "white")

plot_cells(cds, genes = "EGFR", cell_size = 0.3, 
           label_cell_groups = FALSE, cell_stroke = 0.05, 
           scale_to_range = F) +
  facet_wrap(~timepoint) +
  viridis::scale_color_viridis(option = "A") +
  theme_void() +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.2,"line"), 
        legend.key.height = unit(0.5,"line")
        ) +
  guides(color = guide_colorbar(title = "EGFR Exprs",
                                override.aes = list(size=2)))
ggsave("UMAPs/BT112_by_EGFR_exprs_time_facet.png", 
       dpi = 900, width = 2.25, height = 1, bg = "white")


plot_cells(cds, genes = "EGFR", cell_size = 0.3, 
           label_cell_groups = FALSE, cell_stroke = 0.05, 
           scale_to_range = F) +
  facet_wrap(~siRNA_condition) +
  scale_color_manual(values = c(RColorBrewer::brewer.pal(n = 3, name = "Set1"))) +
  theme_void() +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.2,"line"), 
        legend.key.height = unit(0.5,"line")
  ) +
  guides(color = guide_legend(title = "EGFR Exprs",
                                override.aes = list(size=2)))
ggsave("UMAPs/BT112_by_EGFR_exprs_time_siRNA.png", 
       dpi = 900, width = 2.25, height = 1, bg = "white")


plot_cells(cds, color_cells_by = "timepoint", cell_size = 0.5, label_cell_groups = FALSE, cell_stroke = 0.05) +
  facet_wrap(~siRNA_condition, ncol = 1) +
  scale_color_manual(name = "Timepoint",
                     values = viridis::viridis(n = 3)) +
  theme_void() +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.2,"line"), 
        legend.key.height = unit(0.6,"line")) +
  guides(color = guide_legend(title = "Timepoint",
                              override.aes = list(size=2)))
ggsave("UMAPs/BT112_color_by_time_siRNA_condition_facet.png", 
       dpi = 900, width = 1.5, height = 2.25, bg = "white")


# Cell proportions of timepoint
temp_count <- colData(cds) %>% as.data.frame() %>%
  group_by(timepoint, siRNA_condition, replicate) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(timepoint, replicate) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup() %>%
  group_by(timepoint, siRNA_condition) %>%
  summarise(mean_prop = mean(prop))

temp_count <- colData(cds) %>% as.data.frame() %>%
  group_by(timepoint, siRNA_condition) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(prop = count / sum(count)) %>%
  ungroup()

ggplot(temp_count, aes(x = timepoint, y = prop)) +
  geom_point() +
  geom_line(aes(group = siRNA_condition)) +
  facet_wrap(~siRNA_condition) +
  monocle3:::monocle_theme_opts()

# Clustering cells
cds <- cluster_cells(cds, reduction_method = "PCA", resolution = 1e-3, verbose = T)
colData(cds)$Cluster <- clusters(cds, reduction_method = "PCA")

plot_cells(cds, color_cells_by = "Cluster", cell_size = 0.5, label_cell_groups = FALSE, cell_stroke = 0.05) +
  # scale_color_brewer(name = "siRNA Condition", palette = "Set1") +
  # facet_wrap(~siRNA_condition, ncol = 1) +
  scale_color_manual(name = "Cluster",
                     values = c("#95190C", "#610345", "#107E7D")) +
  theme_void() +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.2,"line"), 
        legend.key.height = unit(0.6,"line"),
        legend.direction = "horizontal",
        legend.position = "bottom") +
  guides(color = guide_legend(title = "Cluster",
                              override.aes = list(size=2),
                              byrow = T,
                              nrow = 1))
ggsave("UMAPs/BT112_color_by_cluster.png", 
       dpi = 900, width = 1.5, height = 1.25, bg = "white")

plot_percent_cells_positive(cds_subset = cds[rowData(cds)$gene_short_name %in% c("EGFR", "MKI67")],
                            group_cells_by = "Cluster") +
  scale_fill_manual(values = c("#95190C", "#610345", "#107E7D")) +
  theme(text = element_text(size = 6, color = "black"),
        axis.ticks = element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.5, "mm")) +
  guides(fill = "none")
ggsave("BT112_siRNA_cluster_percent_cells_positive.png",
       dpi = 600, width = 1.25, height = 1.5)

plot_percent_cells_positive(cds_subset = cds[rowData(cds)$gene_short_name %in% c("EGFR", 
                                                                                 "MKI67",
                                                                                 "BRAF",
                                                                                 "ERBB3", 
                                                                                 "TP53",
                                                                                 "KRAS",
                                                                                 "PTEN",
                                                                                 "CDK6",
                                                                                 "HLA-A",
                                                                                 "HLA-B")],
                            group_cells_by = "Cluster",
                            ncol = 2) +
  scale_fill_manual(values = c("#95190C", "#610345", "#107E7D")) +
  theme(text = element_text(size = 6, color = "black"),
        axis.ticks = element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.5, "mm")) +
  guides(fill = "none")

plot_percent_cells_positive(cds_subset = cds[rowData(cds)$gene_short_name %in% c("MKI67")],
                            group_cells_by = "Cluster",
                            ncol = 1) +
  scale_fill_manual(values = c("#95190C", "#610345", "#107E7D")) +
  theme(text = element_text(size = 6, color = "black"),
        axis.ticks = element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.5, "mm")) +
  guides(fill = "none")
ggsave("BT112_siRNA_cluster_percent_cells_positive.png",
       dpi = 600, width = 1.25, height = 1.25)

temp_count <- colData(cds) %>% as.data.frame() %>%
  # filter(timepoint != "72hr") %>%
  group_by(siRNA_condition, Cluster) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(siRNA_condition) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup()

ggplot(temp_count, aes(x = siRNA_condition, y = prop)) +
  geom_bar(stat = "identity",aes(fill = Cluster)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # facet_wrap(~timepoint) +
  monocle3:::monocle_theme_opts()


# Fisher's exact test
test_product_cluster_mat <- colData(cds) %>% as.data.frame() %>%
  select(siRNA_condition, Cluster) %>%
  group_by(siRNA_condition, Cluster) %>%
  summarise(count = n()) %>%
  pivot_wider(id_cols = siRNA_condition, 
              names_from = Cluster, 
              values_from = count, 
              values_fill = 0) %>%
  column_to_rownames(var = "siRNA_condition")
# another way to make this matrix -- different than José's method

final_cluster_enrichment_df <- data.frame()
for (time in c("24hr", "48hr", "72hr"))  {
test_product_cluster_mat <- reshape2::acast(
  colData(cds) %>% as.data.frame() %>% 
    filter(timepoint == time) %>%
    mutate(dummy = 1) %>% 
    mutate(product_dose = siRNA_condition) %>%
    select(product_dose, Cluster, dummy),
  product_dose ~ Cluster,
  value.var = "dummy",
  fun.aggregate = sum,
  fill = 0
)

ntc.counts <- test_product_cluster_mat["siNTC", ]


cluster_enrichment_df_temp <- do.call(rbind, 
                                      lapply(c("siEGFR1", "siEGFR2"),
                                 # lapply(rownames(test_product_cluster_mat), 
                                        function(product_dose){
  do.call(rbind, lapply(1:ncol(test_product_cluster_mat), function(Cluster){
    test <- fisher.test(cbind(c(test_product_cluster_mat[product_dose, Cluster], sum(
      test_product_cluster_mat[product_dose, -Cluster]
    )),
    c(ntc.counts[Cluster], sum(ntc.counts[-Cluster])))
    )
    
    data.frame(
      "product_dose" = product_dose,
      "Cluster" = Cluster,
      "odds_ratio" = unname(test$estimate),
      "p_value" = test$p.value
    )
  }))
}))

cluster_enrichment_df_temp <- cluster_enrichment_df_temp %>%
  mutate(timepoint = time) %>%
  mutate(q_value = p.adjust(p_value, method = "BH"))

final_cluster_enrichment_df <- bind_rows(final_cluster_enrichment_df, cluster_enrichment_df_temp)
}


final_cluster_enrichment_mat <- final_cluster_enrichment_df %>%
  filter(product_dose != "siNTC") %>%
  unite(product_dose, timepoint, col = "siRNA_time", sep = "_") %>%
  select(siRNA_time, Cluster, odds_ratio) %>%
  pivot_wider(id_cols = siRNA_time,
              names_from = Cluster,
              values_from = odds_ratio) %>%
  column_to_rownames("siRNA_time") %>% as.matrix()

col_ha_df <- data.frame("Cluster" = c("1", "2", "3"))
row.names(col_ha_df) <- c("1", "2", "3")
col_ha <- HeatmapAnnotation(df = col_ha_df,
                            col = list(Cluster = setNames(c("#95190C", "#610345", "#107E7D"),
                                                          c("1", "2", "3"))),
                            show_legend = F,
                            show_annotation_name = F)

final_cluster_enrichment_mat_FDR <- final_cluster_enrichment_df %>%
  filter(product_dose != "siNTC") %>%
  unite(product_dose, timepoint, col = "siRNA_time", sep = "_") %>%
  select(siRNA_time, Cluster, q_value) %>%
  pivot_wider(id_cols = siRNA_time,
              names_from = Cluster,
              values_from = q_value) %>%
  column_to_rownames("siRNA_time") %>% as.matrix()

hmcols <- circlize::colorRamp2(breaks = c(0.6, 0.8, 1, 1.2 , 1.4),
                               RColorBrewer::brewer.pal(n = 5, "RdYlBu")[c(5:1)])
hm <- ComplexHeatmap::Heatmap(final_cluster_enrichment_mat, 
                              col = hmcols,
                              name = "Odds Ratio\nto siNTC",
                              cluster_columns = FALSE,
                              cluster_rows = TRUE,
                              row_split = 3,
                              clustering_method_rows = "complete",
                              row_title = NULL,
                              column_title = "Cluster",
                              column_title_side = "bottom",
                              column_title_gp = gpar(fontsize = 8),
                              cell_fun = function(j, i, x, y, width, height, fill) {
                                if(final_cluster_enrichment_mat_FDR[i, j] < 0.05)
                                  grid.text(sprintf("*", final_cluster_enrichment_mat[i, j]),
                                            x, y - 0.1 * height, gp = gpar(fontsize = 10))},
                              bottom_annotation = col_ha,
                              row_names_gp = gpar(fontsize = 8), show_row_names = T,
                              column_names_gp = gpar(fontsize = 8),
                              column_names_rot = 0,
                              heatmap_legend_param = list(title_gp = gpar(fontsize = 8, 
                                                                          fontface = "bold"),
                                                          labels_gp = gpar(fontsize = 6),
                                                          grid_width = unit(3, "mm"),
                                                          legend_height = unit(5, "mm"),
                                                          border = TRUE),
                              rect_gp = gpar(lwd = 0.3),
                              row_dend_gp = gpar(lwd = 0.3),
                              row_dend_width = unit(0.75, "cm")
)

hm_draw <- ComplexHeatmap::draw(hm,
                                heatmap_legend_side = "right", 
                                annotation_legend_side = "right")

pdf("siRNA_cluster_fishers_exact_hm.pdf",
    height = 1.75, width = 2.75, compress = F)
hm_draw <- ComplexHeatmap::draw(hm, 
                                heatmap_legend_side = "right", 
                                annotation_legend_side = "right")
dev.off()

saveRDS(cds, 
        "sci-Plex-EGFRi-siRNA-cds-clustered.RDS")


# ================================================================================
# sci-Plex-GxE adaptive resistance signature
# ================================================================================
calculate_aggregate_expression_score <- function(cds, signature_genes, from_id = FALSE){
  
  if(from_id == TRUE){
    cds_subset = cds[rowData(cds)$id %in% signature_genes,]
  }
  else(cds_subset = cds[rowData(cds)$gene_short_name %in% signature_genes,])
  aggregate_signature_expression = exprs(cds_subset)
  aggregate_signature_expression = t(t(aggregate_signature_expression) / pData(cds_subset)$Size_Factor)
  aggregate_signature_expression = Matrix::colSums(aggregate_signature_expression)
  signature_score = log(aggregate_signature_expression+1)
  return(aggregate_signature_expression)
}


adaptive_up_genes <- read_csv("adaptive_program.csv",
                              col_names = TRUE) %>%
  filter(grepl(pattern = "pregulated", gene_module))

colData(cds)$adaptive_up <- calculate_aggregate_expression_score(cds, 
                                                                 signature_genes = adaptive_up_genes$gene_short_name,
                                                                 from_id = FALSE)


adaptive_plot <- colData(cds) %>% as.data.frame() %>%
  # group_by(cell_line) %>%
  mutate(adaptive_up_scaled = scale(adaptive_up) %>% as.vector())


violin_data_mean <- adaptive_plot %>%
  group_by(Cluster) %>%
  summarise(mean_expression = mean(adaptive_up), 
            median_expression = median(adaptive_up), .groups = "drop")

ggplot(data = adaptive_plot,
       aes(x = Cluster, y = adaptive_up)) +
  geom_violin(aes(fill = Cluster),
              linewidth = 0.4) +
  geom_point(data = violin_data_mean,
             aes(x = Cluster, y = mean_expression),
             size = 0.25,
             color = ifelse(violin_data_mean$Cluster == 2, "white", "black")) +
  ylab("Adaptive Score") +
  scale_fill_manual(values = c("#95190C", "#610345", "#107E7D")) +
  theme(text = element_text(size = 7),
        # axis.text.x = element_blank(),
        axis.ticks.length.y = unit(0.75,"mm"),
        axis.ticks = element_line(linewidth = 0.2),
        # axis.title = element_blank(),
        # axis.ticks.x = element_blank(),
        legend.position = "bottom",
        legend.key.height = unit(4, "mm"),
        legend.key.width = unit(2, "mm"),
        legend.margin = margin(t = -6, l = -5),
        legend.title = element_text(size = 6)
  ) +
  # guides(color = guide_legend(title = "siRNA", nrow = 1,
  #                             override.aes = list(shape = 16, size = 2))) +
  guides(fill = "none") +
  monocle3:::monocle_theme_opts()

ggsave("siRNA_adaptive_by_cluster.png", dpi = 600, height = 1.5, width = 1.5)

diff_test_clean_final <- data.frame()
for (clust in c("1", "2", "3")) {
  score_df_subset <- colData(cds) %>% as.data.frame() %>%
    mutate(cluster_compare = case_when(
      Cluster == clust ~ paste0("Cluster",clust),
      TRUE ~ "Other"
    )) %>%
    mutate(cluster_compare = factor(cluster_compare, levels = c("Other", paste0("Cluster",clust))))
  
  model_to_test <- "adaptive_up ~ cluster_compare"
  # model_to_test <- "adaptive_up ~ cluster_compare"
  
  diff_test_results <- speedglm::speedglm(data = score_df_subset, formula = model_to_test)
  diff_test_results <- broom::tidy(diff_test_results)
  
  diff_test_clean <- data.frame(Cluster = paste0("Cluster",clust),
                                term = diff_test_results$term, 
                                beta = diff_test_results$estimate, 
                                p_value = diff_test_results$p.value)
  
  diff_test_clean_final <- bind_rows(diff_test_clean_final, diff_test_clean)
}

diff_test_clean_final$q_value <- p.adjust(diff_test_clean_final$p_value, "BH")

# ================================================================================
# APM MHC-CI signature
# ================================================================================
APM_MHC_CI_genes <- c("PSMB5",
                      "PSMB6",
                      "PSMB7",
                      "PSMB8",
                      "PSMB9",
                      "PSMB10",
                      "TAP1",
                      "TAP2",
                      "ERAP1",
                      "ERAP2",
                      "CANX",
                      "CALR",
                      "PDIA3",
                      "TAPBP",
                      "B2M",
                      "HLA-A",
                      "HLA-B",
                      "HLA-C")


colData(cds)$APM_MHC_CI_score <- calculate_aggregate_expression_score(cds, signature_genes = APM_MHC_CI_genes)

MHC_plot <- colData(cds) %>% as.data.frame() %>%
  mutate(MHC_up = APM_MHC_CI_score) %>%
  mutate(MHC_up_scaled = scale(APM_MHC_CI_score) %>% as.vector())

violin_data_mean <- MHC_plot %>%
  group_by(Cluster) %>%
  summarise(mean_expression = mean(MHC_up), 
            median_expression = median(MHC_up), .groups = "drop")

ggplot(data = MHC_plot,
       aes(x = Cluster, y = MHC_up)) +
  geom_violin(aes(fill = Cluster),
              linewidth = 0.4) +
  geom_point(data = violin_data_mean,
             aes(x = Cluster, y = mean_expression),
             size = 0.25,
             color = ifelse(violin_data_mean$Cluster == 2, "white", "black")) +
  ylab("APM MHC-CI Score") +
  scale_fill_manual(values = c("#95190C", "#610345", "#107E7D")) +
  theme(text = element_text(size = 7),
        # axis.text.x = element_blank(),
        axis.ticks.length.y = unit(0.75,"mm"),
        axis.ticks = element_line(linewidth = 0.2),
        # axis.title = element_blank(),
        # axis.ticks.x = element_blank(),
        legend.position = "bottom",
        legend.key.height = unit(4, "mm"),
        legend.key.width = unit(2, "mm"),
        legend.margin = margin(t = -6, l = -5),
        legend.title = element_text(size = 6)
  ) +
  # guides(color = guide_legend(title = "siRNA", nrow = 1,
  #                             override.aes = list(shape = 16, size = 2))) +
  guides(fill = "none") +
  monocle3:::monocle_theme_opts()

ggsave("siRNA_MHC-I_by_cluster.png", dpi = 600, height = 1.5, width = 1.5)


