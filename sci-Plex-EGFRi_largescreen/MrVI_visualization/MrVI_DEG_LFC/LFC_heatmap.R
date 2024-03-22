library(monocle3)
library(tidyverse)
library(ComplexHeatmap)

this.dir <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(this.dir))
rm(this.dir)

# ================================================================================
# LFC with cell and DRC annotations
# ================================================================================

heatmap_mat <- read_csv("SharedPDCL_MrVI_multivariateDE_DRC_LFC_lfcANDpvalFiltered.csv") %>%
  dplyr::rename(X1 = 1) %>%
  column_to_rownames(var = "X1")

# heatmap_mat_scaled <- scale(t(scale(t(heatmap_mat))))
heatmap_mat_scaled <- scale(heatmap_mat)
heatmap_mat_scaled[heatmap_mat_scaled > 2] <- 2
heatmap_mat_scaled[heatmap_mat_scaled < -2] <- -2
heatmap_mat_scaled[is.na(heatmap_mat_scaled)] <- 0

heatmap_mat_cluster_labels <- read_csv("SharedPDCL_MrVI_multivariateDE_DRC_HierarchicalClusters_lfcANDpvalFiltered.csv") %>%
  dplyr::rename(response_module = 1, response_module_cluster = 2) %>%
  column_to_rownames(var = "response_module")

color_scheme_rand <- readRDS("rand_color_scheme_response_module_cluster.RDS")

factor_col <- setNames(c(color_scheme_rand[1:max(heatmap_mat_cluster_labels$response_module_cluster)]),
                       heatmap_mat_cluster_labels %>%
                         distinct(response_module_cluster) %>%
                         arrange(response_module_cluster) %>%
                         pull())

col_ha_temp <- data.frame("col_names" = colnames(heatmap_mat)) %>%
  separate(col_names, into = c("Cell Line", "RM"), remove = F) %>%
  column_to_rownames(var = "col_names")

col_ha_df <- bind_cols(heatmap_mat_cluster_labels, col_ha_temp)

color <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
color_scheme_rand_rm <- sample(color, 25)
rm_col <- setNames(c(color_scheme_rand_rm),
                       1:25)

col_ha <- ComplexHeatmap::HeatmapAnnotation(df = col_ha_df,
                                            col = list(`Response Module Cluster` = factor_col,
                                                       `Cell Line` = setNames(c("darkblue", "darkred", "darkgreen"),
                                                                              c("BT112", "BT228", "BT333")),
                                                       RM = rm_col),
                                            annotation_name_gp = gpar(fontsize = 6),
                                            annotation_name_rot = 0,
                                            show_legend = c(FALSE, TRUE, FALSE),
                                            annotation_legend_param = list(
                                              labels_gp = gpar(fontsize = 5),
                                              title_gp = gpar(fontsize = 6,
                                                              fontface = "bold")
                                            ))

hmcols <- colorRampPalette(c("blue","white","red"))(50)
hm <- ComplexHeatmap::Heatmap(matrix = heatmap_mat_scaled,
                              # heatmap_mat %>% as.matrix(),
                              name = "Scaled LFC",
                              heatmap_legend_param = list(labels_gp = gpar(fontsize = 5),
                                                          title_gp = gpar(fontsize = 6,
                                                                          fontface = "bold"),
                                                          legend_height = unit(0.5, "cm"),
                                                          legend_width = unit(0.2, "cm")),
                              cluster_columns = T, cluster_rows = T, 
                              col = hmcols, 
                              # clustering_method_rows = "ward.D2",
                              # row_split = heatmap_mat_cluster_labels %>% column_to_rownames(var = "drug"),
                              column_split = heatmap_mat_cluster_labels,
                              column_gap = unit(0.2, "line"),
                              row_gap = unit(0.2, "line"),
                              # column_title = "Response Cluster",
                              column_title = NULL,
                              # row_title = NULL,
                              show_column_names = F,
                              top_annotation = col_ha,
                              # left_annotation = row_ha,
                              # rect_gp = gpar(col = "black", lwd = 0.3),
                              row_names_gp = gpar(fontsize = 5), 
                              show_row_names = F,
                              show_heatmap_legend = T)

draw(hm, legend_grouping = "original")

png("RM_LFC_filtered_cellline_annotation.png",width=5,height=5,units="in",res=1200)
draw(hm, legend_grouping = "original")
dev.off()

# ================================================================================
# LFC with cell and DRC annotations transposed
# ================================================================================

heatmap_mat <- read_csv("SharedPDCL_MrVI_multivariateDE_DRC_LFC_lfcANDpvalFiltered.csv") %>%
  dplyr::rename(X1 = 1) %>%
  column_to_rownames(var = "X1")

heatmap_mat_scaled <- t(scale(heatmap_mat))
heatmap_mat_scaled[heatmap_mat_scaled > 2] <- 2
heatmap_mat_scaled[heatmap_mat_scaled < -2] <- -2
heatmap_mat_scaled[is.na(heatmap_mat_scaled)] <- 0

heatmap_mat_cluster_labels <- read_csv("SharedPDCL_MrVI_multivariateDE_DRC_HierarchicalClusters_lfcANDpvalFiltered.csv") %>%
  dplyr::rename(response_module = 1, response_module_cluster = 2) %>%
  column_to_rownames(var = "response_module")

color_scheme_rand <- readRDS("rand_color_scheme_response_module_cluster.RDS")

factor_col <- setNames(c(color_scheme_rand[1:max(heatmap_mat_cluster_labels$response_module_cluster)]),
                       heatmap_mat_cluster_labels %>%
                         distinct(response_module_cluster) %>%
                         arrange(response_module_cluster) %>%
                         pull())

col_ha_temp <- data.frame("col_names" = colnames(heatmap_mat)) %>%
  separate(col_names, into = c("Cell Line", "RM"), remove = F) %>%
  column_to_rownames(var = "col_names")

col_ha_df <- bind_cols(heatmap_mat_cluster_labels, col_ha_temp)

color <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
color_scheme_rand_rm <- sample(color, 25)
rm_col <- setNames(c(color_scheme_rand_rm),
                   1:25)

col_ha <- ComplexHeatmap::HeatmapAnnotation(df = col_ha_df,
                                            col = list(`Response Module Cluster` = factor_col,
                                                       `Cell Line` = setNames(c("darkblue", "darkred", "darkgreen"),
                                                                              c("BT112", "BT228", "BT333")),
                                                       RM = rm_col),
                                            annotation_name_gp = gpar(fontsize = 6),
                                            annotation_name_rot = 0,
                                            show_legend = c(FALSE, TRUE, FALSE),
                                            annotation_legend_param = list(
                                              labels_gp = gpar(fontsize = 5),
                                              title_gp = gpar(fontsize = 6,
                                                              fontface = "bold")
                                            ))

hmcols <- colorRampPalette(c("blue","white","red"))(50)
hm <- ComplexHeatmap::Heatmap(matrix = heatmap_mat_scaled,
                              # heatmap_mat %>% as.matrix(),
                              name = "Scaled LFC",
                              heatmap_legend_param = list(labels_gp = gpar(fontsize = 5),
                                                          title_gp = gpar(fontsize = 6,
                                                                          fontface = "bold"),
                                                          legend_height = unit(0.5, "cm"),
                                                          legend_width = unit(0.2, "cm")),
                              cluster_columns = T, cluster_rows = T, 
                              col = hmcols, 
                              # clustering_method_rows = "ward.D2",
                              # row_split = heatmap_mat_cluster_labels %>% column_to_rownames(var = "drug"),
                              row_split = heatmap_mat_cluster_labels,
                              column_gap = unit(0.2, "line"),
                              row_gap = unit(0.2, "line"),
                              # column_title = "Response Cluster",
                              # column_title = NULL,
                              row_title = NULL,
                              show_column_names = F,
                              left_annotation = col_ha,
                              # left_annotation = row_ha,
                              # rect_gp = gpar(col = "black", lwd = 0.3),
                              row_names_gp = gpar(fontsize = 5), 
                              show_row_names = F,
                              show_heatmap_legend = T)

draw(hm, legend_grouping = "original")

png("DRC_LFC_filtered_cellline_DRC_annotation_transposed.png",width=7,height=4,units="in",res=1200)
draw(hm)
dev.off()
