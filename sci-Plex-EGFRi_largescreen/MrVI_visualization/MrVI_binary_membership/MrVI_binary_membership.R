library(monocle3)
library(tidyverse)
library(ComplexHeatmap)

this.dir <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(this.dir))
rm(this.dir)

# ================================================================================
# binary heatmap membership - transposed (sideways) 
# and trying a lower resolution cut for larger drug groups
# ================================================================================

heatmap_mat <- read_csv("SharedPDCL_MrVI_Drug_DRC_BinaryMembership_Matrix_lfcANDpvalFiltered.csv") %>%
  dplyr::rename(X1 = 1) %>%
  column_to_rownames(var = "X1") %>%
  t()

drug_annotations <- read_csv("Annotations_final_RG.csv") %>%
  select(drug, Class = class_broad, Reversible = reversible)

drug_annotation_class <- drug_annotations %>%
  arrange(Class, .locale = "en") %>%
  distinct(Class) %>%
  bind_rows(data.frame("Class" = "control"))

drug_annotations_control <- data.frame("drug" = c("DMSO", "Media", "PBS"),
                                       "Reversible" = c(NA,NA,NA),
                                       "Class" = c("control", "control", "control"))

row_ha_df <- bind_rows(drug_annotations, drug_annotations_control) %>%
  column_to_rownames(var = "drug")

row_ha_df <- row_ha_df[colnames(heatmap_mat),]

class_col <- setNames(c("#DC2127", "#236CAC", "#3DA448", "#863792", "#F16924",
                        "#F3EB30", "#974421", "#EA69A7", "#878888", "#53B993", "#F47751",
                        "#7A8EC1", "#D973AD", "#99CA48", "#FCD520", "#DFB881", "#A4A4A4",
                        "#96C4DD", "#1E65A7", "#262262"),
                      drug_annotation_class %>% pull())

row_ha <- ComplexHeatmap::HeatmapAnnotation(df = row_ha_df,
                                            col = list(Reversible = c("No" = "red", "Yes" = "black"),
                                                       `Class` = class_col),
                                            annotation_legend_param = list(
                                              Reversible = list(labels = c("Covalent", "Reversible"),
                                                                title_gp = gpar(fontsize = 6),
                                                                labels_gp = gpar(fontsize = 5),
                                                                grid_width = unit(0.3, "cm"),
                                                                grid_height = unit(0.3, "cm"),
                                                                ncol = 2),
                                              `Class` = list(title_gp = gpar(fontsize = 6),
                                                                    labels_gp = gpar(fontsize = 5),
                                                                    grid_width = unit(0.3, "cm"),
                                                                    grid_height = unit(0.3, "cm"),
                                                             ncol = 4)),
                                            annotation_name_gp = gpar(fontsize = 6),
                                            annotation_name_rot = 0,
                                            show_legend = c(TRUE, TRUE))

hmcols <- colorRampPalette(c("white","orange"))(2)
lgd <- Legend(title = "Membership", 
              title_gp = gpar(fontsize = 6),
              at = 0:1, 
              labels = c("Absent", "Present"),
              legend_gp = gpar(fill = c("white", "orange"),
                               lwd = 0.3),
              ncol = 1,
              border = "black",
              labels_gp = gpar(fontsize = 5),
              grid_width = unit(0.3, "cm"),
              grid_height = unit(0.3, "cm"))

# initial heatmap to define number of transcriptional drug class (TDCs) 
hm <- ComplexHeatmap::Heatmap(matrix = heatmap_mat,
                              cluster_columns = T, cluster_rows = T, 
                              col = hmcols,
                              # clustering_method_rows = "ward.D2",
                              # column_split = heatmap_mat_cluster_labels %>% column_to_rownames(var = "drug"),
                              column_split = 16,
                              row_split = 9,
                              row_dend_gp = gpar(lwd = 0.5),
                              column_dend_gp = gpar(lwd = 0.5),
                              column_gap = unit(0.2, "line"),
                              row_gap = unit(0.2, "line"),
                              row_title = "TDC",
                              # column_title_side = "bottom",
                              column_names_gp = gpar(fontsize = 5),
                              column_names_rot = 45,
                              row_title_gp = gpar(fontsize = 6),
                              row_names_gp = gpar(fontsize = 5),
                              column_title = NULL,
                              show_column_names = T,
                              top_annotation = row_ha,
                              rect_gp = gpar(col = "black", lwd = 0.3),
                              row_names_side = "right",
                              show_heatmap_legend = F)

# obtaining order of columns for larger TDCs
hm_draw <- ComplexHeatmap::draw(hm, annotation_legend_list = list(lgd), legend_grouping = "original", 
                                heatmap_legend_side = "bottom", annotation_legend_side = "bottom"); col_order_for_color <- column_order(hm_draw)

col_order_heatmap_mat <- data.frame("column" = colnames(heatmap_mat)) %>%
  rownames_to_column(var = "row_numb") 
col_order_heatmap_mat_final <- data.frame()
for(name in 1:length(col_order_for_color)) {
  col_order_heatmap_mat_temp <- col_order_heatmap_mat %>%
    filter(row_numb %in% col_order_for_color[[name]]) %>%
    mutate(`Transcriptional Drug Class` = name)
  col_order_heatmap_mat_final <- bind_rows(col_order_heatmap_mat_final, 
                                           col_order_heatmap_mat_temp)
}
col_order_heatmap_mat_final <- col_order_heatmap_mat_final %>%
  select(drug = column, `Transcriptional Drug Class`)

row_ha_df <- left_join(col_order_heatmap_mat_final, drug_annotations) %>%
  column_to_rownames(var = "drug")
row_ha_df <- row_ha_df[colnames(heatmap_mat),,drop = F]

color_scheme_rand <- readRDS("random_color_scheme_transcriptional_drug_class.RDS")

factor_col <- setNames(c(color_scheme_rand[1:length(col_order_for_color)]),
                       col_order_heatmap_mat_final %>%
                         distinct(`Transcriptional Drug Class`) %>%
                         arrange(`Transcriptional Drug Class`) %>%
                         pull())

row_ha <- ComplexHeatmap::HeatmapAnnotation(df = row_ha_df,
                                            na_col = "grey90",
                                            col = list(`Transcriptional Drug Class` = factor_col,
                                                       Reversible = c("No" = "red", "Yes" = "black", "NA" = "grey90"),
                                                       `Class` = class_col),
                                            annotation_legend_param = list(
                                              Reversible = list(labels = c("Covalent", "Reversible", "Unknown"),
                                                                at = c("No", "Yes", "NA"),
                                                                title_gp = gpar(fontsize = 6),
                                                                labels_gp = gpar(fontsize = 5),
                                                                grid_width = unit(0.3, "cm"),
                                                                grid_height = unit(0.3, "cm"),
                                                                ncol = 1),
                                              `Class` = list(title_gp = gpar(fontsize = 6),
                                                                    labels_gp = gpar(fontsize = 5),
                                                                    grid_width = unit(0.3, "cm"),
                                                                    grid_height = unit(0.3, "cm"),
                                                             nrow = 2)),
                                            annotation_name_gp = gpar(fontsize = 6),
                                            annotation_name_rot = 0,
                                            show_legend = c(FALSE, TRUE, TRUE))

col_split_temp <- col_order_heatmap_mat_final %>% column_to_rownames(var = "drug")
col_split_temp <- col_split_temp[colnames(heatmap_mat),,drop = F]

hm <- ComplexHeatmap::Heatmap(matrix = heatmap_mat,
                              cluster_columns = T, cluster_rows = T, 
                              col = hmcols,
                              # clustering_method_rows = "ward.D2",
                              column_split = col_split_temp,
                              row_split = 9,
                              row_dend_gp = gpar(lwd = 1),
                              column_dend_gp = gpar(lwd = 1),
                              column_gap = unit(0.2, "line"),
                              row_gap = unit(0.2, "line"),
                              row_title = "RM Cluster",
                              column_title_side = "bottom",
                              column_names_gp = gpar(fontsize = 5),
                              column_names_rot = 45,
                              row_title_gp = gpar(fontsize = 6),
                              row_names_gp = gpar(fontsize = 5),
                              column_title = "TDC%s",
                              column_title_gp = gpar(fontsize = 5),
                              column_title_rot = 45,
                              show_column_names = T,
                              top_annotation = row_ha,
                              rect_gp = gpar(col = "black", lwd = 1),
                              row_names_side = "right", 
                              row_dend_width = unit(0.3, "cm"),
                              show_heatmap_legend = F)

hm_draw <- ComplexHeatmap::draw(hm, annotation_legend_list = list(lgd), legend_grouping = "original", 
                                heatmap_legend_side = "bottom", annotation_legend_side = "bottom")

# write_csv(row_ha_df %>% rownames_to_column(var = "drug"),"SharedPDCL_MrVI_Drug_DRC_BinaryMembership_HierarchicalClusters_FDRfilt_TDCs.csv")

png("TDC_binary_membership_heatmap.png",width=12,height=5,units="in",res=1200)
draw(hm, annotation_legend_list = list(lgd),
     heatmap_legend_side = "right", annotation_legend_side = "top")
dev.off()



