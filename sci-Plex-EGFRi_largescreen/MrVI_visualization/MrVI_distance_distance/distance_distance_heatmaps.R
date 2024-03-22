library(monocle3)
library(tidyverse)
library(ComplexHeatmap)

this.dir <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(this.dir))
rm(this.dir)

# ================================================================================
# distance-distance heatmap remade in R
# ================================================================================

heatmap_mat <- read_csv("bt228_mrvi_normalized_dist_df_namedosesample_attention.csv") %>%
  dplyr::rename(X1 = 1) %>%
  column_to_rownames(var = "X1")
heatmap_mat[heatmap_mat > 2] <- 2

colnames_temp <- data.frame("X1" = colnames(heatmap_mat)) %>%
  separate(col = X1, into = c("drug", "dose"), sep = "_") %>%
  mutate(dose = case_when(
    !drug %in% c("DMSO", "Media", "PBS", "Panitumumab") ~ as.double(dose) - 1,
    TRUE ~ as.double(dose))) %>%
  unite(col = "X1", sep = "_", drug, dose) %>%
  pull()

colnames(heatmap_mat) <- colnames_temp
rownames(heatmap_mat) <- colnames_temp

heatmap_mat_cluster_labels <- read_csv("bt228_normalized_attention_drug_cluster_labels.csv") %>%
  dplyr::select(`Drug Cluster` = 2) %>%
  mutate(drug = colnames_temp) %>%
  separate(col = drug, into = c(NA, "Log10 Dose"), sep = "_", remove = F) %>%
  column_to_rownames(var = "drug")

color <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
color_scheme_rand <- sample(color, 20)

saveRDS(color_scheme_rand, "BT228_color_scheme_rand_20240124.RDS")
color_scheme_rand <- readRDS("BT228_color_scheme_rand_20240124.RDS")
factor_col <- setNames(c(color_scheme_rand),
                       heatmap_mat_cluster_labels %>%
                         distinct(`Drug Cluster`) %>%
                         arrange(`Drug Cluster`) %>%
                         pull())

row_ha <- ComplexHeatmap::rowAnnotation(df = heatmap_mat_cluster_labels,
                                        col = list(`Drug Cluster` = factor_col,
                                                   `Log10 Dose` = setNames(c("grey30", viridis::magma(n = 4)), 
                                                                           c("0", "1", "2", "3", "4"))),
                                        annotation_legend_param = list(
                                          `Log10 Dose` = list(title_gp = gpar(fontsize = 6),
                                                            labels_gp = gpar(fontsize = 5),
                                                            grid_width = unit(0.3, "cm"),
                                                            grid_height = unit(0.3, "cm"),
                                                            direction = "horizontal",
                                                            nrow = 1)),
                                        simple_anno_size = unit(0.3, "cm"),
                                        show_annotation_name = FALSE,
                                        show_legend = c(FALSE, TRUE))

col_ha <- ComplexHeatmap::HeatmapAnnotation(df = heatmap_mat_cluster_labels,
                                        col = list(`Drug Cluster` = factor_col,
                                                   `Log10 Dose` = setNames(c("grey30", viridis::magma(n = 4)), 
                                                                           c("0", "1", "2", "3", "4"))), 
                                        show_annotation_name = TRUE,
                                        simple_anno_size = unit(0.3, "cm"),
                                        annotation_name_gp = gpar(fontsize = 6),
                                        # annotation_name_rot = 45,
                                        show_legend = c(FALSE, FALSE))

# hmcols <- colorRampPalette(c("white","orange"))(2)
hmcols <- viridis::plasma(n = 50)

hm <- ComplexHeatmap::Heatmap(matrix = heatmap_mat %>% as.matrix(),
                              cluster_columns = T, cluster_rows = T, 
                              col = hmcols,
                              name = "Distance", 
                              heatmap_legend_param = list(title_gp = gpar(fontsize = 6),
                                                          legend_height = unit(0.1, "cm"),
                                                          legend_width = unit(1.5, "cm"),
                                                          labels_gp = gpar(fontsize = 5),
                                                          at = c(0,0.5,1,1.5,2),
                                                          labels = c("0", "0.5", "1", "1.5", ">2"), 
                                                          direction = "horizontal"),
                              top_annotation = col_ha,
                              column_split = heatmap_mat_cluster_labels %>% select(`Drug Cluster`),
                              column_gap = unit(0, "line"),
                              row_split = heatmap_mat_cluster_labels %>% select(`Drug Cluster`),
                              row_gap = unit(0, "line"),
                              column_title = NULL,
                              column_title_side = "bottom",
                              column_title_gp = gpar(fontsize = 6),
                              column_names_gp = gpar(fontsize = 5),
                              row_title = NULL,
                              show_column_names = F,
                              show_row_names = F,
                              left_annotation = row_ha,
                              # rect_gp = gpar(col = "black", lwd = 0.3),
                              row_names_gp = gpar(fontsize = 5), 
                              row_names_side = "right",
                              row_dend_gp = gpar(lwd = 0.5),
                              column_dend_gp = gpar(lwd = 0.5),
                              show_heatmap_legend = T)

ComplexHeatmap::draw(hm, heatmap_legend_side = "bottom", 
                     annotation_legend_side = "bottom")

png("BT228_distance_distance_plasma_legend_bottom.png",width=4,height=4,units="in",res=1200)
draw(hm, heatmap_legend_side = "bottom", 
     annotation_legend_side = "bottom")
dev.off()

# ================================================================================
# distance-distance heatmap remade in R (sampled)
# ================================================================================

heatmap_mat <- read_csv("bt228_mrvi_normalized_dist_df_namedosesample_attention.csv") %>%
  dplyr::rename(X1 = 1) %>%
  column_to_rownames(var = "X1")
heatmap_mat[heatmap_mat > 2] <- 2

colnames_temp <- data.frame("X1" = colnames(heatmap_mat)) %>%
  separate(col = X1, into = c("drug", "dose"), sep = "_") %>%
  mutate(dose = case_when(
    !drug %in% c("DMSO", "Media", "PBS", "Panitumumab") ~ as.double(dose) - 1,
    TRUE ~ as.double(dose))) %>%
  unite(col = "X1", sep = "_", drug, dose) %>%
  pull()

colnames(heatmap_mat) <- colnames_temp
rownames(heatmap_mat) <- colnames_temp

heatmap_mat_cluster_labels <- read_csv("bt228_normalized_attention_drug_cluster_labels.csv") %>%
  dplyr::select(`Drug Cluster` = 2) %>%
  mutate(drug = colnames_temp) %>%
  separate(col = drug, into = c(NA, "Log10 Dose"), sep = "_", remove = F) 

heatmap_mat_cluster_labels_sampled <- heatmap_mat_cluster_labels %>%
  group_by(`Drug Cluster`) %>%
  slice_sample(n = 10) %>%
  column_to_rownames(var = "drug")

# color <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
# color_scheme_rand <- sample(color, 20)

# saveRDS(color_scheme_rand, "BT228_color_scheme_rand_20240124.RDS")
color_scheme_rand <- readRDS("BT228_color_scheme_rand_20240124.RDS")
factor_col <- setNames(c(color_scheme_rand),
                       heatmap_mat_cluster_labels_sampled %>%
                         distinct(`Drug Cluster`) %>%
                         arrange(`Drug Cluster`) %>%
                         pull())

row_ha <- ComplexHeatmap::rowAnnotation(df = heatmap_mat_cluster_labels_sampled,
                                        col = list(`Drug Cluster` = factor_col,
                                                   `Log10 Dose` = setNames(c("grey30", viridis::magma(n = 4)), 
                                                                           c("0", "1", "2", "3", "4"))),
                                        annotation_legend_param = list(
                                          `Log10 Dose` = list(title_gp = gpar(fontsize = 6),
                                                              labels_gp = gpar(fontsize = 5),
                                                              grid_width = unit(0.3, "cm"),
                                                              grid_height = unit(0.3, "cm"),
                                                              direction = "horizontal",
                                                              nrow = 1)),
                                        simple_anno_size = unit(0.3, "cm"),
                                        show_annotation_name = FALSE,
                                        show_legend = c(FALSE, TRUE))

col_ha <- ComplexHeatmap::HeatmapAnnotation(df = heatmap_mat_cluster_labels_sampled,
                                            col = list(`Drug Cluster` = factor_col,
                                                       `Log10 Dose` = setNames(c("grey30", viridis::magma(n = 4)), 
                                                                               c("0", "1", "2", "3", "4"))), 
                                            show_annotation_name = TRUE,
                                            simple_anno_size = unit(0.3, "cm"),
                                            annotation_name_gp = gpar(fontsize = 6),
                                            # annotation_name_rot = 45,
                                            show_legend = c(FALSE, FALSE))

# hmcols <- colorRampPalette(c("white","orange"))(2)
hmcols <- viridis::plasma(n = 50)

heatmap_mat_sampled <- heatmap_mat[rownames(heatmap_mat_cluster_labels_sampled),
                           rownames(heatmap_mat_cluster_labels_sampled)]

hm <- ComplexHeatmap::Heatmap(matrix = heatmap_mat_sampled %>% as.matrix(),
                              cluster_columns = T, cluster_rows = T, 
                              col = hmcols,
                              name = "Distance", 
                              heatmap_legend_param = list(title_gp = gpar(fontsize = 6),
                                                          legend_height = unit(0.1, "cm"),
                                                          legend_width = unit(1.5, "cm"),
                                                          labels_gp = gpar(fontsize = 5),
                                                          at = c(0,0.5,1,1.5,2),
                                                          labels = c("0", "0.5", "1", "1.5", ">2"), 
                                                          direction = "horizontal"),
                              top_annotation = col_ha,
                              column_split = heatmap_mat_cluster_labels_sampled %>% select(`Drug Cluster`),
                              column_gap = unit(0, "line"),
                              row_split = heatmap_mat_cluster_labels_sampled %>% select(`Drug Cluster`),
                              row_gap = unit(0, "line"),
                              column_title = NULL,
                              column_title_side = "bottom",
                              column_title_gp = gpar(fontsize = 6),
                              column_names_gp = gpar(fontsize = 5),
                              row_title = NULL,
                              show_column_names = F,
                              show_row_names = F,
                              left_annotation = row_ha,
                              # rect_gp = gpar(col = "black", lwd = 0.3),
                              row_names_gp = gpar(fontsize = 5), 
                              row_names_side = "right",
                              row_dend_gp = gpar(lwd = 0.5),
                              column_dend_gp = gpar(lwd = 0.5),
                              show_heatmap_legend = T)

ComplexHeatmap::draw(hm, heatmap_legend_side = "bottom", 
                     annotation_legend_side = "bottom")

png("BT228_distance_distance_plasma_legend_bottom_sampled_10.png",width=4,height=4,units="in",res=1200)
draw(hm, heatmap_legend_side = "bottom", 
     annotation_legend_side = "bottom")
dev.off()

# ================================================================================
# distance-distance heatmap remade in R (sampled) - BT112
# ================================================================================

heatmap_mat <- read_csv("bt112_mrvi_normalized_dist_df_namedosesample_attention.csv") %>%
  dplyr::rename(X1 = 1) %>%
  column_to_rownames(var = "X1")
heatmap_mat[heatmap_mat > 2] <- 2

colnames_temp <- data.frame("X1" = colnames(heatmap_mat)) %>%
  separate(col = X1, into = c("drug", "dose"), sep = "_") %>%
  mutate(dose = case_when(
    !drug %in% c("DMSO", "Media", "PBS", "Panitumumab") ~ as.double(dose) - 1,
    TRUE ~ as.double(dose))) %>%
  unite(col = "X1", sep = "_", drug, dose) %>%
  pull()

colnames(heatmap_mat) <- colnames_temp
rownames(heatmap_mat) <- colnames_temp

heatmap_mat_cluster_labels <- read_csv("bt112_normalized_attention_drug_cluster_labels.csv") %>%
  dplyr::select(`Drug Cluster` = 2) %>%
  mutate(drug = colnames_temp) %>%
  separate(col = drug, into = c(NA, "Log10 Dose"), sep = "_", remove = F) 

heatmap_mat_cluster_labels_sampled <- heatmap_mat_cluster_labels %>%
  group_by(`Drug Cluster`) %>%
  slice_sample(n = 10) %>%
  column_to_rownames(var = "drug")

color <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
color_scheme_rand <- sample(color, 23)

saveRDS(color_scheme_rand, "BT112_color_scheme_rand_20240124.RDS")
# color_scheme_rand <- readRDS("BT228_color_scheme_rand_20240124.RDS")
factor_col <- setNames(c(color_scheme_rand),
                       heatmap_mat_cluster_labels_sampled %>%
                         distinct(`Drug Cluster`) %>%
                         arrange(`Drug Cluster`) %>%
                         pull())

row_ha <- ComplexHeatmap::rowAnnotation(df = heatmap_mat_cluster_labels_sampled,
                                        col = list(`Drug Cluster` = factor_col,
                                                   `Log10 Dose` = setNames(c("grey30", viridis::magma(n = 4)), 
                                                                           c("0", "1", "2", "3", "4"))),
                                        annotation_legend_param = list(
                                          `Log10 Dose` = list(title_gp = gpar(fontsize = 6),
                                                              labels_gp = gpar(fontsize = 5),
                                                              grid_width = unit(0.3, "cm"),
                                                              grid_height = unit(0.3, "cm"),
                                                              direction = "horizontal",
                                                              nrow = 1)),
                                        simple_anno_size = unit(0.3, "cm"),
                                        show_annotation_name = FALSE,
                                        show_legend = c(FALSE, TRUE))

col_ha <- ComplexHeatmap::HeatmapAnnotation(df = heatmap_mat_cluster_labels_sampled,
                                            col = list(`Drug Cluster` = factor_col,
                                                       `Log10 Dose` = setNames(c("grey30", viridis::magma(n = 4)), 
                                                                               c("0", "1", "2", "3", "4"))), 
                                            show_annotation_name = TRUE,
                                            simple_anno_size = unit(0.3, "cm"),
                                            annotation_name_gp = gpar(fontsize = 6),
                                            # annotation_name_rot = 45,
                                            show_legend = c(FALSE, FALSE))

# hmcols <- colorRampPalette(c("white","orange"))(2)
hmcols <- viridis::plasma(n = 50)

heatmap_mat_sampled <- heatmap_mat[rownames(heatmap_mat_cluster_labels_sampled),
                                   rownames(heatmap_mat_cluster_labels_sampled)]

hm <- ComplexHeatmap::Heatmap(matrix = heatmap_mat_sampled %>% as.matrix(),
                              cluster_columns = T, cluster_rows = T, 
                              col = hmcols,
                              name = "Distance", 
                              heatmap_legend_param = list(title_gp = gpar(fontsize = 6),
                                                          legend_height = unit(0.1, "cm"),
                                                          legend_width = unit(1.5, "cm"),
                                                          labels_gp = gpar(fontsize = 5),
                                                          at = c(0,0.5,1,1.5,2),
                                                          labels = c("0", "0.5", "1", "1.5", ">2"),
                                                          # at = c(0,1,2,3),
                                                          # labels = c("0", "1", "2", ">3"), 
                                                          direction = "horizontal"),
                              top_annotation = col_ha,
                              column_split = heatmap_mat_cluster_labels_sampled %>% select(`Drug Cluster`),
                              column_gap = unit(0, "line"),
                              row_split = heatmap_mat_cluster_labels_sampled %>% select(`Drug Cluster`),
                              row_gap = unit(0, "line"),
                              column_title = NULL,
                              column_title_side = "bottom",
                              column_title_gp = gpar(fontsize = 6),
                              column_names_gp = gpar(fontsize = 5),
                              row_title = NULL,
                              show_column_names = F,
                              show_row_names = F,
                              left_annotation = row_ha,
                              # rect_gp = gpar(col = "black", lwd = 0.3),
                              row_names_gp = gpar(fontsize = 5), 
                              row_names_side = "right",
                              row_dend_gp = gpar(lwd = 0.5),
                              column_dend_gp = gpar(lwd = 0.5),
                              show_heatmap_legend = T)

ComplexHeatmap::draw(hm, heatmap_legend_side = "bottom", 
                     annotation_legend_side = "bottom")

png("BT112_distance_distance_plasma_legend_bottom_sampled_10.png",
    width=4,height=4,units="in",res=1200)
draw(hm, heatmap_legend_side = "bottom", 
     annotation_legend_side = "bottom")
dev.off()

# ================================================================================
# distance-distance heatmap remade in R (sampled) - BT333
# ================================================================================

heatmap_mat <- read_csv("bt333_mrvi_normalized_dist_df_namedosesample_attention.csv") %>%
  dplyr::rename(X1 = 1) %>%
  column_to_rownames(var = "X1")
heatmap_mat[heatmap_mat > 3] <- 3

colnames_temp <- data.frame("X1" = colnames(heatmap_mat)) %>%
  separate(col = X1, into = c("drug", "dose"), sep = "_") %>%
  mutate(dose = case_when(
    !drug %in% c("DMSO", "Media", "PBS", "Panitumumab") ~ as.double(dose) - 1,
    TRUE ~ as.double(dose))) %>%
  unite(col = "X1", sep = "_", drug, dose) %>%
  pull()

colnames(heatmap_mat) <- colnames_temp
rownames(heatmap_mat) <- colnames_temp

heatmap_mat_cluster_labels <- read_csv("bt333_normalized_attention_drug_cluster_labels.csv") %>%
  dplyr::select(`Drug Cluster` = 2) %>%
  mutate(drug = colnames_temp) %>%
  separate(col = drug, into = c(NA, "Log10 Dose"), sep = "_", remove = F) 

heatmap_mat_cluster_labels_sampled <- heatmap_mat_cluster_labels %>%
  group_by(`Drug Cluster`) %>%
  slice_sample(n = 10) %>%
  column_to_rownames(var = "drug")

color <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
color_scheme_rand <- sample(color, 25)

saveRDS(color_scheme_rand, "BT333_color_scheme_rand_20240124.RDS")
# color_scheme_rand <- readRDS("BT228_color_scheme_rand_20240124.RDS")
factor_col <- setNames(c(color_scheme_rand),
                       heatmap_mat_cluster_labels_sampled %>%
                         distinct(`Drug Cluster`) %>%
                         arrange(`Drug Cluster`) %>%
                         pull())

row_ha <- ComplexHeatmap::rowAnnotation(df = heatmap_mat_cluster_labels_sampled,
                                        col = list(`Drug Cluster` = factor_col,
                                                   `Log10 Dose` = setNames(c("grey30", viridis::magma(n = 4)), 
                                                                           c("0", "1", "2", "3", "4"))),
                                        annotation_legend_param = list(
                                          `Log10 Dose` = list(title_gp = gpar(fontsize = 6),
                                                              labels_gp = gpar(fontsize = 5),
                                                              grid_width = unit(0.3, "cm"),
                                                              grid_height = unit(0.3, "cm"),
                                                              direction = "horizontal",
                                                              nrow = 1)),
                                        simple_anno_size = unit(0.3, "cm"),
                                        show_annotation_name = FALSE,
                                        show_legend = c(FALSE, TRUE))

col_ha <- ComplexHeatmap::HeatmapAnnotation(df = heatmap_mat_cluster_labels_sampled,
                                            col = list(`Drug Cluster` = factor_col,
                                                       `Log10 Dose` = setNames(c("grey30", viridis::magma(n = 4)), 
                                                                               c("0", "1", "2", "3", "4"))), 
                                            show_annotation_name = TRUE,
                                            simple_anno_size = unit(0.3, "cm"),
                                            annotation_name_gp = gpar(fontsize = 6),
                                            # annotation_name_rot = 45,
                                            show_legend = c(FALSE, FALSE))

# hmcols <- colorRampPalette(c("white","orange"))(2)
hmcols <- viridis::plasma(n = 50)

heatmap_mat_sampled <- heatmap_mat[rownames(heatmap_mat_cluster_labels_sampled),
                                   rownames(heatmap_mat_cluster_labels_sampled)]

hm <- ComplexHeatmap::Heatmap(matrix = heatmap_mat_sampled %>% as.matrix(),
                              cluster_columns = T, cluster_rows = T, 
                              col = hmcols,
                              name = "Distance", 
                              heatmap_legend_param = list(title_gp = gpar(fontsize = 6),
                                                          legend_height = unit(0.1, "cm"),
                                                          legend_width = unit(1.5, "cm"),
                                                          labels_gp = gpar(fontsize = 5),
                                                          # at = c(0,0.5,1,1.5,2),
                                                          # labels = c("0", "0.5", "1", "1.5", ">2"), 
                                                          at = c(0,1,2,3),
                                                          labels = c("0", "1", "2", ">3"),
                                                          direction = "horizontal"),
                              top_annotation = col_ha,
                              column_split = heatmap_mat_cluster_labels_sampled %>% select(`Drug Cluster`),
                              column_gap = unit(0, "line"),
                              row_split = heatmap_mat_cluster_labels_sampled %>% select(`Drug Cluster`),
                              row_gap = unit(0, "line"),
                              column_title = NULL,
                              column_title_side = "bottom",
                              column_title_gp = gpar(fontsize = 6),
                              column_names_gp = gpar(fontsize = 5),
                              row_title = NULL,
                              show_column_names = F,
                              show_row_names = F,
                              left_annotation = row_ha,
                              # rect_gp = gpar(col = "black", lwd = 0.3),
                              row_names_gp = gpar(fontsize = 5), 
                              row_names_side = "right",
                              row_dend_gp = gpar(lwd = 0.5),
                              column_dend_gp = gpar(lwd = 0.5),
                              show_heatmap_legend = T)

ComplexHeatmap::draw(hm, heatmap_legend_side = "bottom", 
                     annotation_legend_side = "bottom")

png("BT333_distance_distance_plasma_legend_bottom_sampled_10.png",
    width=4,height=4,units="in",res=1200)
draw(hm, heatmap_legend_side = "bottom", 
     annotation_legend_side = "bottom")
dev.off()
