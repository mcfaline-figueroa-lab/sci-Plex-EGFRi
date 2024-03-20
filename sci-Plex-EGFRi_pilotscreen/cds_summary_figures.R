library(tidyverse)
library(monocle3)

cds <- readRDS("cds_pilot_screen.RDS")
correct_cells <- readRDS("GBM_cells_filtered.rds")
cds <- cds[,correct_cells]

saveRDS(cds, "cds_pilot_screen.RDS")

# ==============================================================================
# Cells per cell line
# ==============================================================================
cells_total <- colData(cds) %>% as.data.frame() %>% 
  summarise(`Cell Count` = n(),
            `Median UMI` = median(n.umi),
            `Median Hash UMI` = median(hash_umis), .groups = "drop") %>%
  mutate(`Cell Line` = "Total") %>%
  select(`Cell Line`, everything())

cells_per_conditions <- colData(cds) %>% as.data.frame() %>%
  group_by(Cell.Line) %>%
  summarise(`Cell Count` = n(),
            `Median UMI` = median(n.umi),
            `Median Hash UMI` = median(hash_umis), .groups = "drop") %>%
  dplyr::rename(`Cell Line` = Cell.Line)

cells_per_conditions_merge <- bind_rows(cells_per_conditions,
                                        cells_total)

# ==============================================================================
# Cells/DEGs per drug
# ==============================================================================

cells_per_conditions <- colData(cds) %>% as.data.frame() %>%
  group_by(treatment) %>%
  summarise(`Cell Count` = n()) %>%
  dplyr::rename(Drug = treatment)

diff_test_result <- read_csv("supplementary_table_S1_pilot_screen_dose_DEG_fit_models_result.csv")

DEG_filtered <- diff_test_result %>%
  filter(grepl("log", term) & q_value < 0.001 & abs(normalized_effect)>0.05) %>%
  pull(id) %>%
  unique()

DEG_filtered_per_drug_up <- diff_test_result %>%
  filter(term == "log(dose + 0.1)", 
         id %in% DEG_filtered, 
         normalized_effect > 0, 
         q_value < 0.001, 
         !is.na(normalized_effect)) %>%
  group_by(treatment) %>%
  summarise(count = n(), .groups = "drop")

DEG_filtered_per_drug_dn <- diff_test_result %>%
  filter(term == "log(dose + 0.1)", 
         id %in% DEG_filtered, 
         normalized_effect < 0, 
         q_value < 0.001, 
         !is.na(normalized_effect)) %>%
  group_by(treatment) %>%
  summarise(count = n(), .groups = "drop")

DEG_filtered_per_drug_tot <- diff_test_result %>%
  filter(term == "log(dose + 0.1)", 
         id %in% DEG_filtered, 
         q_value < 0.001, 
         !is.na(normalized_effect)) %>%
  group_by(treatment) %>%
  summarise(count = n(), .groups = "drop")

# ==============================================================================
# Volcano plots stacked
# ==============================================================================

ggplot(diff_test_result %>% filter(term == "log(dose + 0.1)"),
       aes(x = normalized_effect, y = -log10(q_value), 
           color = q_value < 0.001 & abs(normalized_effect) > 0.05)) +
  geom_point(size = 0.2) +
  geom_vline(xintercept = 0.05, linetype = "dashed", size = 0.1) +
  geom_vline(xintercept = -0.05, linetype = "dashed", size = 0.1) +
  monocle3:::monocle_theme_opts() +
  facet_wrap(~factor(Cell.Line, levels = c("A172", "BT333", "T98G", "U87MG")), 
             ncol = 2,
             scales = "free_y") +
  xlab("Normalized Effect") +
  ylab("-Log10 FDR") +
  scale_color_manual("FDR < 0.1%", values = c("TRUE" = "firebrick3", "FALSE" = "black")) +
  theme(text = element_text(size = 6),
        legend.position = "bottom",
        legend.margin = margin(t = -3, r = 0, b = 0, l = 0, unit = "mm")) +
        # legend.spacing.y = unit(0.1, "cm")) +
  guides(color = guide_legend(title = "FDR < 0.1%", byrow = TRUE, override.aes = list(size=3)))
ggsave("/Users/rossgiglio/Documents/McFaline-Figueroa_Lab/Experiment/4lines_combined/sup_fig_1/DEG_volcano_facet_line.png", 
       width = 3, height = 2, dpi = 600)

ggplot(diff_test_result %>% filter(term == "log(dose + 0.1)"),
       aes(x = normalized_effect, y = -log10(q_value), 
           color = q_value < 0.001 & abs(normalized_effect) > 0.05)) +
  geom_point(size = 0.2) +
  geom_vline(xintercept = 0.05, linetype = "dashed", size = 0.1) +
  geom_vline(xintercept = -0.05, linetype = "dashed", size = 0.1) +
  monocle3:::monocle_theme_opts() +
  facet_wrap(~factor(Cell.Line, levels = c("A172", "BT333", "T98G", "U87MG")), 
             nrow = 1,
             scales = "free_y") +
  xlab("Normalized Effect") +
  ylab("-Log10 FDR") +
  scale_color_manual("FDR < 0.1%", values = c("TRUE" = "firebrick3", "FALSE" = "black")) +
  theme(text = element_text(size = 10),
        strip.text = element_text(size = 9),
        legend.position = "bottom",
        legend.margin = margin(t = -3, r = 0, b = 0, l = 0, unit = "mm")) +
  # legend.spacing.y = unit(0.1, "cm")) +
  guides(color = guide_legend(title = "FDR < 0.1%", byrow = TRUE, override.aes = list(size=3)))
ggsave("DEG_volcano_facet_line.png", 
       width = 6, height = 2, dpi = 600)

# ==============================================================================
# Covariance matrices
# ==============================================================================
deg_list <- diff_test_result %>% 
  filter(grepl(pattern = "log", term), 
         q_value < 0.001) %>%
  pull(id) %>%
  unique()

deg_precorr <- list()
deg_corr <- list()
for (cell_type in c("A172", "T98G", "U87MG", "BT333")){

  deg_precorr[[cell_type]] <- diff_test_result %>%
    filter(Cell.Line == cell_type) %>%
    filter(grepl("log", term),id %in% deg_list) %>%
    # mutate(condition = paste0(Cell.Line,"_",treatment)) %>%
    dplyr::select(condition = treatment, id, normalized_effect) %>%
    tidyr::spread(key = condition, value = normalized_effect) %>%
    as.data.frame()
  
  row.names(deg_precorr[[cell_type]]) <- deg_precorr[[cell_type]]$id
  deg_precorr[[cell_type]]$id <- NULL
  
  deg_corr[[cell_type]] <- cor(deg_precorr[[cell_type]], method = "pearson")
  
  hmcols = colorRampPalette(c("white","red"))(50)
  hm <- ComplexHeatmap::Heatmap(deg_corr[[cell_type]],
                          col = hmcols,
                          row_title = NULL, 
                          column_title = paste0(cell_type), column_title_gp = gpar(fontsize = 10),
                          name = "Pearson\nCorrelation (r)",
                          rect_gp = gpar(col = "black", lwd = 1),
                          row_dend_gp = gpar(lwd = 1),
                          column_dend_gp = gpar(lwd = 1),
                          column_split = 4,
                          row_split = 4,
                          row_dend_width = unit(0.5, "cm"),
                          column_dend_height = unit(0.5, "cm"),
                          heatmap_legend_param = list(border = "black",
                                                      title_gp = gpar(fontsize = 6),
                                                      title_position = "leftcenter",
                                                      at = seq(0,1,0.25), 
                                                      labels = seq(0,1,0.25),
                                                      labels_gp = gpar(fontsize = 5),
                                                      direction = "horizontal",
                                                      legend_height = unit(0.1, "cm"),
                                                      legend_width = unit(2.5, "cm")
                                                      ),
                          row_names_gp = gpar(fontsize = 7),
                          column_names_gp = gpar(fontsize = 7)
                          )
  
  png(paste0(cell_type,
             "_DEG_pearson_correlation_legend.png"),
      width=2.2,height=2,units="in",res=900)
  ComplexHeatmap::draw(hm, heatmap_legend_side = "bottom")
  dev.off()
  
}

# ==============================================================================
# Intersection Upset plots and z-scored aggregate signature score heatmap 
# ==============================================================================
library(ComplexHeatmap)

upset.list <- list()
for (drug in c("Afatinib", "Brigatinib", "CUDC-101", "Neratinib", "Osimertinib")) {
  temp <- diff_test_result %>%
    filter(term == "log(dose + 0.1)", treatment == drug, id %in% DEG_filtered, normalized_effect > 0, q_value < 0.001, !is.na(normalized_effect)) %>%
    # filter(term == "log10dose_pseudo", treatment == drug, id %in% DEG_filtered, normalized_effect > 0, q_value < 0.001, !is.na(normalized_effect)) %>%
    select(gene_short_name) %>%
    distinct() %>%
    pull()
  upset.list[[drug]] <- temp
}

intersection <- ComplexHeatmap::make_comb_mat(upset.list, mode = "distinct")

intersect_plot <- ComplexHeatmap::UpSet(intersection[comb_size(intersection) >= 30], 
                                        comb_order = order(comb_size(intersection[comb_size(intersection) >= 30])),
                                        pt_size = unit(0.1, "cm"), lwd = unit(0.5, "mm"), 
                                        row_names_gp = gpar(fontsize = 6),
                                        # set_order = order(set_name(A172)),
                                        top_annotation = upset_top_annotation(intersection[comb_size(intersection) >= 30],  
                                                                              # ylim = c(0, 500), 
                                                                              add_numbers = TRUE, 
                                                                              numbers_gp = gpar(fontsize = 5), 
                                                                              height = unit(0.75, "cm"), 
                                                                              axis_param = list(gp = gpar(fontsize = 6)),
                                                                              annotation_name_gp = gpar(fontsize = 6)),
                                        right_annotation = upset_right_annotation(intersection[comb_size(intersection) >= 30], 
                                                                                  add_numbers = FALSE, 
                                                                                  width = unit(0.75, "cm"),
                                                                                  axis_param = list(gp = gpar(fontsize = 6)),
                                                                                  annotation_name_gp = gpar(fontsize = 6), 
                                                                                  show_annotation_name = T
                                        ),
                                        #column_title = " DEG Intersection")
)
intersect_plot

png(file="DEG_up_intersection_plot.png", 
    width=3,height=1.5,units="in",
    res = 900)
ComplexHeatmap::draw(intersect_plot)
dev.off()

# transposed
intersect_plot <- ComplexHeatmap::UpSet(t(intersection[comb_size(intersection) >= 30]), 
                                        comb_order = order(comb_size(t(intersection[comb_size(intersection) >= 30]))),
                                        pt_size = unit(0.1, "cm"), lwd = unit(0.5, "mm"), 
                                        column_names_gp = gpar(fontsize = 6),
                                        column_names_rot = 45,
                                        # set_order = order(set_name(A172)),
                                        top_annotation = upset_top_annotation(t(intersection[comb_size(intersection) >= 30]),  
                                                                              # ylim = c(0, 500), 
                                                                              add_numbers = F, 
                                                                              numbers_gp = gpar(fontsize = 5), 
                                                                              height = unit(0.75, "cm"), 
                                                                              axis_param = list(gp = gpar(fontsize = 6)),
                                                                              annotation_name_gp = gpar(fontsize = 6)),
                                        right_annotation = upset_right_annotation(t(intersection[comb_size(intersection) >= 30]),
                                                                                  title = "# DEGs",
                                                                                  add_numbers = TRUE, 
                                                                                  numbers_gp = gpar(fontsize = 5), 
                                                                                  width = unit(0.75, "cm"),
                                                                                  axis_param = list(gp = gpar(fontsize = 6)),
                                                                                  annotation_name_gp = gpar(fontsize = 6), 
                                                                                  show_annotation_name = T
                                        ),
                                        #column_title = " DEG Intersection")
)
intersect_plot

png(file="DEG_up_intersection_plot_transposed.png", 
    width=1.8,height=3,units="in",
    res = 900)
ComplexHeatmap::draw(intersect_plot)
dev.off()

# signature heatmap
source("calculate_aggreg_expression.R")
gene_sets <- readRDS("pilot_screen_EGFRi_intersection_sigantures.RDS")

cds.list <- list()
for (cell_type in c("A172", "T98G", "U87MG", "BT333")) {
  cds.list[[cell_type]] <- cds[,colData(cds)$Cell.Line == cell_type & 
                                 !colData(cds)$treatment %in% c("Water", "Puromycin", "EAI045")]
}

col_info_merge <- colData(cds) %>% 
  as.data.frame() %>% 
  select(cell, Cell.Line, treatment, dose) %>%
  unite(col = "id", sep = "_", Cell.Line, treatment, dose)

sig_mat <- data.frame()
# sig_mat <- list()
for (cell_type in c("A172", "T98G", "U87MG", "BT333")) {
  sig_mat_cell <- data.frame()
  for (name in names(gene_sets)) {
    temp <- calculate_aggregate_expression_score(cds = cds.list[[cell_type]],
                                                 signature_genes = gene_sets[[name]],
                                                 from_id = FALSE) %>% as.data.frame() %>% dplyr::rename(V1 = 1)
    temp$cell <- row.names(temp)
    temp_median <- left_join(temp, col_info_merge, by = c("cell" = "cell"))
    temp_median <- temp_median %>% group_by(id) %>% 
      summarise(med_score = median(V1)) %>% 
      mutate(signature = name) 
    sig_mat_cell <- bind_rows(sig_mat_cell, temp_median)
    
  }
  sig_mat_cell <- pivot_wider(data = sig_mat_cell, 
                              names_from = signature, 
                              values_from = med_score,
                              id_cols = id) %>%
    column_to_rownames("id")
  agg_scaled <- t(scale(t(scale(sig_mat_cell))))
  agg_scaled[agg_scaled > 2] <- 2
  agg_scaled[agg_scaled < -2] <- -2
  agg_scaled[is.na(agg_scaled)] <- 0
  agg_scaled <- as.data.frame(agg_scaled)
  sig_mat <- bind_rows(sig_mat, agg_scaled)
  # sig_mat[[cell_type]] <- agg_scaled
  
}

rm(temp, agg_scaled, name, cell_type, cds.list, cds)


sig_mat_wide <- sig_mat %>% t() 

col_anno <- data.frame(x = colnames(sig_mat_wide)) %>%
  mutate(row = x) %>%
  separate(x, into = c("Cell Line", "Treatment", "Dose"), sep = "_") %>%
  mutate(Dose = as.numeric(Dose)) %>%
  mutate(`Log10(Dose)` = log10(Dose + 0.1))
row.names(col_anno) <- col_anno$row
col_anno <- select(col_anno, `Cell Line`, Drug = Treatment, `Log10 Dose` = `Log10(Dose)`)

treat_var <- setNames(RColorBrewer::brewer.pal(n = 6, "Accent"), 
                      c("Afatinib", "Brigatinib", "CUDC-101", "Neratinib", "Osimertinib", "DMSO"))

dose_var <- circlize::colorRamp2(c(-1,0,1,2,3,4), 
                       viridis::magma(n = 6))
cell_var <- setNames(RColorBrewer::brewer.pal(n = 4, "Dark2"),
                     c("A172", "T98G", "U87MG", "BT333"))

anno_colors = list(Drug = treat_var, `Log10 Dose` = dose_var, `Cell Line` = cell_var)

anno_complex <- HeatmapAnnotation(df = col_anno,
                                  col = anno_colors, 
                                  annotation_name_side = "right", 
                                  annotation_name_gp = gpar(fontsize = 6),
                                  annotation_legend_param = list(`Log10 Dose` = list(at = c(-1,0,1,2,3,4),
                                                                                     labels = c("-1","0","1","2","3","4"),
                                                                                     title_gp = gpar(fontsize = 6),
                                                                                     labels_gp = gpar(fontsize = 5),
                                                                                     border = "black"),
                                                                 `Cell Line` = list(title_gp = gpar(fontsize = 6),
                                                                                    labels_gp = gpar(fontsize = 5),
                                                                                    border = "black"),
                                                                 Drug = list(title_gp = gpar(fontsize = 6),
                                                                             labels_gp = gpar(fontsize = 5),
                                                                             border = "black"))
                                  )

hm <- ComplexHeatmap::Heatmap(matrix = sig_mat_wide, name = "Median Signature\nExpression (z-scored)",
                        clustering_method_columns = "ward.D2", clustering_method_rows = "ward.D2",
                        show_column_names = F,
                        row_names_gp = gpar(fontsize = 6),
                        column_title = "Cell Line, Drug, Dose", column_title_side = "bottom", column_title_gp = gpar(fontsize = 6),
                        row_title = "Signature (DEG Intersection)", row_title_gp = gpar(fontsize = 6),
                        col = colorRampPalette(c("navy","white","red"))(35),
                        heatmap_legend_param = list(title_gp = gpar(fontsize = 6),
                                                    title_position = "leftcenter",
                                                    labels_gp = gpar(fontsize = 5),
                                                    direction = "horizontal",
                                                    legend_width = unit(2.5, "cm"),
                                                    legend_height = unit(0.1, "cm"),
                                                    border = "black"
                                                    ),
                        top_annotation = anno_complex,
                        column_dend_height = unit(0.4, "cm"),
                        row_dend_width = unit(0.4, "cm"),
                        row_split = 6,
                        column_split = 6)

png(file="aggregate_signature_score_z-scored_heatmap.png", 
    width=4,height=4,units="in",
    res = 900)
draw(hm, 
     heatmap_legend_side = "top", 
     annotation_legend_side = "right",
     legend_grouping = "adjusted")
dev.off()

# ==============================================================================
# Overlap with MEK adaptive resistance (sci-Plex-GxE)
# ==============================================================================

adaptive_program_up <- read_csv("sci-Plex-GxE_adaptive_program.csv",
                                col_names = TRUE) %>%
  filter(grepl("pregulated", gene_module))

gene_sets <- readRDS("pilot_screen_EGFRi_intersection_sigantures.RDS")

overlap_func <- function(a, b) {
  intersection = length(intersect(a, b))
  minimum = min(length(a), length(b))
  union = length(a) + length(b) - intersection
  # return (intersection/union) # jaccard coefficient
  return (intersection/minimum)
}

overlap_res <- data.frame()
gene_set_size <- data.frame()
for (set in names(gene_sets)) {
  overlap_temp <- data.frame("set" = set, 
                             "overlap_up" = overlap_func(adaptive_program_up$gene_short_name, 
                                                         gene_sets[[set]])
                            )
  overlap_res <- bind_rows(overlap_res, overlap_temp)
  gene_set_size_temp <- data.frame("set" = set,
                                   "genes_in_sig" = length(gene_sets[[set]]))
  gene_set_size <- bind_rows(gene_set_size, gene_set_size_temp)
}

overlap_res_heatmap <- overlap_res %>%
  # mutate(difference = jaccard_up - jaccard_down) %>%
  select(set, overlap_up) %>%
  arrange(desc(overlap_up)) %>%
  column_to_rownames(var = "set") %>%
  as.matrix()
gene_set_size_ordered <- gene_set_size %>%
  mutate(genes_in_sig = case_when(
    genes_in_sig > 300 ~ 300,
    TRUE ~ genes_in_sig
  )) %>%
  column_to_rownames(var = "set")
gene_set_size_ordered <- gene_set_size_ordered[row.names(overlap_res_heatmap),]
row_bar <- rowAnnotation(`Set Size` = anno_barplot(gene_set_size_ordered, 
                                                 axis_param = list(gp = gpar(fontsize = 6),
                                                                   side = "top",
                                                                   labels_rot = 45,
                                                                   at = c(0,100,200,300),
                                                                   labels = c("0", "100", "200", ">300"))
                                                 ), 
                         annotation_name_gp = gpar(fontsize = 6),
                         annotation_name_side = "top")


overlap_hm <- ComplexHeatmap::Heatmap(name = "Overlap\nCoefficient", matrix = overlap_res_heatmap,
                                      col = colorRampPalette(c("white","purple1"))(35), 
                                      border = T, rect_gp = gpar(col = "white"), 
                                      show_column_names = F,
                                      row_names_gp = gpar(fontsize = 6), 
                                      heatmap_legend_param = list(labels_gp =  gpar(fontsize = 6),
                                                                  title_gp = gpar(fontsize = 6), 
                                                                  title_position = "leftcenter",
                                                                  gp = gpar(fontsize = 6),
                                                                  border = "black", 
                                                                  direction = "horizontal",
                                                                  legend_height = unit(0.1, "cm"),
                                                                  legend_width = unit(2.5, "cm")), 
                                      
                                      column_names_gp = gpar(fontsize = 6), 
                                      column_title_gp = gpar(fontsize = 6),
                                      cluster_rows = FALSE, 
                                      right_annotation = row_bar, 
                                      cell_fun = function(j, i, x, y, width, height, fill) {
                                        grid.text(sprintf("%.2f", overlap_res_heatmap[i, j]), x, y, gp = gpar(fontsize = 6))},
                                      width = 2, height = 2
)

png(file="overlap_adapative_signature_heatmap.png",
    width = 2.3, height = 3, unit = "in", res = 900)
draw(overlap_hm, 
     heatmap_legend_side = "bottom")
dev.off()

# ========================================================================================================
# Comparing CUDC/Osimert with proliferation
# ========================================================================================================

source("calculate_aggreg_expression.R")
gene_sets <- readRDS("pilot_screen_EGFRi_intersection_sigantures.RDS")

cudc_osimert_sig <- gene_sets[["CUDC/Osimert"]]

colData(cds)$cudc_osimert_sig <- calculate_aggregate_expression_score(cds, 
                                                                  signature_genes = cudc_osimert_sig)

ggplot(colData(cds) %>% as.data.frame() %>% filter(dose == 10000), aes(x = proliferation_index, y = cudc_osimert_sig)) +
  geom_point(aes(color = treatment, group = treatment)) +
  facet_grid(Cell.Line~treatment) +
  monocle3:::monocle_theme_opts()

cudc_osimert_sig_prolif_cor_df <- colData(cds) %>% as.dataframe() %>%
  select(Cell.Line, treatment, dose, cudc_osimert_sig, proliferation_index) %>%
  filter(dose %in% c(0,10000))

cor_result <- data.frame()
for (line in c("A172", "BT333", "T98G", "U87MG")) {
  for (drug in unique(cudc_osimert_sig_prolif_cor_df$treatment)) {
    temp_cor_df <- cudc_osimert_sig_prolif_cor_df %>%
      filter(Cell.Line == line & treatment == drug) %>%
      select(cudc_osimert_sig, proliferation_index)
    
    cor_result_temp <- cor.test(temp_cor_df$cudc_osimert_sig, 
                                temp_cor_df$proliferation_index, 
                                method = "spearman")
    cor_result_temp_row <- data.frame("cell_line" = line,
                                      "drug" = drug,
                                      "spearman_coeff" = cor_result_temp[["estimate"]]["rho"],
                                      "p_value" = cor_result_temp[["p.value"]], row.names = NULL)
    cor_result <- bind_rows(cor_result, cor_result_temp_row)
  }
}

cor_result_hm <- cor_result %>%
  select(-p_value) %>%
  pivot_wider(id_cols = cell_line, names_from = drug, values_from = spearman_coeff) %>%
  column_to_rownames("cell_line") %>%
  select("DMSO", "Afatinib", "Brigatinib", "CUDC-101", "EAI045", "Neratinib", "Osimertinib") %>%
  as.matrix()


hmcols <- circlize::colorRamp2(breaks = seq(-1,1,0.1), 
                               colors = colorRampPalette(c("darkblue","white","darkred"))(21))
cudc_osimert_sig_prolif_cor_hm <- ComplexHeatmap::Heatmap(name = "Spearman\nCoefficient", 
                                                          matrix = cor_result_hm,
                                      col = hmcols, 
                                      border = T, rect_gp = gpar(col = "white"), 
                                      column_names_gp = gpar(fontsize = 6),
                                      row_names_gp = gpar(fontsize = 6), 
                                      heatmap_legend_param = list(labels_gp =  gpar(fontsize = 6),
                                                                  title_gp = gpar(fontsize = 6), 
                                                                  title_position = "topleft",
                                                                  gp = gpar(fontsize = 6),
                                                                  border = "black", 
                                                                  # direction = "horizontal",
                                                                  legend_height = unit(1, "cm"),
                                                                  legend_width = unit(0.1, "cm"),
                                                                  at = c(-1,-0.5, 0, 0.5, 1),
                                                                  labels = c("-1", "-0.5", "0", "0.5", "1")), 
                                      
                                      column_title_gp = gpar(fontsize = 6),
                                      cluster_rows = FALSE, cluster_columns = FALSE,
                                      # right_annotation = row_bar, 
                                      cell_fun = function(j, i, x, y, width, height, fill) {
                                        grid.text(sprintf("%.2f", cor_result_hm[i, j]), x, y, gp = gpar(fontsize = 6))},
                                      width = 2, height = 2
)

cudc_osimert_sig_prolif_cor_hm

png(file="cudc_osim_sig_proliferation_cor_heatmap.png",
    width = 3, height = 1.2, unit = "in", res = 900)
draw(cudc_osimert_sig_prolif_cor_hm, 
     heatmap_legend_side = "right")
dev.off()

ggplot(cudc_osimert_sig_prolif_cor_df %>%
         filter(Cell.Line == "BT333" & treatment == "Osimertinib"), 
       aes(x = proliferation_index, y = cudc_osimert_sig)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", linewidth = 0.1, color = "blue", se = F) +
  scale_x_continuous(breaks = seq(0,3, 0.5)) +
  xlab("Proliferation Index") +
  ylab("CUDC/Osimert\nAgg Expression") +
  annotate("text", x = 2.5, y = 2.75, label = "rho = -0.05", size = 1.75) +
  theme(text = element_text(size = 6),
        axis.line = element_line(linewidth = 0.05),
        axis.ticks.length = unit(0.05, "mm")) +
  monocle3:::monocle_theme_opts()
ggsave("cudc_osim_sig_proliferation_cor_BT333_osimert_example.png",
       dpi = 900, width = 1.75, height = 1)

# pseudobulked because of sparsity
cudc_osimert_sig_pseudobulk <- colData(cds) %>% as.data.frame() %>% filter(dose %in% c(0, 10000)) %>%
  group_by(Cell.Line, treatment, well_ID) %>%
  summarise(mean_cudc_osimert_sig = mean(cudc_osimert_sig),
            mean_proliferation_index = mean(proliferation_index))

ggplot(cudc_osimert_sig_pseudobulk, aes(x = mean_proliferation_index,
                                        y = mean_cudc_osimert_sig)) +
  geom_point(aes(color = treatment, group = treatment)) +
  facet_grid(~ treatment) +
  monocle3:::monocle_theme_opts()

cor_result <- data.frame()
for (drug in unique(cudc_osimert_sig_pseudobulk$treatment)) {
  temp_cor_df <- cudc_osimert_sig_pseudobulk %>%
    filter(treatment == drug) %>%
    select(mean_cudc_osimert_sig, mean_proliferation_index)
  
  cor_result_temp <- cor.test(temp_cor_df$mean_cudc_osimert_sig, 
                              temp_cor_df$mean_proliferation_index, 
                              method = "pearson")
  cor_result_temp_row <- data.frame("drug" = drug,
                                    "spearman_coeff" = cor_result_temp[["estimate"]]["cor"],
                                    "p_value" = cor_result_temp[["p.value"]], row.names = NULL)
  cor_result <- bind_rows(cor_result, cor_result_temp_row)
}

temp_cor_df <- cudc_osimert_sig_pseudobulk %>%
  select(mean_cudc_osimert_sig, mean_proliferation_index)

cor_result_temp <- cor.test(temp_cor_df$mean_cudc_osimert_sig, 
                            temp_cor_df$mean_proliferation_index, 
                            method = "pearson")
cor_result_temp_row <- data.frame("drug" = drug,
                                  "spearman_coeff" = cor_result_temp[["estimate"]]["cor"],
                                  "p_value" = cor_result_temp[["p.value"]], row.names = NULL)

cor_result_hm <- cor_result %>%
  select(-p_value) %>%
  pivot_wider(names_from = drug, values_from = spearman_coeff) %>%
  select("DMSO", "Afatinib", "Brigatinib", "CUDC-101", "EAI045", "Neratinib", "Osimertinib") %>%
  mutate("Overall" = -0.3393) %>%
  as.matrix()


hmcols <- circlize::colorRamp2(breaks = seq(-1,1,0.1), 
                               colors = colorRampPalette(c("darkblue","white","darkred"))(21))
cudc_osimert_sig_prolif_cor_hm <- ComplexHeatmap::Heatmap(name = "Pearson\nCoefficient", 
                                                          matrix = cor_result_hm,
                                                          col = hmcols, 
                                                          border = T, rect_gp = gpar(col = "white"), 
                                                          column_names_gp = gpar(fontsize = 6),
                                                          row_names_gp = gpar(fontsize = 6), 
                                                          heatmap_legend_param = list(labels_gp =  gpar(fontsize = 6),
                                                                                      title_gp = gpar(fontsize = 6), 
                                                                                      title_position = "topleft",
                                                                                      gp = gpar(fontsize = 6),
                                                                                      border = "black", 
                                                                                      # direction = "horizontal",
                                                                                      legend_height = unit(1, "cm"),
                                                                                      legend_width = unit(0.1, "cm"),
                                                                                      at = c(-1,-0.5, 0, 0.5, 1),
                                                                                      labels = c("-1", "-0.5", "0", "0.5", "1")), 
                                                          
                                                          column_title_gp = gpar(fontsize = 6),
                                                          cluster_rows = FALSE, cluster_columns = FALSE,
                                                          # right_annotation = row_bar, 
                                                          cell_fun = function(j, i, x, y, width, height, fill) {
                                                            grid.text(sprintf("%.2f", cor_result_hm[i, j]), x, y, gp = gpar(fontsize = 6))},
                                                          width = 2, height = 2
)

cudc_osimert_sig_prolif_cor_hm

png(file="cudc_osim_sig_proliferation_cor_heatmap_pseudobulk.png",
    width = 3, height = 1, unit = "in", res = 900)
draw(cudc_osimert_sig_prolif_cor_hm, 
     heatmap_legend_side = "right")
dev.off()

ggplot(cudc_osimert_sig_pseudobulk %>% filter(treatment %in% c("CUDC-101", "Osimertinib")), 
       aes(x = mean_proliferation_index,y = mean_cudc_osimert_sig)) +
  geom_point(aes(shape = Cell.Line), size = 0.5) +
  geom_smooth(method = "lm", linewidth = 0.1, color = "blue", se = F) +
  facet_grid(~ treatment) + 
  xlab("Mean Proliferation Index") +
  ylab("Mean CUDC/Osimert\nAgg Expression") +
  scale_shape_manual(values = c(1,2,6,5)) +
  theme(text = element_text(size = 6),
        axis.line = element_line(linewidth = 0.05),
        axis.ticks.length = unit(0.05, "mm"),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.2, "cm"),
        legend.box.spacing = unit(c(0,-0.35,0,0), "mm"),
        legend.position = "bottom") +
  guides(shape = guide_legend(title = "Cell Line", override.aes = list(size = 2), byrow = TRUE, nrow = 1)) +
  monocle3:::monocle_theme_opts()
ggsave("cudc_osim_sig_proliferation_cor_BT333_osimert_example_pseudobulk.png",
       dpi = 900, width = 2.75, height = 1.5)
