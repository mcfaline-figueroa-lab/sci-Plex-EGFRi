library(tidyverse)
library(monocle3)
library(ComplexHeatmap)

this.dir <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(this.dir))
rm(this.dir)

#============================================================================================
# Lollipop collection of DRC signatures across drugs 
#============================================================================================

calculate_aggregate_expression_score <- function(cds, signature_genes, from_id = FALSE){
  
  if(from_id == TRUE){
    cds_subset = cds[rowData(cds)$id %in% signature_genes,]
  }
  else(cds_subset = cds[rowData(cds)$gene_short_name %in% signature_genes,])
  aggregate_signature_expression = exprs(cds_subset)
  aggregate_signature_expression = t(t(aggregate_signature_expression) / pData(cds_subset)$Size_Factor)
  aggregate_signature_expression = Matrix::colSums(aggregate_signature_expression)
  aggregate_signature_expression = log(aggregate_signature_expression+1)
  return(aggregate_signature_expression)
}

cds <- readRDS("cds_large_screen.RDS")

# Read in DEG lists for each DRC 
files <- (Sys.glob("DRC_up_signature_lfcANDpval_filtered_01252024/*"))
response_group_list <- list()
for (file in files) {
  temp_line <- as.data.frame("x" = file) %>%
    separate(file, into = c(NA, "file"), sep = "DRC_up_signature_lfcANDpval_filtered_01252024/") %>%
    separate(file, into = c("line", "signature", NA, NA), sep = "_") %>%
    mutate(line = toupper(line)) %>%
    unite(col = "line_signature", sep = "_", line, signature) %>%
    pull()
  
  response_group_list[[temp_line]] <- read_csv(file) %>%
    dplyr::select(genes = 2) %>%
    pull()
}
rm(file, files, temp_line)

col_info_merge <- colData(cds) %>% 
  as.data.frame() %>% 
  mutate(dose = case_when(
    drug == "Panitumumab" ~ dose * 10,
    TRUE ~ dose
  )) %>%
  select(cell_line, drug, dose) %>%
  unite(col = "id", sep = "_", cell_line, drug, dose)

# Calculate signatures from the DRC DEGs 
sig_mat <- data.frame()
for (name in names(response_group_list)) {
  line <- str_split_i(string = name, pattern = "_", i = 1)
  temp <- calculate_aggregate_expression_score(cds[,colData(cds)$cell_line == line],
                                               signature_genes = response_group_list[[name]],
                                               from_id = FALSE) %>% as.data.frame() %>% dplyr::rename(score = 1)
  temp$cell_ID <- row.names(temp)

  temp <- left_join(temp, 
                    col_info_merge %>% rownames_to_column(var = "cell_ID"), 
                    by = c("cell_ID" = "cell_ID")) %>%
    separate(id, 
             into = c("cell_line", "drug", "dose"), 
             sep = "_", 
             remove = F)

  temp_median <- temp %>% 
    group_by(id, cell_line, drug, dose) %>%
    summarise(med_score = mean(score), 
              .groups = "drop") %>%
    mutate(signature = name)
  
  temp_DMSO <- temp_median %>%
    ungroup() %>%
    filter(grepl("DMSO", id)) %>%
    select(cell_line, med_score_DMSO = med_score)
  
  if (temp_DMSO$med_score_DMSO[1] != 0) {
    
    temp_percent <- left_join(temp, temp_DMSO, by = c("cell_line" = "cell_line")) %>%
      mutate(above = case_when(
        score >= med_score_DMSO ~ 1,
        TRUE ~ 0
      )) %>%
      group_by(id, med_score_DMSO) %>%
      summarise(percent_over_DMSO = sum(above) / n(),
                .groups = "drop")
    
    temp_final <- left_join(temp_median, 
                            temp_percent, 
                            by = c("id" = "id")) %>%
      mutate(signature = name) %>%
      mutate(med_score_norm_DMSO = med_score/med_score_DMSO)
    sig_mat <- bind_rows(sig_mat, temp_final)
  }
}

# Create dataframe for plotting 
sig_mat_plot <- sig_mat %>% 
  mutate(id_plot = factor(id)) %>%
  mutate(rank = factor(drug))

# Insert TDC annotations (from MrVI binary membership script)
drug_cluster <- read.csv('SharedPDCL_MrVI_Drug_DRC_BinaryMembership_HierarchicalClusters_FDRfilt_TDCs.csv') %>%
  dplyr::rename(drug = 1, drug_cluster = 2)
sig_mat_plot_test <- left_join(sig_mat_plot, drug_cluster, by = c("drug" = "drug"))

# Order data
sig_mat_plot_test <- sig_mat_plot_test %>%
  filter(dose == "10000" | dose == "0") %>%
  mutate(id_plot = factor(signature)) %>%
  arrange(drug_cluster) %>%
  mutate(drug = fct_reorder(drug, drug_cluster)) %>%
  mutate(med_score_norm_DMSO_plot = log2(med_score_norm_DMSO)) %>%
  mutate(med_score_norm_DMSO_plot = case_when(
    med_score_norm_DMSO_plot >= 0.2 ~ 0.2,
    med_score_norm_DMSO_plot <= -0.2 ~ -0.2,
    TRUE ~ med_score_norm_DMSO_plot
  )) %>%
  mutate(line_length = percent_over_DMSO*8) %>%
  mutate(line_length = case_when(
    line_length <= 0.9 ~ 0.9,
    TRUE ~ line_length
  ))

# saveRDS(sig_mat_plot_test, "signature_matrix_for_plotting.RDS")

# Create plot
color_scheme_rand <- readRDS("random_color_scheme_transcriptional_drug_class.RDS")

ggplot(sig_mat_plot_test %>% arrange(drug), aes(x = id_plot, y = line_length)) +
  # facet_wrap(~drug, ncol = 12) + 
  facet_wrap(~drug, ncol = 8) +
  geom_segment(
    aes(x = id_plot, xend = id_plot, y = 0, yend = line_length, color = cell_line), 
    size = 0.2
  ) + 
  scale_color_manual(values = c("grey70", "grey50", "black")) +
  ggnewscale::new_scale_color() +
  geom_rect(
    aes(xmin = 0, xmax = 65, ymin = 0, ymax = 0.8, fill = as.character(drug_cluster)),
    show.legend = F
  ) +
  geom_point(aes(color = med_score_norm_DMSO_plot,
                 size = med_score_norm_DMSO_plot), 
             alpha = 0.5, stroke = 0.1) +
  # geom_text(aes(label = id_plot), size = 0.5) +
  scale_fill_manual(values = color_scheme_rand, 
                    breaks = seq(1,16,1)) +
  # viridis::scale_color_viridis(discrete = T, option = "D") +
  viridis::scale_color_viridis(discrete = F, option = "plasma", end = 0.9
  ) +
  scale_x_discrete(expand = c(0,0)) +
  scale_size_continuous(
    range = c(0.6,1.2)
  ) +
  theme_void() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.clip = "on",
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 6),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.justification = "center",
        legend.spacing.x = unit(0.5, "mm"),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.2, "cm"),
        panel.spacing.y = unit(0, "mm")) +
  guides(color= guide_legend(title = expression(bold(paste(Log[2], " (",frac("Mean Score"["10uM drug"], "Mean Score"["DMSO"]), ")   ")))),
         size= guide_legend(title = expression(bold(paste(Log[2], " (",frac("Mean Score"["10uM drug"], "Mean Score"["DMSO"]), ")   "))),
                            override.aes = list(size = c(2,2.5,3,3.5,4)))) + 
  coord_polar()

ggsave("largescreen_signature_by_drug_lollipop_UpGenes_normtoDMSO_scale_size_cell_line_long.png", 
       dpi = 900, height = 9, width = 6, bg = "white")

ggplot(sig_mat_plot_test %>% arrange(drug) %>% filter(drug == "AZ5104"), aes(x = id_plot, y = percent_over_DMSO)) +
  geom_segment(
    aes(x = id_plot, xend = id_plot, y = 0, yend = percent_over_DMSO, color = cell_line), 
    size = 0.4, show.legend = F
  ) + 
  scale_color_manual(values = c("grey70", "grey50", "black")) +
  ggnewscale::new_scale_color() +
  geom_point(aes(color = med_score_norm_DMSO_plot,
                 size = med_score_norm_DMSO_plot), 
             alpha = 0.5, stroke = 0.1, show.legend = F) +
  scale_fill_manual(values = color_scheme_rand, 
                    breaks = seq(1,16,1)) +
  viridis::scale_color_viridis(discrete = F, option = "plasma", end = 0.9
  ) +
  scale_size_continuous(
    range = c(0.6,1.2)
  ) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_blank(),
        text = element_text(size = 6),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.justification = "center",
        legend.spacing.x = unit(0.5, "mm"),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.2, "cm"), 
        panel.border = element_rect(fill = NULL, color = "white")) +
  guides(color= guide_legend(title = expression(bold(paste(Log[2], " (",frac("Mean Score"["drug"], "Mean Score"["DMSO"]), ")   ")))),
         size= guide_legend(title = expression(bold(paste(Log[2], " (",frac("Mean Score"["drug"], "Mean Score"["DMSO"]), ")   "))),
                            override.aes = list(size = c(2,2.5,3,3.5,4)))) +
  monocle3:::monocle_theme_opts() +
  coord_flip()
ggsave(paste0("larger_groups/AZ5104_legend_linear.png"), 
       dpi = 900, height = 2.5, width = 0.5, bg = "white")

# ===========================================================================================
# correlating DRC scores for DG4 and DG5 (osimertinib vs AZ5104)
# ===========================================================================================
sig_mat_plot_test <- readRDS("signature_matrix_for_plotting.RDS")

cudc_print <- sig_mat_plot_test %>% filter(drug == "CUDC-101") %>% arrange(id_plot)
write_csv(cudc_print, "CUDC_DRC_lollipop_scores.csv")

AG555_print <- sig_mat_plot_test %>% filter(drug == "AG555") %>% arrange(id_plot)
write_csv(AG555_print, "AG555_DRC_lollipop_scores.csv")

sig_corr <- sig_mat_plot_test %>%
  filter(drug %in% c("Osimertinib", "HS-10296")) %>%
  select(drug, signature, med_score_norm_DMSO) %>%
  pivot_wider(id_cols = signature, values_from = med_score_norm_DMSO, names_from = drug) %>%
  column_to_rownames(var = "signature")
cor_test_res <- cor(x = sig_corr, method = "pearson")

sig_corr_plot <- sig_corr %>% 
  rownames_to_column(var = "signature") %>%
  mutate(diff = `HS-10296` - Osimertinib) %>% 
  arrange(desc(abs(diff))) %>%
  mutate(rank = row_number()) %>%
  mutate(drc_label = case_when(
    rank %in% c(1,2,3,4) ~ signature
  )) %>%
  mutate(drc_label = str_replace_all(drc_label, pattern = "DRC", "RM"))

ggplot(sig_corr_plot, aes(x = Osimertinib, y = `HS-10296`)) +
  geom_point(alpha = 0.5, stroke = 0.01, size = 1) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = 2, linewidth = 0.1) +
  ggrepel::geom_text_repel(aes(label = drc_label), size = 2, 
                           min.segment.length = 0.05, 
                           segment.size = 0.2, max.overlaps = 3, force_pull = 3) +
  scale_y_continuous(breaks = seq(0.3,1.8,0.3)) +
  scale_x_continuous(breaks = seq(0.3,1.8,0.3)) +
  annotate(x = 0.4, y = 1.6, 
           geom = "text", 
           label = paste0("r = ", round(cor_test_res[1,2], 3)),
           size = 2) +
  theme(text = element_text(size = 6),
        axis.ticks = element_line(linewidth = unit(0.2, "mm"))) +
  monocle3:::monocle_theme_opts()
ggsave("osimertinib_AZ5104_compare/osimertinib_HS-10296_DRC_correlation.png",
       dpi = 600, height = 2, width = 2)

sig_corr <- sig_mat_plot_test %>%
  filter(drug %in% c("Osimertinib", "AZ5104")) %>%
  select(drug, signature, med_score_norm_DMSO) %>%
  pivot_wider(id_cols = signature, values_from = med_score_norm_DMSO, names_from = drug) %>%
  column_to_rownames(var = "signature")
cor_test_res <- cor(x = sig_corr, method = "pearson")

sig_corr_plot <- sig_corr %>% 
  rownames_to_column(var = "signature") %>%
  mutate(diff = AZ5104 - Osimertinib) %>% 
  arrange(desc(abs(diff))) %>%
  mutate(rank = row_number()) %>%
  mutate(drc_label = case_when(
    rank %in% c(1,2,3,4) ~ signature
  )) %>%
  mutate(drc_label = str_replace_all(drc_label, pattern = "DRC", "RM"))

ggplot(sig_corr_plot, aes(x = Osimertinib, y = AZ5104)) +
  geom_point(alpha = 0.5, stroke = 0.01, size = 1) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = 2, linewidth = 0.1) +
  ggrepel::geom_text_repel(aes(label = drc_label), size = 2, 
                           min.segment.length = 0.05, 
                           segment.size = 0.2, max.overlaps = 2, force_pull = 3) +
  scale_y_continuous(breaks = seq(0.3,1.8,0.3)) +
  scale_x_continuous(breaks = seq(0.3,1.8,0.3)) +
  annotate(x = 0.4, y = 1.6, 
           geom = "text", 
           label = paste0("r = ", round(cor_test_res[1,2], 3)),
           size = 2) +
  theme(text = element_text(size = 6),
        axis.ticks = element_line(linewidth = unit(0.2, "mm"))) +
  monocle3:::monocle_theme_opts()
ggsave("osimertinib_AZ5104_compare/osimertinib_AZ5104_DRC_correlation.png",
       dpi = 600, height = 2, width = 2)

sig_corr <- sig_mat_plot_test %>%
  filter(drug %in% c("Osimertinib", "Puromycin")) %>%
  select(drug, signature, med_score_norm_DMSO) %>%
  pivot_wider(id_cols = signature, values_from = med_score_norm_DMSO, names_from = drug) %>%
  column_to_rownames(var = "signature")
cor_test_res <- cor(x = sig_corr, method = "pearson")

sig_corr_plot <- sig_corr %>% 
  rownames_to_column(var = "signature") %>%
  mutate(diff = Puromycin - Osimertinib) %>% 
  arrange(desc(abs(diff))) %>%
  mutate(rank = row_number()) %>%
  mutate(drc_label = case_when(
    rank %in% c(1,2,3,4) ~ signature
  )) %>%
  mutate(drc_label = str_replace_all(drc_label, pattern = "DRC", "RM"))

ggplot(sig_corr_plot, aes(x = Osimertinib, y = Puromycin)) +
  geom_point(alpha = 0.5, stroke = 0.01, size = 1) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = 2, linewidth = 0.1) +
  # ggrepel::geom_text_repel(aes(label = drc_label), size = 1, 
  #                          min.segment.length = 0.05, 
  #                          segment.size = 0.2, max.overlaps = 2, force_pull = 3) +
  scale_y_continuous(breaks = seq(0.3,1.8,0.3)) +
  scale_x_continuous(breaks = seq(0.3,1.8,0.3)) +
  annotate(x = 0.4, y = 1.6, 
           geom = "text", 
           label = paste0("r = ", round(cor_test_res[1,2], 3)),
           size = 2) +
  theme(text = element_text(size = 6),
        axis.ticks = element_line(linewidth = unit(0.2, "mm"))) +
  monocle3:::monocle_theme_opts()
ggsave("osimertinib_AZ5104_compare/osimertinib_puromycin_DRC_correlation.png",
       dpi = 600, height = 2, width = 2)

sig_corr_final <- data.frame()
for (i in unique(sig_mat_plot_test$drug)) {
  for (j in unique(sig_mat_plot_test$drug)) {
    sig_corr_df_temp <- sig_mat_plot_test %>%
      filter(drug %in% c(i, j)) %>%
      select(drug, signature, med_score_norm_DMSO) %>%
      pivot_wider(id_cols = signature, values_from = med_score_norm_DMSO, names_from = drug) %>%
      column_to_rownames(var = "signature")
    cor_test_temp <- cor(x = , method = "pearson")
  }
}

sig_corr_df_temp <- sig_mat_plot_test %>%
  select(drug, signature, med_score_norm_DMSO) %>%
  pivot_wider(id_cols = signature, values_from = med_score_norm_DMSO, names_from = drug) %>%
  column_to_rownames(var = "signature")
cor_test_temp <- cor(sig_corr_df_temp, method = "pearson")


hmcols <- circlize::colorRamp2(breaks = seq(-1,1,0.1), 
                               colors = colorRampPalette(c("darkblue","white","darkred"))(21))
sig_cor_hm <- ComplexHeatmap::Heatmap(name = "Pearson\nCoefficient", 
                                                          matrix = cor_test_temp,
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
                                                          # cell_fun = function(j, i, x, y, width, height, fill) {
                                                          #   grid.text(sprintf("%.2f", cor_result_hm[i, j]), x, y, gp = gpar(fontsize = 6))},
                                                          width = 2, height = 2
)

sig_cor_hm

png(file="osimertinib_AZ5104_compare/heatmap_RM_correlation.png",
    width = 6, height = 6, unit = "in", res = 900)
draw(sig_cor_hm, 
     heatmap_legend_side = "right")
dev.off()

