library(tidyverse)
library(monocle3)
library(ComplexHeatmap)

cds <- readRDS("cds_large_screen.RDS")

# ================================================================================
# sci-Plex-GxE adaptive resistance signature
# ================================================================================
adaptive_up_genes <- read_csv("adaptive_program.csv",
                              col_names = TRUE) %>%
  filter(grepl(pattern = "pregulated", gene_module))

colData(cds)$adaptive_up <- calculate_aggregate_expression_score(cds, 
                                                                 signature_genes = adaptive_up_genes$gene_short_name,
                                                                 from_id = FALSE)

# Heatmaps of signature by dose
mean_adaptive_up <- colData(cds) %>% 
  as.data.frame() %>% 
  # group_by(cell_line) %>%
  # mutate(adaptive_up = scale(adaptive_up)) %>%
  group_by(cell_line, drug, dose) %>%
  summarise(mean_adaptive = mean(adaptive_up)) %>%
  group_by(cell_line) %>%
  mutate(mean_adaptive = scale(mean_adaptive)) %>%
  mutate(dose = case_when(
    drug == "Panitumumab" ~ as.character(log10(dose * 10)),
    drug %in% c("DMSO", "PBS", "Media") ~ drug,
    TRUE ~ as.character(log10(dose))
  )) %>%
  unite(col = "cell_line_drug",cell_line, drug, sep = "_") %>%
  pivot_wider(names_from = "cell_line_drug", values_from = mean_adaptive, id_cols = dose) %>%
  column_to_rownames("dose") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "cell_line_drug") %>%
  separate(cell_line_drug, into = c("cell_line", "drug"), remove = FALSE, sep = "_")

colha_df <- mean_adaptive_up %>%
  filter(!drug %in% c("DMSO", "PBS", "Media")) %>%
  select(cell_line_drug, cell_line, drug) %>%
  left_join(read_csv("Annotations_final_RG.csv", col_names = TRUE), by = c("drug" = "drug")) %>%
  # mutate(`Select Drug` = case_when(
  #   drug == "Osimertinib" ~ "Osimertinib",
  #   TRUE ~ NA
  # )) %>%
  select(cell_line_drug, cell_line, reversible) %>%
  column_to_rownames(var = "cell_line_drug")

mean_adaptive_up_controls <- mean_adaptive_up %>%
  filter(drug %in% c("DMSO", "PBS", "Media"))

mean_adaptive_up <- mean_adaptive_up %>%
  mutate(DMSO = case_when(
    cell_line == "BT112" ~ mean_adaptive_up_controls %>% filter(cell_line_drug == "BT112_DMSO") %>% select(DMSO) %>% pull(),
    cell_line == "BT228" ~ mean_adaptive_up_controls %>% filter(cell_line_drug == "BT228_DMSO") %>% select(DMSO) %>% pull(),
    cell_line == "BT333" ~ mean_adaptive_up_controls %>% filter(cell_line_drug == "BT333_DMSO") %>% select(DMSO) %>% pull()
  ),
  PBS = case_when(
    cell_line == "BT112" ~ mean_adaptive_up_controls %>% filter(cell_line_drug == "BT112_PBS") %>% select(PBS) %>% pull(),
    cell_line == "BT228" ~ mean_adaptive_up_controls %>% filter(cell_line_drug == "BT228_PBS") %>% select(PBS) %>% pull(),
    cell_line == "BT333" ~ mean_adaptive_up_controls %>% filter(cell_line_drug == "BT333_PBS") %>% select(PBS) %>% pull()),
  Media = case_when(
    cell_line == "BT112" ~ mean_adaptive_up_controls %>% filter(cell_line_drug == "BT112_Media") %>% select(Media) %>% pull(),
    cell_line == "BT228" ~ mean_adaptive_up_controls %>% filter(cell_line_drug == "BT228_Media") %>% select(Media) %>% pull(),
    cell_line == "BT333" ~ mean_adaptive_up_controls %>% filter(cell_line_drug == "BT333_Media") %>% select(Media) %>% pull())) %>%
  filter(!drug %in% c("DMSO", "PBS", "Media")) %>%
  select(-cell_line, -drug) %>%
  pivot_longer(cols = c("1","2","3","4", "PBS", "Media"), 
               names_to = "dose", 
               values_to = "agg_score") %>%
  # mutate(norm_agg_score = case_when(
  #   is.na(agg_score) ~ NA,
  #   TRUE ~ log2((agg_score + 1)/(DMSO + 1))
  # )) %>%
  # mutate(DMSO = 0) %>%
  mutate(norm_agg_score = case_when(
    is.na(agg_score) ~ NA,
    TRUE ~ agg_score
  )) %>%
  select(-agg_score) %>%
  pivot_wider(id_cols = c(cell_line_drug, DMSO), 
              names_from = dose, 
              values_from = norm_agg_score) %>%
  column_to_rownames(var = "cell_line_drug") %>%
  t()

colha <- ComplexHeatmap::HeatmapAnnotation(df = colha_df %>% select(`Cell Line` = cell_line, Reversible = reversible),
                                           col = list(`Cell Line` = setNames(RColorBrewer::brewer.pal(name = "Accent", n = 3), 
                                                                             c("BT112", "BT228", "BT333")),
                                                      Reversible = c("No" = "Red", "Yes" = "Blue")),
                                           # annotation_name_gp= gpar(fontsize = 12),
                                           annotation_legend_param = list(title_gp = gpar(fontface = "bold", fontsize = 6),
                                                                          labels_gp = gpar(fontsize = 6)))

ComplexHeatmap::Heatmap(matrix = mean_adaptive_up[c("DMSO","Media","PBS","1","2","3","4"),], 
                        name = "Adaptive Resistance Score",
                        # heatmap_legend_param = list(title = expression(bold(paste(Log[2], " (",frac(Prolif["cell,drug,dose"] + 1, Prolif["cell,DMSO"] + 1), ")   "))),
                        #                             title_gp = gpar(fontface = "bold", fontsize = 6),
                        #                             labels_gp = gpar(fontsize = 6)),
                        # col = hmcols,
                        show_column_names = TRUE,
                        column_names_gp = gpar(fontsize = 6),
                        column_split = 6,
                        cluster_columns = TRUE, 
                        cluster_rows = FALSE,
                        row_title = "Dose")
# top_annotation = colha)

saveRDS(mean_adaptive_up, 
        "adaptive_resistance_upgenes_zscored.RDS")


# ================================================================================
# Proliferation Index
# ================================================================================

mean_prolif <- colData(cds) %>% 
  as.data.frame() %>% 
  # group_by(cell_line, drug, dose) %>%
  # summarise(mean_prolif = mean(proliferation_index)) %>%
  group_by(cell_line, drug, dose) %>%
  summarise(mean_prolif = mean(proliferation_index)) %>%
  group_by(cell_line) %>%
  mutate(mean_prolif = scale(mean_prolif)) %>%
  mutate(dose = case_when(
    drug == "Panitumumab" ~ as.character(log10(dose * 10)),
    drug %in% c("DMSO", "PBS", "Media") ~ drug,
    TRUE ~ as.character(log10(dose))
  )) %>%
  unite(col = "cell_line_drug",cell_line, drug, sep = "_") %>%
  pivot_wider(names_from = "cell_line_drug", values_from = mean_prolif, id_cols = dose) %>%
  column_to_rownames("dose") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "cell_line_drug") %>%
  separate(cell_line_drug, into = c("cell_line", "drug"), remove = FALSE, sep = "_")

colha_df <- mean_prolif %>%
  filter(!drug %in% c("DMSO", "PBS", "Media")) %>%
  select(cell_line_drug, cell_line, drug) %>%
  left_join(read_csv("Annotations_final_RG.csv", col_names = TRUE), by = c("drug" = "drug")) %>%
  mutate(`Select Drug` = case_when(
    drug == "Osimertinib" ~ "Osimertinib", # to highlight an agent
    TRUE ~ NA
  )) %>%
  select(cell_line_drug, cell_line, reversible, `Select Drug`) %>%
  column_to_rownames(var = "cell_line_drug")

mean_prolif_controls <- mean_prolif %>%
  filter(drug %in% c("DMSO", "PBS", "Media"))

mean_prolif <- mean_prolif %>%
  mutate(DMSO = case_when(
    cell_line == "BT112" ~ mean_prolif_controls %>% filter(cell_line_drug == "BT112_DMSO") %>% select(DMSO) %>% pull(),
    cell_line == "BT228" ~ mean_prolif_controls %>% filter(cell_line_drug == "BT228_DMSO") %>% select(DMSO) %>% pull(),
    cell_line == "BT333" ~ mean_prolif_controls %>% filter(cell_line_drug == "BT333_DMSO") %>% select(DMSO) %>% pull()
  ),
  PBS = case_when(
    cell_line == "BT112" ~ mean_prolif_controls %>% filter(cell_line_drug == "BT112_PBS") %>% select(PBS) %>% pull(),
    cell_line == "BT228" ~ mean_prolif_controls %>% filter(cell_line_drug == "BT228_PBS") %>% select(PBS) %>% pull(),
    cell_line == "BT333" ~ mean_prolif_controls %>% filter(cell_line_drug == "BT333_PBS") %>% select(PBS) %>% pull()),
  Media = case_when(
    cell_line == "BT112" ~ mean_prolif_controls %>% filter(cell_line_drug == "BT112_Media") %>% select(Media) %>% pull(),
    cell_line == "BT228" ~ mean_prolif_controls %>% filter(cell_line_drug == "BT228_Media") %>% select(Media) %>% pull(),
    cell_line == "BT333" ~ mean_prolif_controls %>% filter(cell_line_drug == "BT333_Media") %>% select(Media) %>% pull())) %>%
  filter(!drug %in% c("DMSO", "PBS", "Media")) %>%
  select(-cell_line, -drug) %>%
  pivot_longer(cols = c("1","2","3","4", "PBS", "Media"), names_to = "dose", values_to = "prolif") %>%
  # mutate(norm_prolif = case_when(
  #   is.na(prolif) ~ NA,
  #   TRUE ~ log2((prolif + 1)/(DMSO + 1))
  # )) %>%
  # mutate(DMSO = 0) %>%
  mutate(norm_prolif = case_when(
    is.na(prolif) ~ NA,
    TRUE ~ prolif
  )) %>%
  select(-prolif) %>%
  pivot_wider(id_cols = c(cell_line_drug, DMSO), names_from = dose, values_from = norm_prolif) %>%
  column_to_rownames(var = "cell_line_drug") %>%
  t()

colha <- ComplexHeatmap::HeatmapAnnotation(df = colha_df %>% select(`Cell Line` = cell_line, Reversible = reversible, `Select Drug`),
                                           col = list(`Cell Line` = setNames(RColorBrewer::brewer.pal(name = "Accent", n = 3), 
                                                                             c("BT112", "BT228", "BT333")),
                                                      Reversible = c("No" = "Red", "Yes" = "Blue"),
                                                      `Select Drug` = setNames(viridis::magma(n = 4)[2], colha_df %>% select(`Select Drug`) %>% drop_na() %>% unique() %>% pull())),
                                           # annotation_name_gp= gpar(fontsize = 12),
                                           annotation_legend_param = list(title_gp = gpar(fontface = "bold", fontsize = 6),
                                                                          labels_gp = gpar(fontsize = 6)))
hmcols = colorRampPalette(viridis::plasma(n=10))(50)
ComplexHeatmap::Heatmap(matrix = mean_prolif[c("DMSO","Media","PBS","1","2","3","4"),], 
                        name = "prolif\nNorm to DMSO",
                        col = hmcols,
                        heatmap_legend_param = list(title = "Scaled Mean Expression",
                                                    title_gp = gpar(fontface = "bold", fontsize = 6),
                                                    labels_gp = gpar(fontsize = 6)),
                        # col = hmcols,
                        show_column_names = FALSE,
                        column_split = 4,
                        cluster_columns = TRUE, 
                        cluster_rows = FALSE,
                        row_title = "Dose",)
                        # top_annotation = colha)
dir.create("proliferation_heatmaps")
saveRDS(mean_prolif[c("DMSO","Media","PBS","1","2","3","4"),] %>% t(), 
        file = "proliferation_index_transposed_heatmap_data.RDS")


# ================================================================================
# correlation between MKI67 beta and number of DEGs
# ================================================================================
treatment_diff_test_results <- read_csv("supplementary_table_S4_large_EGFRi_screen_dose_DEG_fit_models_signif_result.csv")

sig_genes_results <- treatment_diff_test_results %>%
  filter(term %in% c("log(dose + 0.1)"),
         q_value < 0.01, 
         abs(normalized_effect) > 0.05)

DEGs_per_PDCL_compound <- sig_genes_results %>%
  group_by(cell_line, drug) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(log10_count = log10(count))

MKI_per_PDCL_compound <- treatment_diff_test_results %>%
  filter(term %in% c("log(dose + 0.1)")) %>%
  filter(gene_short_name == "MKI67") %>%
  select(MKI67_norm_beta = normalized_effect, cell_line, drug)

corr_plot <- left_join(DEGs_per_PDCL_compound, MKI_per_PDCL_compound)

corr_final <- data.frame()
for (i in c("BT112", "BT228", "BT333")) {
  corr_plot_temp <- filter(corr_plot, cell_line == i)
  temp <- cor(x = corr_plot_temp$log10_count, corr_plot_temp$MKI67_norm_beta, method = "pearson")
  corr_final <- bind_rows(corr_final,
                          data.frame("cell_line" = i, "pearson_coeff" = temp))
}

for (i in c("BT112", "BT228", "BT333")) {
  corr_plot_temp <- filter(corr_plot, cell_line == i) %>%
    mutate(diff = log10_count - 1.5 - abs(MKI67_norm_beta)) %>% 
    arrange(desc(abs(diff))) %>%
    mutate(rank = row_number()) %>%
    mutate(drug_label = case_when(
      rank %in% c(1,2,3,4) ~ drug
    ))
  
  ggplot(data = corr_plot_temp,
         aes(x = log10_count, y = MKI67_norm_beta)) +
    # facet_wrap(~cell_line) +
    geom_point(size = 0.1) +
    geom_smooth(method = "lm", se = F, linetype = 2, linewidth = 0.3) +
    annotate("text", x = 1.75, y = -0.09, label = paste0("r = ", round(corr_final %>% 
                                                                         filter(cell_line == i) %>%
                                                                         pull(pearson_coeff), 2)),
             size = 1.5) +
    xlab("Log10 (# DEGs)") +
    theme(text = element_text(size = 6),
          axis.ticks = element_line(linewidth = 0.1),
          axis.title.y = element_blank(),
          axis.ticks.length = unit(0.3, "mm"),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    monocle3:::monocle_theme_opts()
  ggsave(filename = paste0("MKI_global_corr/MKI67_DEGs_corr_", i, ".png"),
         dpi = 600, width = 0.9, height = 0.75)
  
}

# ================================================================================
# APM machinery expression and enrichment by Transcriptional Drug Class (output of MrVI scripts)
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

APM_MHC_CII_genes <- c("CIITA", 
                       "HLA-DRA", 
                       "HLA-DRB1", 
                       "HLA-DRB5", 
                       "HLA-DRB6", 
                       "HLA-DQA1", 
                       "HLA-DQB1", 
                       "HLA-DPA1", 
                       "HLA-DMA", 
                       "HLA-DMB", 
                       "HLA-DOA")

colData(cds)$APM_MHC_CI_score <- calculate_aggregate_expression_score(cds, signature_genes = APM_MHC_CI_genes)
colData(cds)$APM_MHC_CII_score <- calculate_aggregate_expression_score(cds, signature_genes = APM_MHC_CII_genes)

transcript_drug_groups <- read_csv("supplementary_table_S5_MrVI_TDC_drug_membership.csv") %>%
  select(drug, drug_group = `Drug Cluster`)

temp_col_data <- colData(cds) %>% as.data.frame() %>%
  left_join(transcript_drug_groups, by = c("drug" = "drug"))
rownames(temp_col_data) <- temp_col_data$cell_ID

colData(cds) <- DataFrame(temp_col_data)

glm_fit_df <- data.frame()
for (line in c("BT112", "BT228", "BT333")) {
  glm_test <- glm(APM_MHC_CI_score ~ drug_group, 
                  family = "gaussian",
                  data = temp_col_data %>% 
                    filter(cell_line == line & dose %in% c(0, 10000)) %>%
                    mutate(drug_group = factor(drug_group, 
                                               levels = c(16, 1:15))))
  
  temp_rows <- coef(summary(glm_test)) %>% as.data.frame() %>% 
    rownames_to_column(var = "drug_group") %>%
    mutate(cell_line = line)
  
  glm_fit_df <- bind_rows(glm_fit_df, temp_rows)
}

glm_fit_df_corrected <- glm_fit_df %>%
  dplyr::rename(beta_coeff = Estimate, 
                p_value = `Pr(>|t|)`) %>%
  mutate(q_value = p.adjust(p_value, method = "BH")) %>%
  mutate(plot_text = case_when(
    drug_group == "drug_group15" ~ "TDC15",
    TRUE ~ NA
  ))

ggplot(glm_fit_df_corrected %>% filter(drug_group != "(Intercept)"),
       aes(x = beta_coeff, y = -log10(q_value), 
           color = q_value < 0.05 & abs(beta_coeff) > 0.05)) +
  geom_point(size = 0.4) +
  geom_vline(xintercept = 0.05, linetype = "dashed", size = 0.1) +
  geom_vline(xintercept = -0.05, linetype = "dashed", size = 0.1) +
  geom_hline(yintercept = 0.05, linetype = "dashed", size = 0.1) +
  ggrepel::geom_text_repel(inherit.aes = F, 
                           min.segment.length = 0.05, 
                           segment.size = 0.1,
                           aes(x = beta_coeff, 
                               y = -log10(q_value), 
                               label = plot_text), size = 1.75, vjust = 1.1, hjust = 1) +
  monocle3:::monocle_theme_opts() +
  facet_wrap(~factor(cell_line, levels = c("BT112", "BT228", "BT333")), 
             nrow = 3,
             scales = "free") +
  xlab("Beta Coefficient") +
  ylab("-Log10 FDR") +
  scale_color_manual("FDR < 0.1%", values = c("TRUE" = "firebrick3", "FALSE" = "black")) +
  theme(text = element_text(size = 6),
        strip.text = element_text(size = 8),
        legend.position = "bottom",
        legend.key.height = unit(0.1, "cm"),
        legend.key.width = unit(0.25, "cm"),
        legend.margin = margin(t = -3, r = 0, b = 0, l = 0, unit = "mm")) +
  # legend.spacing.y = unit(0.1, "cm")) +
  guides(color = guide_legend(title = "FDR < 1%", byrow = TRUE, nrow = 2, override.aes = list(size=2)))
ggsave("APM_machinery_MCH-CI_enrichment_by_TDC.png",
       dpi = 600,
       height = 3.5, width = 1.2)

# MCH-CII
glm_fit_df <- data.frame()
for (line in c("BT112", "BT228", "BT333")) {
  glm_test <- glm(APM_MHC_CII_score ~ drug_group, 
                  family = "gaussian",
                  data = temp_col_data %>% 
                    filter(cell_line == line & dose %in% c(0, 10000)) %>%
                    mutate(drug_group = factor(drug_group, 
                                               levels = c(16, 1:15))))
  
  temp_rows <- coef(summary(glm_test)) %>% as.data.frame() %>% 
    rownames_to_column(var = "drug_group") %>%
    mutate(cell_line = line)
  
  glm_fit_df <- bind_rows(glm_fit_df, temp_rows)
}

glm_fit_df_corrected <- glm_fit_df %>%
  dplyr::rename(beta_coeff = Estimate, 
                p_value = `Pr(>|t|)`) %>%
  mutate(q_value = p.adjust(p_value, method = "BH")) %>%
  mutate(plot_text = case_when(
    beta_coeff > 0.3 ~ str_remove(string = drug_group, pattern = "drug_group"),
    TRUE ~ NA
  ))

ggplot(glm_fit_df_corrected %>% filter(drug_group != "(Intercept)"),
       aes(x = beta_coeff, y = -log10(q_value), 
           color = q_value < 0.05 & abs(beta_coeff) > 0.05)) +
  geom_point(size = 0.4) +
  geom_vline(xintercept = 0.05, linetype = "dashed", size = 0.1) +
  geom_vline(xintercept = -0.05, linetype = "dashed", size = 0.1) +
  geom_hline(yintercept = 0.05, linetype = "dashed", size = 0.1) +
  ggrepel::geom_text_repel(inherit.aes = F, 
                           min.segment.length = 0.05, 
                           segment.size = 0.1,
                           aes(x = beta_coeff, 
                               y = -log10(q_value), 
                               label = plot_text), size = 1.75, vjust = 1.1, hjust = 1) +
  monocle3:::monocle_theme_opts() +
  facet_wrap(~factor(cell_line, levels = c("BT112", "BT228", "BT333")), 
             nrow = 3,
             scales = "free") +
  xlab("Beta Coefficient") +
  ylab("-Log10 FDR") +
  scale_color_manual("FDR < 0.1%", values = c("TRUE" = "firebrick3", "FALSE" = "black")) +
  theme(text = element_text(size = 6),
        strip.text = element_text(size = 8),
        legend.position = "bottom",
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.margin = margin(t = -3, r = 0, b = 0, l = 0, unit = "mm")) +
  # legend.spacing.y = unit(0.1, "cm")) +
  guides(color = guide_legend(title = "FDR < 1%", byrow = TRUE, override.aes = list(size=2)))
ggsave("APM_machinery_MCH-CII_enrichment_by_TDC.png",
       dpi = 600,
       height = 3, width = 1.5)




