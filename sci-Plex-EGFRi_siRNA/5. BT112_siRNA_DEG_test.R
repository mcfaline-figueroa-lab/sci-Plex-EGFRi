library(tidyverse)
library(monocle3)
library(ComplexHeatmap)
this.dir <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(this.dir))
rm(this.dir)

cds <- readRDS("sci-Plex-EGFRi-siRNA-cds.RDS")

colData(cds)$siRNA_condition <- factor(colData(cds)$siRNA_condition, 
                                       levels = c("siNTC", "siEGFR1", "siEGFR2"))

expressed_genes <-  rowData(cds)[rowData(cds)$num_cells_expressed >= ncol(cds)*0.05,]$id

treatment_diff_test.list <- list()
for (time in c("24hr", "48hr", "72hr")) {
  treatment_diff_test.list[[time]] <- fit_models(cds = cds[expressed_genes,
                                                           colData(cds)$timepoint == time], 
                                                 model_formula_str = "~ siRNA_condition + replicate",
                                                 expression_family = "quasipoisson",
                                                 cores = 1, 
                                                 verbose = TRUE)
  treatment_diff_test.list[[time]] <- coefficient_table(treatment_diff_test.list[[time]]) %>% 
    dplyr::select(-model,-model_summary)
  treatment_diff_test.list[[time]]$siRNA_condition <- stringr::str_remove(treatment_diff_test.list[[time]]$term,
                                                                          "siRNA_condition")
  treatment_diff_test.list[[time]]$timepoint <- time
  
  message("Finished timepoint ", time)
}

saveRDS(treatment_diff_test.list, "BT112_siRNA_treatment_diff_test_list.RDS")

treatment_diff_test.list <- readRDS("BT112_siRNA_treatment_diff_test_list.RDS")
for (time in c("24hr", "48hr", "72hr")) {
  treatment_diff_test.list[[time]]$q_value <- p.adjust(treatment_diff_test.list[[time]]$p_value, 
                                                       method = "BH")
}

treatment_diff_test <- do.call("rbind", treatment_diff_test.list)

treatment_diff_test_siRNA <- treatment_diff_test %>%
  filter(grepl(pattern = "siRNA_condition", term)) %>%
  filter(q_value <= 0.05)

treatment_diff_test_siRNA_plot <- treatment_diff_test %>%
  filter(grepl(pattern = "siRNA_condition", term)) %>%
  mutate(plot_color = case_when(
    abs(normalized_effect) >= 0.05 & q_value <= 0.05 ~ "red",
    TRUE ~ "black"
  )) %>%
  mutate(siRNA_condition = factor(siRNA_condition, levels = c("siNTC", "siEGFR1", "siEGFR2"))) %>%
  mutate(plot_q_value = case_when(
    -log10(p_value) >= 30 ~ 30,
    TRUE ~ -log10(q_value)
  ))

library(ggrepel)
library(ggpp)
ggplot(treatment_diff_test_siRNA_plot %>% 
         mutate(plot_gene = case_when(
           gene_short_name == "EGFR" ~ gene_short_name,
           TRUE ~ NA
         )) %>%
         mutate(plot_color = case_when(
           !is.na(plot_gene) ~ "red",
           TRUE ~ "black"
         )) %>%
         arrange(plot_color) %>%
         mutate(plot_color = factor(plot_color)),
       aes(x = normalized_effect, y = plot_q_value)) +
  facet_grid(timepoint~siRNA_condition) +
  geom_point(aes(color = plot_color),
             size = 0.2,
             show.legend = F) +
  ggrepel::geom_text_repel(aes(label = plot_gene), 
                           size = 1.7,
                           min.segment.length = unit(0, 'lines'),
                           position = position_nudge_to(x = -1, y = 6),
                           segment.size = 0.1
  ) +
  geom_vline(xintercept = 0.1, linetype = 2, linewidth = 0.1) +
  geom_vline(xintercept = -0.1, linetype = 2, linewidth = 0.1) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.1) +
  scale_color_manual(values = c("black", "red")) +
  xlab("Normalized Effect") +
  ylab("-Log10 FDR") +
  theme(text = element_text(size = 6),
        axis.ticks.length = unit(0.5, "mm"),
        axis.ticks = element_line(linewidth = 0.1)) +
  monocle3:::monocle_theme_opts()
ggsave("BT112_siRNA_volcano_plot_EGFR_highlight.png",
       dpi = 600, height = 3, width = 2)

ggplot(treatment_diff_test_siRNA_plot %>% 
         mutate(plot_gene = case_when(
           q_value <= 0.05 ~ gene_short_name,
           TRUE ~ NA
         )),
       aes(x = normalized_effect, y = plot_q_value)) +
  facet_grid(timepoint~siRNA_condition) +
  geom_point(aes(color = plot_color),
             size = 0.2,
             show.legend = F) +
  ggrepel::geom_text_repel(aes(label = plot_gene), 
                           size = 1.3, min.segment.length = unit(0, "line"),
                           segment.size = 0.1, max.overlaps = 30
  ) +
  geom_vline(xintercept = 0.1, linetype = 2, linewidth = 0.1) +
  geom_vline(xintercept = -0.1, linetype = 2, linewidth = 0.1) +
  scale_color_manual(values = c("black", "red")) +
  xlab("Normalized Effect") +
  ylab("-Log10 FDR") +
  theme(text = element_text(size = 6),
        axis.ticks.length = unit(0.5, "mm"),
        axis.ticks = element_line(linewidth = 0.1)) +
  monocle3:::monocle_theme_opts()
ggsave("BT112_siRNA_volcano_plot_top_highlight.png",
       dpi = 600, height = 3, width = 2)

for (input_gene in c("ERAP1", "ERRFI1", "ERBB3", "CCNB1", "CCDC84")) {
  ggplot(treatment_diff_test_siRNA_plot %>% 
           mutate(plot_gene = case_when(
             gene_short_name == input_gene ~ gene_short_name,
             TRUE ~ NA
           )) %>%
           mutate(plot_color = case_when(
             !is.na(plot_gene) ~ "red",
             TRUE ~ "black"
           )) %>%
           arrange(plot_color) %>%
           mutate(plot_color = factor(plot_color)),
         aes(x = normalized_effect, y = plot_q_value)) +
    facet_wrap(~siRNA_condition) +
    geom_point(aes(color = plot_color),
               size = 0.2,
               show.legend = F) +
    ggrepel::geom_text_repel(aes(label = plot_gene), 
                             size = 1.7,
                             min.segment.length = unit(0, 'lines'),
                             segment.size = 0.1, 
                             position = position_nudge_to(x = -1, y = 5)
    ) +
    geom_vline(xintercept = 0.1, linetype = 2, linewidth = 0.1) +
    geom_vline(xintercept = -0.1, linetype = 2, linewidth = 0.1) +
    geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.1) +
    scale_color_manual(values = c("black", "red")) +
    xlab("Normalized Effect") +
    ylab("-Log10 FDR") +
    theme(text = element_text(size = 6),
          axis.ticks.length = unit(0.5, "mm"),
          axis.ticks = element_line(linewidth = 0.1)) +
    monocle3:::monocle_theme_opts()
  ggsave(paste0("BT112_siRNA_volcano_plot_",input_gene,"_highlight.png"),
         dpi = 600, height = 2.5, width = 3.5)
}


ggplot(treatment_diff_test_siRNA_plot %>% 
         mutate(plot_gene = case_when(
           gene_short_name %in% c("EGFR", "ERBB2", "ERBB3", "ERBB4") ~ gene_short_name,
           TRUE ~ NA
         )) %>%
         mutate(plot_color = case_when(
           !is.na(plot_gene) ~ "red",
           TRUE ~ "black"
         )) %>%
         arrange(plot_color) %>%
         mutate(plot_color = factor(plot_color)),
       aes(x = normalized_effect, y = plot_q_value)) +
  facet_wrap(~siRNA_condition) +
  geom_point(aes(color = plot_color),
             size = 0.2,
             show.legend = F) +
  ggrepel::geom_text_repel(aes(label = plot_gene), 
                           size = 1.7,
                           min.segment.length = unit(0, 'lines'),
                           segment.size = 0.1,
                           position = position_nudge_center(x = 0.5, center_x = 0,
                                                            y = 5, center_y = 0)
  ) +
  geom_vline(xintercept = 0.1, linetype = 2, linewidth = 0.1) +
  geom_vline(xintercept = -0.1, linetype = 2, linewidth = 0.1) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.1) +
  scale_color_manual(values = c("black", "red")) +
  xlab("Normalized Effect") +
  ylab("-Log10 FDR") +
  theme(text = element_text(size = 6),
        axis.ticks.length = unit(0.5, "mm"),
        axis.ticks = element_line(linewidth = 0.1)) +
  monocle3:::monocle_theme_opts()
ggsave("BT112_siRNA_volcano_plot_ERBB_highlight.png",
       dpi = 600, height = 2.5, width = 3.5)

ggplot(treatment_diff_test_siRNA_plot %>% 
         filter(timepoint == "48hr") %>%
         mutate(plot_gene = case_when(
           gene_short_name %in% c("EGFR", "BRAF", "PTEN", "SHC3", "SOS1",
                                  "IFI6", "NDRG1", "MAP3K1", "AURKB", "MKI67",
                                  "E2F1", "JUN", "EGR3") &
             q_value < 0.05 ~ gene_short_name,
           TRUE ~ NA
         )) %>%
         mutate(plot_color = case_when(
           !is.na(plot_gene) ~ "red",
           TRUE ~ "grey60"
         )) %>%
         arrange(plot_color) %>%
         mutate(plot_color = factor(plot_color)),
       aes(x = normalized_effect, y = plot_q_value)) +
  facet_wrap(~siRNA_condition, nrow = 1) +
  geom_point(aes(color = plot_color),
             size = 0.2,
             show.legend = F) +
  ggrepel::geom_label_repel(aes(label = plot_gene), 
                            size = 1.4, min.segment.length = unit(0, "line"),
                            segment.size = 0.3, max.overlaps = 15, 
                            force_pull = 0.05,
                            box.padding = 0.5, 
                            label.padding = 0.08
                            # segment.curvature = -0.1
                            # segment.ncp = 3,
                            # segment.angle = 20
  ) +
  geom_vline(xintercept = 0.1, linetype = 2, linewidth = 0.1) +
  geom_vline(xintercept = -0.1, linetype = 2, linewidth = 0.1) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.1) +
  scale_color_manual(values = c("grey60", "red")) +
  xlab("Normalized Effect") +
  ylab("-Log10 FDR") +
  theme(text = element_text(size = 7),
        axis.ticks.length = unit(0.5, "mm"),
        axis.ticks = element_line(linewidth = 0.1)) +
  monocle3:::monocle_theme_opts()
ggsave("BT112_48hr_siRNA_volcano_plot_select_gene_highlight.png",
       dpi = 600, height = 1.5, width = 2)

# Correlating siRNA DEGs within timepoints
for (time in c("24hr", "48hr", "72hr")) {
  cor_test_df <- treatment_diff_test_siRNA_plot %>%
    # filter(siRNA_condition == "siEGFR1") %>%
    filter(timepoint %in% time) %>%
    filter(id %in% treatment_diff_test_siRNA$id) %>%
    # select(id, normalized_effect, timepoint) %>%
    select(id, normalized_effect, siRNA_condition) %>%
    pivot_wider(names_from = siRNA_condition,
                values_from = normalized_effect,
                id_cols = id) %>%
    # pivot_wider(names_from = timepoint,
    #             values_from = normalized_effect,
    #             id_cols = id) %>%
    tibble::column_to_rownames(var = "id")
  
  cor_test <- cor(cor_test_df, method = "pearson")
  
  t <- ggplot(cor_test_df, aes(x = siEGFR1, y = siEGFR2)) +
    geom_point(size = 0.1) +
    geom_smooth(method = "lm", linewidth = 0.4, se = F) +
    annotate(geom = "text", x = 0.6, y = -1, 
             label = paste0("r = ",as.character(round(cor_test[1,2], 3))),
             size = 1.25) +
    ggtitle(paste0(time)) +
    theme(axis.ticks.length = unit(0.6, "mm"),
          axis.ticks = element_line(linewidth = 0.1),
          text = element_text(size = 6),
          plot.title = element_text(size = 6, hjust = 0.5)) +
    monocle3:::monocle_theme_opts()
  assign(paste0("plot_",time), t)

}

library(cowplot)
p <- cowplot::plot_grid(plotlist = list(plot_24hr, plot_48hr, plot_72hr), nrow = 1)

ggsave("BT112_siRNA_beta_coeff_cor_filt_union.png", plot = p, dpi = 600, height = 1.5, width = 3)


# Correlating DEGs between timepoints
for (siRNA in c("siEGFR1", "siEGFR2")) {
  for (case in c("1", "2", "3")) {
    
    filt <- case_when(
      case == "1" ~ c("24hr", "48hr"),
      case == "2" ~ c("48hr", "72hr"),
      case == "3" ~ c("24hr", "72hr")
    )
    
    title_case <- case_when(
      case == "1" ~ "48hr vs. 24hr",
      case == "2" ~ "72hr vs. 48hr",
      case == "3" ~ "72hr vs. 24hr"
    )
    
    x_ax <- case_when(
      case == "1" ~ "24hr",
      case == "2" ~ "48hr",
      case == "3" ~ "24hr"
    )
    
    y_ax <- case_when(
      case == "1" ~ "48hr",
      case == "2" ~ "72hr",
      case == "3" ~ "72hr"
    )
    
    cor_test_df <- treatment_diff_test_siRNA_plot %>%
      filter(siRNA_condition == siRNA) %>%
      filter(timepoint %in% filt) %>%
      filter(id %in% treatment_diff_test_siRNA$id) %>%
      select(id, normalized_effect, timepoint) %>%
      pivot_wider(names_from = timepoint,
                  values_from = normalized_effect,
                  id_cols = id) %>%
      tibble::column_to_rownames(var = "id") %>%
      dplyr::rename(X1 = 1, X2 = 2)
    
    cor_test <- cor(cor_test_df, method = "pearson")
    
    t <- ggplot(cor_test_df, aes(x = X1, y = X2)) +
      geom_point(size = 0.1) +
      geom_smooth(method = "lm", linewidth = 0.4, se = F) +
      annotate(geom = "text", x = ifelse(siRNA == "siEGFR1", 0.60, 0.73), y = -1, 
               label = paste0("r = ",as.character(round(cor_test[1,2], 3))),
               size = 1.25) +
      ggtitle(paste0(title_case)) +
      xlab(paste0(x_ax)) +
      ylab(paste0(y_ax)) +
      theme(axis.ticks.length = unit(0.6, "mm"),
            axis.ticks = element_line(linewidth = 0.1),
            text = element_text(size = 6),
            plot.title = element_text(size = 6, hjust = 0.5)) +
      monocle3:::monocle_theme_opts()
    assign(paste0("plot_",case), t)
  }
  
  p <- cowplot::plot_grid(plotlist = list(plot_1, plot_2, plot_3), nrow = 1)
  
  ggsave(paste0("BT112_",siRNA,"_timepoint_beta_coeff_cor_filt_union.png"), 
         plot = p, dpi = 600, height = 1.5, width = 3)
}




# ================================================================================
# Intersecting DEGs between siRNA for a "signature"
# ================================================================================

treatment_diff_test_siRNA_union <- treatment_diff_test %>%
  filter(grepl(pattern = "siRNA_condition", term)) %>%
  filter(q_value <= 0.05) %>%
  mutate(direction = case_when(
    normalized_effect < 0 ~ "Downregulated",
    normalized_effect > 0 ~ "Upregulated"
  ))


upset.list <- list()
for (time in c("24hr", "48hr", "72hr")) {
  for (siRNA in c("siEGFR1", "siEGFR2")) {
    upset.list[[paste0(siRNA,"_", time, "_up")]] <- treatment_diff_test_siRNA_union %>%
      filter(siRNA_condition == siRNA & timepoint == time & direction == "Upregulated") %>%
      pull(id)
    upset.list[[paste0(siRNA,"_", time,"_down")]] <- treatment_diff_test_siRNA_union %>%
      filter(siRNA_condition == siRNA & timepoint == time & direction == "Downregulated") %>%
      pull(id)
  }
}


intersection <- ComplexHeatmap::make_comb_mat(upset.list,
                                              mode = "distinct")

intersect_plot <- ComplexHeatmap::UpSet(intersection[comb_size(intersection) > 5], pt_size = unit(2, "mm"), lwd = 0.9,
                                        comb_order = order(comb_size(intersection[comb_size(intersection) > 5])), 
                                        # set_order = c("WT_up",
                                        #               "WT_down",
                                        #               "KO_up",
                                        #               "KO_down"),
                                        top_annotation = upset_top_annotation(intersection[comb_size(intersection) > 5],
                                                                              # ylim = c(0, 500),
                                                                              add_numbers = T,
                                                                              numbers_rot = 0,
                                                                              numbers_gp = gpar(fontsize = 6),
                                                                              height = unit(1, "cm"),
                                                                              axis_param = list(gp = gpar(fontsize = 6)),
                                                                              annotation_name_gp = gpar(fontsize = 6),
                                                                              show_annotation_name = T),
                                        right_annotation = upset_right_annotation(intersection[comb_size(intersection) > 5],
                                                                                  add_numbers = FALSE,
                                                                                  width = unit(1, "cm"),
                                                                                  axis_param = list(gp = gpar(fontsize = 6)),
                                                                                  annotation_name_gp = gpar(fontsize = 6),
                                                                                  show_annotation_name = T),
                                        row_names_gp = gpar(fontsize = 6, fontface = "bold"), 
                                        column_title = "siRNA DEG Intersection (FDR < 0.05)",
                                        column_title_gp = gpar(fontsize = 8, fontface = "bold"))

pdf("siRNA_DEG_direction_intersection.pdf",width=5, height=3, compress = F)
draw(intersect_plot)
dev.off()

saveRDS(upset.list, "siRNA_DEG_within_time_direction_list.RDS")

