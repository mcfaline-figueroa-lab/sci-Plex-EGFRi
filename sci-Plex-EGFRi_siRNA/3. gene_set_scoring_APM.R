library(devtools)
library(parallel)
library(ggplot2)
library(ggridges)
library(gridExtra)
library(dplyr)
library(tidyr)
library(pheatmap)
library(piano)
library(monocle3)
this.dir <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(this.dir))
rm(this.dir)

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

cds <- readRDS("sci-Plex-EGFRi-siRNA-cds.RDS")
colData(cds)$cell <- colData(cds)$cell_ID

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

# ================================================================================
# Violin plots of siRNA
# ================================================================================
MHC_plot <- colData(cds) %>% as.data.frame() %>%
  mutate(adaptive_up = scale(APM_MHC_CI_score) %>% as.vector())

library(gghalves)
split_violin_data <- MHC_plot %>%
  dplyr::mutate(siRNA = siRNA_condition) %>%
  mutate(`siRNA` = factor(siRNA, 
                          levels = c("siNTC", "siEGFR1", "siEGFR2"))) %>%
  mutate(x_axis = case_when(
    `siRNA` %in% c("siNTC", "siEGFR1") ~ 0,
    TRUE ~ 1
  )) %>%
  mutate(siRNA_split = case_when(`siRNA` != "siNTC" ~ "siEGFR",
                                 TRUE ~ "siNTC")) %>%
  bind_rows(MHC_plot %>% 
              filter(`siRNA_condition` == "siNTC") %>% 
              mutate(x_axis = 1) %>% 
              mutate(siRNA_split = "siNTC")
  ) %>%
  mutate(siRNA_condition = factor(siRNA_condition, levels = c("siNTC", "siEGFR1", "siEGFR2")))

split_violin_data_mean <- split_violin_data %>%
  group_by(siRNA_split, x_axis, timepoint) %>%
  summarise(mean_expression = mean(adaptive_up), median_expression = median(adaptive_up), .groups = "drop")

ggplot() +
  geom_half_violin(data = split_violin_data %>%
                     mutate(siRNA_split = factor(siRNA_split, levels = c("siNTC", "siEGFR"))) %>%
                     mutate(adaptive_up = case_when(
                       adaptive_up > 3 ~ 3,
                       TRUE ~ adaptive_up
                     )),
                   aes(x = as.factor(x_axis), 
                       y = adaptive_up,
                       split = siRNA_split, 
                       fill = `siRNA_condition`),
                   position = "identity",
                   show.legend = F, 
                   # nudge = -.19,
                   scale = "area",
                   linewidth = 0.2) +
  facet_wrap(~timepoint, scales = "free_x") +
  geom_point(data = split_violin_data_mean, 
             aes(x = as.factor(x_axis), y = mean_expression), 
             size = 0.25,
             position = position_nudge(x = ifelse(split_violin_data_mean$siRNA_split == "siNTC", 
                                                  -0.1, 0.1))) +
  geom_point(data = split_violin_data %>%
               mutate(adaptive_up = case_when(
                 adaptive_up > 3 ~ 3,
                 TRUE ~ adaptive_up
               )), 
             aes(x = as.factor(x_axis), y = adaptive_up,
                 color = `siRNA_condition`), size = 0, stroke = 0) +
  # geom_hline(yintercept = 0, linetype = 2, linewidth = 0.1) +
  # scale_y_log10() +
  ylab("Scaled Adaptive Score") +
  scale_color_manual(name = "siRNA", 
                     values = RColorBrewer::brewer.pal(name = "Dark2",n = 3),
                     breaks = c("siNTC", "siEGFR1", "siEGFR2"),
                     labels = c("siNTC", "siEGFR1\nexon 9", "siEGFR2\nexon 2")) +
  scale_fill_manual(name = "siRNA", 
                    values = RColorBrewer::brewer.pal(name = "Dark2",n = 3)) +
  theme(text = element_text(size = 7),
        axis.text.x = element_blank(),
        axis.ticks.length.y = unit(0.75,"mm"),
        axis.ticks.y = element_line(linewidth = 0.2),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        legend.key.height = unit(4, "mm"),
        legend.key.width = unit(2, "mm"),
        legend.margin = margin(t = -6, l = -5),
        legend.title = element_text(size = 6)
  ) +
  guides(color = guide_legend(title = "siRNA", nrow = 1,
                              override.aes = list(shape = 16, size = 2))) +
  monocle3:::monocle_theme_opts()
ggsave("BT112_siRNA_APM_MHC-CI_violin_plots_mean.png",
       dpi = 600, height = 1.5, width = 2.25)

