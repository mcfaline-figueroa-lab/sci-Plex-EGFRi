library(tidyverse)
library(monocle3)
library(ComplexHeatmap)
this.dir <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(this.dir))
rm(this.dir)

cds <- readRDS("sci-Plex-EGFRi-siRNA-cds.RDS")
colData(cds)$cell_line <- "BT112"

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
# ================================================================================
# sci-Plex-GxE adaptive resistance signature
# ================================================================================
adaptive_up_genes <- read_csv("adaptive_program.csv",
                              col_names = TRUE) %>%
  filter(grepl(pattern = "pregulated", gene_module))

colData(cds)$adaptive_up <- calculate_aggregate_expression_score(cds, 
                                                                 signature_genes = adaptive_up_genes$gene_short_name,
                                                                 from_id = FALSE)

# ================================================================================
# Violin plots comparing to NTC by timepoint
# ================================================================================
adaptive_plot <- colData(cds) %>% as.data.frame() %>%
  group_by(cell_line) %>%
  mutate(adaptive_up_scaled = scale(adaptive_up) %>% as.vector())

library(gghalves)
split_violin_data <- adaptive_plot %>%
  dplyr::mutate(siRNA = siRNA_condition) %>%
  mutate(`siRNA` = factor(siRNA, 
                                 levels = c("siNTC", "siEGFR1", "siEGFR2"))) %>%
  mutate(x_axis = case_when(
    `siRNA` %in% c("siNTC", "siEGFR1") ~ 0,
    TRUE ~ 1
  )) %>%
  mutate(siRNA_split = case_when(`siRNA` != "siNTC" ~ "siEGFR",
                                 TRUE ~ "siNTC")) %>%
  bind_rows(adaptive_plot %>% 
              filter(`siRNA_condition` == "siNTC") %>% 
              mutate(x_axis = 1) %>% 
              mutate(siRNA_split = "siNTC")
  ) %>%
  mutate(siRNA_condition = factor(siRNA_condition, levels = c("siNTC", "siEGFR1", "siEGFR2"))) %>%
  filter(timepoint == "48hr")

split_violin_data_mean <- split_violin_data %>%
  group_by(cell_line, siRNA_split, x_axis, timepoint) %>%
  summarise(mean_expression = mean(adaptive_up), median_expression = median(adaptive_up), .groups = "drop")

ggplot() +
  geom_half_violin(data = split_violin_data %>%
                     mutate(siRNA_split = factor(siRNA_split, levels = c("siNTC", "siEGFR"))),
                   aes(x = as.factor(x_axis), 
                       y = adaptive_up,
                       split = siRNA_split, 
                       fill = `siRNA_condition`),
                   position = "identity",
                   show.legend = F, 
                   scale = "area",
                   linewidth = 0.2) +
  facet_wrap(~timepoint, scales = "free_x") +
  geom_point(data = split_violin_data_mean, 
             aes(x = as.factor(x_axis), y = mean_expression), 
             size = 0.4,
             position = position_nudge(x = ifelse(split_violin_data_mean$siRNA_split == "siNTC", 
                                                  -0.1, 0.1))) +
  geom_point(data = split_violin_data, 
             aes(x = as.factor(x_axis), y = adaptive_up,
                 color = `siRNA_condition`), size = 0, stroke = 0) +
  ylab("Scaled Adaptive Score") +
  scale_color_manual(name = "siRNA", 
                     values = RColorBrewer::brewer.pal(name = "Dark2",n = 3),
                     breaks = c("siNTC", "siEGFR1", "siEGFR2"),
                     labels = c("siNTC", "siEGFR1\nexon 9", "siEGFR2\nexon 2")) +
  scale_fill_manual(name = "siRNA", 
                    values = RColorBrewer::brewer.pal(name = "Dark2",n = 3)) +
  scale_y_continuous(limits = c(5.9, 6.8)) +
  theme(text = element_text(size = 7),
        axis.text.x = element_blank(),
        axis.ticks.length.y = unit(0.75,"mm"),
        axis.ticks.y = element_line(linewidth = 0.2),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line = element_line(colour = 'black', linewidth = 0.1),
        legend.position = "bottom",
        legend.key.height = unit(4, "mm"),
        legend.key.width = unit(2, "mm"),
        legend.margin = margin(t = -6, l = -5),
        legend.title = element_text(size = 6)
  ) +
  guides(color = "none") +
  monocle3:::monocle_theme_opts()
ggsave("BT112_siRNA_adaptive_score_violin_plots_mean_48hr.png",
       dpi = 600, height = 1.2, width = 1.2)

