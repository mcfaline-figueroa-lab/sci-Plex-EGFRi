library(tidyverse)
library(monocle3)
library(ComplexHeatmap)
source("calculate_aggreg_expression.R")
this.dir <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(this.dir))
rm(this.dir)

cds <- readRDS("sci-Plex-EGFRi-GFfree-cds.RDS")

# ================================================================================
# sci-Plex-GxE adaptive resistance signature
# ================================================================================
adaptive_up_genes <- read_csv("adaptive_program.csv",
                              col_names = TRUE) %>%
  filter(grepl(pattern = "pregulated", gene_module))

colData(cds)$adaptive_up <- calculate_aggregate_expression_score(cds, 
                                                                 signature_genes = adaptive_up_genes$gene_short_name,
                                                                 from_id = FALSE)

# fitted line plots to show drugs over dose
adaptive_plot <- colData(cds) %>% as.data.frame() %>%
  group_by(cell_line) %>%
  mutate(adaptive_up = scale(adaptive_up)) %>%
  filter(drug %in% c("Osimertinib", "MTX-211")) %>%
  filter(cell_line == "BT333") %>%
  mutate(dose1 = log10(dose + 1))

ggplot(adaptive_plot %>%
         mutate(GF_condition = factor(GF_condition, 
                                      levels = c("NoGF", "FGF-alone", "Complete"))) %>%
         arrange(desc(GF_condition)),
       aes(x = dose1, y = adaptive_up, color = GF_condition, group = GF_condition)) +
  facet_wrap(~drug, nrow = 2, scales = "free_x") +
  geom_smooth(se = F, method = "glm", formula = "y ~ log(x)") +
  scale_color_manual(name = "GF Condition", 
                     values = RColorBrewer::brewer.pal(name = "Set1",n = 5)[c(5,4,3)]) +
  ylab("Scaled Adaptive Score") +
  xlab("Log10 Dose") +
  guides(color = guide_legend(ncol = 1,
                              override.aes = list(size = 2))) +
  theme(text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 10),
        legend.key.height = unit(4, "mm"),
        legend.key.width = unit(4, "mm"),
        legend.margin = margin(t = -5, l = -18),
        legend.position = "bottom") +
  monocle3:::monocle_theme_opts()
ggsave("BT333_adaptive_score_line.png",
       dpi = 600, height = 3.75, width = 2)


# Correlating adaptive responses between media conditions for each PDCL_drug
drugs <- unique(colData(cds)$drug)[unique(colData(cds)$drug) != "DMSO"]

cor_df <- data.frame()
for (PDCL in c("BT112", "BT333")) {
  
  for (treatment in drugs) {
    
    for (media in c("FGF-alone", "NoGF")) {
      temp_df <- mean_adaptive_up %>%
        ungroup() %>%
        filter(cell_line == PDCL) %>%
        filter(drug %in%  c("DMSO", treatment)) %>%
        filter(GF_condition %in% c("Complete", media)) %>%
        select(GF_condition, dose, mean_adaptive) %>%
        pivot_wider(id_cols = dose,
                    names_from = GF_condition,
                    values_from = mean_adaptive,
                    values_fill = NA) %>%
        drop_na() %>%
        tibble::column_to_rownames(var = "dose")
      
      temp_cor <- cor(temp_df, method = "pearson")
      temp_row <- data.frame("cell_line" = PDCL,
                             "GF_condition" = media,
                             "drug" = treatment,
                             "pearson" = temp_cor[1,2])
      cor_df <- bind_rows(cor_df, temp_row)
      
      ggplot(temp_df %>% dplyr::rename(X1 = 1, X2 = 2), aes(x = X1, y = X2)) +
        geom_point() +
        geom_smooth(method = "lm", se = F) +
        xlab("Complete") +
        ylab(paste0(media)) +
        annotate(geom = "text",
                 x = min(temp_df[,1]), y = max(temp_df[,2]),
                 label = paste0("r = ",round(temp_cor[1,2],3))) +
        monocle3:::monocle_theme_opts()
    }
    
  }
  
}


cor_df_BT112 <- cor_df %>%
  filter(cell_line == "BT112") %>%
  select(-cell_line) %>%
  pivot_wider(id_cols = GF_condition, 
              names_from = drug,
              values_from = pearson) %>%
  tibble::column_to_rownames(var = "GF_condition") %>%
  as.matrix()

Heatmap(matrix = cor_df_BT112,
        col = circlize::colorRamp2(breaks = c(-1,0,1),
                                   colors = c("blue", "white", "red")))

cor_df_hm <- cor_df %>%
  unite(col = "PDCL_GF", sep = "_", cell_line, GF_condition, remove = T) %>%
  pivot_wider(id_cols = PDCL_GF, 
              names_from = drug,
              values_from = pearson) %>%
  tibble::column_to_rownames(var = "PDCL_GF") %>%
  as.matrix()

library(ComplexHeatmap)
hm <- Heatmap(matrix = cor_df_hm,
        name = "Pearson\nCorrelation",
        col = circlize::colorRamp2(breaks = c(-1,0,1),
                                   colors = c("blue", "white", "red")),
        border = T, border_gp = gpar(lwd = 0.4),
        rect_gp = gpar(lwd = 0.4),
        cluster_rows = F,
        cluster_columns = T,
        clustering_method_columns = "complete",
        row_split = c(rep("BT112", 2), rep("BT333",2)),
        row_title_rot = 90, row_title_gp = gpar(fontsize = 8),
        show_row_names = F,
        column_names_rot = 45, column_names_gp = gpar(fontsize = 6),
        heatmap_legend_param = list(title_gp = gpar(fontsize = 7),
                                    labels_gp = gpar(fontsize = 6),
                                    border = T,
                                    legend_height = unit(0., "cm"),
                                    grid_width = unit(0.2, "cm")
                                    ),
        column_dend_gp = gpar(lwd = 0.75),
        column_dend_height = unit(5, "mm")
        )

hm

png("combined_GFcondition_resistance_correlation_hm_cell_line_stacked.png",width=2.5,height=1.5,
    units="in",res=1200)
ComplexHeatmap::draw(hm, 
     heatmap_legend_side = "right", annotation_legend_side = "right", legend_grouping = "original")
dev.off()

# ================================================================================
# DMSO cells alone - violin plots
# ================================================================================
adaptive_plot <- colData(cds) %>% as.data.frame() %>%
  group_by(cell_line) %>%
  mutate(adaptive_up = scale(adaptive_up)) %>%
  filter(drug %in% c("DMSO"))

library(gghalves)
split_violin_data <- adaptive_plot %>%
  dplyr::mutate(`GF Condition` = GF_condition) %>%
  mutate(`GF Condition` = factor(`GF Condition`, 
                                 levels = c("Complete", "FGF-alone", "NoGF"))) %>%
  mutate(x_axis = case_when(
    `GF Condition` %in% c("Complete", "FGF-alone") ~ 0,
    TRUE ~ 1
  )) %>%
  mutate(siRNA_split = case_when(`GF Condition` != "Complete" ~ "siEGFR",
                                 TRUE ~ "siNTC")) %>%
  bind_rows(adaptive_plot %>% 
              filter(`GF_condition` == "Complete") %>% 
              mutate(x_axis = 1) %>% 
              mutate(siRNA_split = "siNTC") %>%
              mutate(`GF Condition` = "Complete"))

split_violin_data_mean <- split_violin_data %>%
  group_by(cell_line, siRNA_split, x_axis) %>%
  summarise(mean_expression = mean(adaptive_up), median_expression = median(adaptive_up), .groups = "drop")

ggplot() +
  geom_half_violin(data = split_violin_data %>%
                     mutate(siRNA_split = factor(siRNA_split, levels = c("siNTC", "siEGFR"))),
                   aes(x = as.factor(x_axis), 
                       y = adaptive_up,
                       split = siRNA_split, 
                       fill = `GF Condition`),
                   position = "identity",
                   show.legend = F, 
                   # nudge = -.19,
                   scale = "width",
                   linewidth = 0.2) +
  facet_wrap(~cell_line, scales = "free_y") +
  geom_point(data = split_violin_data_mean, 
             aes(x = as.factor(x_axis), y = median_expression), 
             size = 0.4,
             position = position_nudge(x = ifelse(split_violin_data_mean$siRNA_split == "siNTC", 
                                                  -0.1, 0.1))) +
  geom_point(data = split_violin_data, 
             aes(x = as.factor(x_axis), y = adaptive_up,
                 color = `GF Condition`), size = 0, stroke = 0) +
  # scale_y_log10() +
  ylab("Scaled Adaptive Score") +
  scale_color_manual(name = "GF Condition", 
                     values = RColorBrewer::brewer.pal(name = "Set1",n = 5)[3:5]) +
  scale_fill_manual(name = "GF Condition", 
                     values = RColorBrewer::brewer.pal(name = "Set1",n = 5)[3:5]) +
  theme(text = element_text(size = 6),
        axis.text.x = element_blank(),
        axis.ticks.length.y = unit(0.75,"mm"),
        axis.ticks.y = element_line(linewidth = 0.2),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        legend.key.height = unit(4, "mm"),
        legend.key.width = unit(4, "mm"),
        legend.margin = margin(t = -6, l = -5),
        legend.title = element_text(size = 7)
  ) +
  # guides(color = guide_legend(title = "GF Condition", ncol = 1,
  #                             override.aes = list(shape = 16, size = 3))) +
  guides(color = "none") +
  monocle3:::monocle_theme_opts()
ggsave("DMSO_alone_GF_condition_violin_plots.png",
       dpi = 600, height = 1.25, width = 2.5)
