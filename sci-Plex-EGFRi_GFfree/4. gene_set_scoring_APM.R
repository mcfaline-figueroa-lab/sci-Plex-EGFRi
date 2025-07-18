library(devtools)
library(parallel)
library(ggplot2)
library(ggridges)
library(gridExtra)
library(plyr)
library(dplyr)
library(tidyr)
library(pheatmap)
library(piano)
library(monocle3)
this.dir <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(this.dir))
rm(this.dir)

source("cell_cycle.R")
cc.genes = readRDS("cc.genes.RDS")

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

cds <- readRDS("sci-Plex-EGFRi-GFfree-cds.RDS")
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

# ===================================================================================
# Bulking drugs together by TDC4 or other
# ===================================================================================

gene_set_score_diff_test.list <- list()

for(PDCL in c("BT112", "BT333")) {
  for (condition in c("Complete", "FGF-alone", "NoGF")) {
    
    for (combined in c("TDC4", "Other")) {
      if (combined == "TDC4") {
        score_df_subset <- colData(cds) %>% as.data.frame() %>%
          filter(cell_line == PDCL &
                   GF_condition == condition &
                   drug %in% c("DMSO", "Tyrphostin9", "AG-18", 
                               "AG-490", "AG-494", "AG-555", "TyrphostinAG-528"
                   ))
        
        model_to_test <-
          paste0("APM_MHC_CI_score ~ log(dose + 0.1)")
        
        diff_test_results <-
          speedglm::speedglm(data = score_df_subset,
                             formula = model_to_test)
        diff_test_results <- broom::tidy(diff_test_results)
        
        gene_set_score_diff_test.list[[PDCL]][[condition]][[combined]] <-
          data.frame(
            cell_line = PDCL,
            signature = "APM_MHC_CI_score",
            GF_condition = condition,
            drug = "TDC4",
            term = diff_test_results$term,
            beta = diff_test_results$estimate,
            p_value = diff_test_results$p.value
          )
        } else {
            score_df_subset <- colData(cds) %>% as.data.frame() %>%
              filter(cell_line == PDCL &
                       GF_condition == condition &
                       !drug %in% c("Tyrphostin9", "AG-18", 
                                   "AG-490", "AG-494", "AG-555", 
                                   "TyrphostinAG-528"
                       ))
            
            model_to_test <-
              paste0("APM_MHC_CI_score ~ log(dose + 0.1)")
            
            diff_test_results <-
              speedglm::speedglm(data = score_df_subset,
                                 formula = model_to_test)
            diff_test_results <- broom::tidy(diff_test_results)
            
            gene_set_score_diff_test.list[[PDCL]][[condition]][[combined]] <-
              data.frame(
                cell_line = PDCL,
                signature = "APM_MHC_CI_score",
                GF_condition = condition,
                drug = "Other",
                term = diff_test_results$term,
                beta = diff_test_results$estimate,
                p_value = diff_test_results$p.value)
          }
      }
    message("Finished ", PDCL, " ", condition)
    }



  }


gene_set_score_diff_test_results.list <- gene_set_score_diff_test.list

for(PDCL in c("BT112", "BT333")){
  
  for(condition in c("Complete", "FGF-alone", "NoGF")){
    
    gene_set_score_diff_test_results.list[[PDCL]][[condition]] <- do.call("rbind",
                                                                          gene_set_score_diff_test.list[[PDCL]][[condition]])
    
  }
  
}

for(PDCL in c("BT112", "BT333")){
  
  gene_set_score_diff_test_results.list[[PDCL]] <- do.call("rbind",gene_set_score_diff_test_results.list[[PDCL]])
  gene_set_score_diff_test_results.list[[PDCL]]$q_value <- p.adjust(gene_set_score_diff_test_results.list[[PDCL]]$p_value, 
                                                                    method = "BH")
  
}

gene_set_score_diff_test_results <- do.call("rbind",gene_set_score_diff_test_results.list)
saveRDS(gene_set_score_diff_test_results, 
        file = "APM_gene_set_score_diff_test_results.RDS")

gene_set_score_diff_test_results <- readRDS("APM_gene_set_score_diff_test_results.RDS")
gene_set_score_for_plot <- gene_set_score_diff_test_results %>% filter(grepl(x = term, pattern = "dose"))

ggplot(gene_set_score_for_plot %>%
         mutate(plot_label = case_when(
           drug == "TDC4" ~ paste0(drug),
           # drug == "TDC4" ~ paste0(GF_condition,"_",drug),
           TRUE ~ NA
         )) %>%
         mutate(GF_condition_plot = case_when(
           q_value < 0.05 & drug == "TDC4" ~ GF_condition,
           TRUE ~ NA
         )) %>%
         filter(cell_line == "BT112"),
       aes(x = beta, y = -log10(q_value))) +
  # facet_wrap(~cell_line, scales = "free", nrow = 2) +
  geom_point(aes(color = GF_condition_plot),
             size = 1.5,
             show.legend = F) +
  # geom_point(aes(color = ifelse(drug == "TDC4", TRUE, FALSE)),
  #            size = 1.25,
  #            show.legend = F) +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 5, "Set1")[3:5]) +
  scale_x_continuous(limits = c(-0.2, 0.3)) +
  xlab("Beta Coefficient") +
  ylab("-Log10 FDR") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.25) +
  geom_vline(xintercept = 0.05, linetype = 2, linewidth = 0.25) +
  geom_vline(xintercept = -0.05, linetype = 2, linewidth = 0.25) +
  ggrepel::geom_text_repel(aes(label = plot_label), size = 2, max.overlaps = 5) +
  theme(text = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.ticks = element_line(linewidth = unit(0.35, "mm")),
        axis.ticks.length = unit(0.75, "mm")) +
  monocle3:::monocle_theme_opts()
ggsave("BT112_GF_TDC_APM-MHC-CI_glm_results_stacked.png",
       dpi = 600, width = 1.5, height = 1.25)

ggplot(gene_set_score_for_plot %>%
         mutate(plot_label = case_when(
           drug == "TDC4" ~ paste0(drug),
           # drug == "TDC4" ~ paste0(GF_condition,"_",drug),
           TRUE ~ NA
         )) %>%
         mutate(GF_condition_plot = case_when(
           q_value < 0.05 & drug == "TDC4" ~ GF_condition,
           TRUE ~ NA
         )) %>%
         filter(cell_line == "BT333"),
       aes(x = beta, y = -log10(q_value))) +
  # facet_wrap(~cell_line, scales = "free", nrow = 2) +
  geom_point(aes(color = GF_condition_plot),
             size = 1.5,
             show.legend = F) +
  # geom_point(aes(color = ifelse(drug == "TDC4", TRUE, FALSE)),
  #            size = 1.25,
  #            show.legend = F) +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 5, "Set1")[3:5]) +
  scale_x_continuous(limits = c(-0.1, 0.15)) +
  xlab("Beta Coefficient") +
  ylab("-Log10 FDR") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.25) +
  geom_vline(xintercept = 0.05, linetype = 2, linewidth = 0.25) +
  geom_vline(xintercept = -0.05, linetype = 2, linewidth = 0.25) +
  ggrepel::geom_text_repel(aes(label = plot_label), size = 2, max.overlaps = 5) +
  theme(text = element_text(size = 7),
        axis.title = element_text(size = 8),
        axis.ticks = element_line(linewidth = unit(0.35, "mm")),
        axis.ticks.length = unit(0.75, "mm")) +
  monocle3:::monocle_theme_opts()
ggsave("BT333_GF_TDC_APM-MHC-CI_glm_results_stacked.png",
       dpi = 600, width = 1.5, height = 1.25)

# ===================================================================================
# TDC DEG test to see beta coefficients that vary with dose
# ===================================================================================

colData(cds)$TDC <- case_when(
  colData(cds)$drug == "DMSO" ~ "DMSO",
  colData(cds)$drug %in% c("Tyrphostin9", "AG-18", 
                           "AG-490", "AG-494", "AG-555", "TyrphostinAG-528") ~ "TDC4",
  TRUE ~ "Other TDC"
)

gene_set_score_diff_test.list <- list()

for(PDCL in c("BT112", "BT333")) {
  for (condition in c("Complete", "FGF-alone", "NoGF")) {
    
    for (treatment in c("TDC4", "Other TDC")) {
      cds_subset <- cds[expressed_genes, 
                        colData(cds)$cell_line == PDCL &
                          colData(cds)$GF_condition == condition]
      
      model_to_test <-
        paste0("APM_MHC_CI_score ~ log(dose + 0.1)")
      
      diff_test_results <-
        speedglm::speedglm(data = score_df_subset,
                           formula = model_to_test)
      diff_test_results <- broom::tidy(diff_test_results)
      
      gene_set_score_diff_test.list[[PDCL]][[condition]][[treatment]] <-
        data.frame(
          cell_line = PDCL,
          signature = "APM_MHC_CI_score",
          GF_condition = condition,
          drug = treatment,
          term = diff_test_results$term,
          beta = diff_test_results$estimate,
          p_value = diff_test_results$p.value
        )
      
    }
    
    message("Finished ", PDCL, " ", condition)
  }
}

# ================================================================================
# Violin plots of DMSO or Tyrphostin9 cells alone
# ================================================================================
APM_plot <- colData(cds) %>% as.data.frame() %>%
  group_by(cell_line) %>%
  mutate(APM_up = scale(APM_MHC_CI_score)) %>%
  filter(drug %in% c("DMSO"))

library(gghalves)
split_violin_data <- APM_plot %>%
  dplyr::mutate(`GF Condition` = GF_condition) %>%
  mutate(`GF Condition` = factor(`GF Condition`, 
                                 levels = c("Complete", "FGF-alone", "NoGF"))) %>%
  mutate(x_axis = case_when(
    `GF Condition` %in% c("Complete", "FGF-alone") ~ 0,
    TRUE ~ 1
  )) %>%
  mutate(siRNA_split = case_when(`GF Condition` != "Complete" ~ "siEGFR",
                                 TRUE ~ "siNTC")) %>%
  bind_rows(APM_plot %>% 
              filter(`GF_condition` == "Complete") %>% 
              mutate(x_axis = 1) %>% 
              mutate(siRNA_split = "siNTC") %>%
              mutate(`GF Condition` = "Complete"))

split_violin_data_mean <- split_violin_data %>%
  dplyr::group_by(cell_line, siRNA_split, x_axis) %>%
  dplyr::summarise(mean_expression = mean(APM_up), median_expression = median(APM_up), .groups = "drop")

ggplot() +
  geom_half_violin(data = split_violin_data %>%
                     mutate(siRNA_split = factor(siRNA_split, levels = c("siNTC", "siEGFR"))),
                   aes(x = as.factor(x_axis), 
                       y = APM_up,
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
             aes(x = as.factor(x_axis), y = APM_up,
                 color = `GF Condition`), size = 0, stroke = 0) +
  # scale_y_log10() +
  ylab("Scaled APM MHC-CI Score") +
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


APM_plot <- colData(cds) %>% as.data.frame() %>%
  group_by(cell_line) %>%
  mutate(APM_up = scale(APM_MHC_CII_score)) %>%
  filter(drug %in% c("DMSO")) %>%
  mutate(APM_up = case_when(
    APM_up > 4 ~ 4,
    TRUE ~ APM_up
  ))

split_violin_data <- APM_plot %>%
  dplyr::mutate(`GF Condition` = GF_condition) %>%
  mutate(`GF Condition` = factor(`GF Condition`, 
                                 levels = c("Complete", "FGF-alone", "NoGF"))) %>%
  mutate(x_axis = case_when(
    `GF Condition` %in% c("Complete", "FGF-alone") ~ 0,
    TRUE ~ 1
  )) %>%
  mutate(siRNA_split = case_when(`GF Condition` != "Complete" ~ "siEGFR",
                                 TRUE ~ "siNTC")) %>%
  bind_rows(APM_plot %>% 
              filter(`GF_condition` == "Complete") %>% 
              mutate(x_axis = 1) %>% 
              mutate(siRNA_split = "siNTC") %>%
              mutate(`GF Condition` = "Complete"))

split_violin_data_mean <- split_violin_data %>%
  dplyr::group_by(cell_line, siRNA_split, x_axis) %>%
  dplyr::summarise(mean_expression = mean(APM_up), 
                   median_expression = median(APM_up), 
                   .groups = "drop")

ggplot() +
  geom_half_violin(data = split_violin_data %>%
                     mutate(siRNA_split = factor(siRNA_split, levels = c("siNTC", "siEGFR"))),
                   aes(x = as.factor(x_axis), 
                       y = APM_up,
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
             aes(x = as.factor(x_axis), y = APM_up,
                 color = `GF Condition`), size = 0, stroke = 0) +
  # scale_y_log10() +
  ylab("Scaled APM MHC-CII Score") +
  scale_y_continuous(breaks = seq(0,4,1), 
                     labels = c("0", "1", "2", "3", ">4")) +
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
ggsave("DMSO_alone_GF_condition_violin_plots_APM_MHC-CII.png",
       dpi = 600, height = 1.25, width = 2.5)
