library(tidyverse)
library(monocle3)
library(ComplexHeatmap)
this.dir <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(this.dir))
rm(this.dir)

# DRCs signify response modules (RMs)
cds <- readRDS("cds_large_screen.RDS")

DRC_labels <- read_csv("AllCells_DRC_labels_01122024.csv", col_select = 2) %>%
  dplyr::rename(drug_response_class = 1)

colData(cds)$drug_response_class <- DRC_labels

# ===========================================================================================
# Plotting number of DEGs on the proliferation heatmap
# ===========================================================================================
gene_list <- read_csv("hvg_5percent_genes.csv", col_select = 2) %>%
  pull()

drc_labels_loop <- arrange(DRC_labels, drug_response_class) %>% distinct() %>% pull() 

diff_test.list <- list()
for (line in c("BT112", "BT228", "BT333")) {
  for (drc in drc_labels_loop) {
    if (grepl(pattern = line, drc) == TRUE) {
      colData(cds)$DEG_test <- case_when(
        colData(cds)$drug == "DMSO" ~ "DMSO",
        colData(cds)$drug_response_class == drc ~ "test_group",
        TRUE ~ NA
      )
      
      diff_test.list[[line]][[drc]] <- fit_models(cds[rowData(cds)$gene_short_name %in% gene_list,
                                                      colData(cds)$cell_line == line], 
                                                  # expression_family = "mixed-negbinomial",
                                                  model_formula_str = "~ DEG_test + replicate")
      
      diff_test.list[[line]][[drc]] <- coefficient_table(diff_test.list[[line]][[drc]]) %>% 
        dplyr::select(-model,-model_summary)
      
      diff_test.list[[line]][[drc]]$drug_response_class <- rep(drc,
                                                               nrow(diff_test.list[[line]][[drc]]))
      diff_test.list[[line]][[drc]]$cell_line <- rep(line, 
                                                     nrow(diff_test.list[[line]][[drc]]))
      
      message(paste("Processed", drc))
    }
  }
}

# when DEG is run for each cell line and each drc
diff_test_results_temp.list <- list()
for (line in c("BT112", "BT228", "BT333")) {
  diff_test_results_temp.list[[line]] <- do.call("rbind", diff_test.list[[line]])
}

diff_test_results <- do.call("rbind", diff_test_results_temp.list)
diff_test_results$q_value <- p.adjust(diff_test_results$p_value, method = "BH")

diff_test_results_filtered <- diff_test_results %>%
  filter(term == "DEG_testtest_group") %>%
  filter(q_value <= 0.05) %>%
  filter(abs(normalized_effect) >= 0.1)

rm(diff_test.list)

ggplot(diff_test_results %>% 
         filter(term == "DEG_testtest_group") %>%
         mutate(signif = case_when(
           # q_value <= 0.05 & abs(normalized_effect) >= 0.05 ~ "Yes",
           q_value <= 0.05 ~ "Yes",
           TRUE ~ "No"
         )), 
       aes(x = normalized_effect, y = -log10(q_value))) +
  facet_wrap(~drug_response_class, nrow = 8) +
  geom_point(aes(color = signif), alpha = 0.5, stroke = 0.01, size = 0.6) +
  scale_color_manual(name = "Significant", values = c("black", "red")) +
  scale_x_continuous(limits = c(-10,10)) +
  scale_y_continuous(limits = c(0,50)) +
  monocle3:::monocle_theme_opts()


for (line in c("BT112", "BT228", "BT333")) {
  diff_test_result_temp <- diff_test_results %>%
    filter(term == "DEG_testtest_group" & cell_line == line) %>%
    select(gene_short_name, drug_response_class, q_value) %>%
    distinct(gene_short_name, drug_response_class, .keep_all = T) %>%
    pivot_wider(id_cols = gene_short_name,
                names_from = drug_response_class,
                values_from = q_value)
  write_csv(x = diff_test_result_temp, 
            file = paste0(line, "_DRC_GLM_DEG_qvalue_matrix_20240125.csv"), )
}

q_value_matrix <- diff_test_results %>% 
  filter(term == "DEG_testtest_group") %>%
  mutate(signif = case_when(
    # q_value <= 0.05 & abs(normalized_effect) >= 0.05 ~ "Yes",
    q_value <= 0.05 ~ "Yes",
    TRUE ~ "No"
  )) %>%
  filter(signif == "Yes") %>%
  group_by(drug_response_class) %>%
  summarise(count = n())

