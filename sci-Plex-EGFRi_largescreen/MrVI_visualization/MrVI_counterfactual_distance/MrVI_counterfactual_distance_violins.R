library(monocle3)
library(tidyverse)
library(ComplexHeatmap)

this.dir <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(this.dir))
rm(this.dir)

# ================================================================================
# single cell distance-distance violins
# ================================================================================

TDC_df <- read_csv("SharedPDCL_MrVI_Drug_DRC_BinaryMembership_HierarchicalClusters_FDRfilt_TDCs.csv")

counter_dist_df <- read_csv("BT112_SingleCellCounterfactualSampleDistance.csv") %>%
  select(X1 = 1, counter_factual_distance = 2) %>%
  separate(col = X1, into = c("drug", "dose"), sep = "_") %>%
  mutate(dose = case_when(
    !drug %in% c("DMSO", "Media", "PBS", "Panitumumab") ~ as.double(dose) - 1,
    TRUE ~ as.double(dose))) %>%
  unite(col = "drug_dose", sep = "_", drug, dose, remove = FALSE) %>%
  mutate(dose = as.factor(dose))

plot1 <- ggplot(counter_dist_df %>% mutate(counter_factual_distance = log2(counter_factual_distance + 1)) %>%
         left_join(TDC_df) %>%
         mutate(drug = fct_reorder(drug, `Drug Cluster`)),        
       aes(x = dose, y = counter_factual_distance)) +
  facet_wrap(~drug) +
  geom_violin(aes(fill = dose), linewidth = 0.1) +
  stat_summary(fun = "mean", size = 0.01) +
  xlab("Log10 (Dose)") +
  ylab("Log2 (Counterfactual Distance + 1)") + 
  viridis::scale_fill_viridis(option = "magma", discrete = TRUE) +
  guides(fill = guide_legend(title = "Log10 Dose", override.aes = list(size = 2))) +
  theme(text =element_text(size = 6),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.5, "cm"),
        axis.line = element_line(linewidth = 0.5),
        axis.ticks = element_line(linewidth = 0.1)) +
  monocle3:::monocle_theme_opts()
ggsave(plot = plot1, filename = "single_cell_distance_violin_plots_BT112.png",
       dpi = 900, width = 7.5, height = 11)


counter_dist_df <- read_csv("BT228_SingleCellCounterfactualSampleDistance.csv") %>%
  select(X1 = 1, counter_factual_distance = 2) %>%
  separate(col = X1, into = c("drug", "dose"), sep = "_") %>%
  mutate(dose = case_when(
    !drug %in% c("DMSO", "Media", "PBS", "Panitumumab") ~ as.double(dose) - 1,
    TRUE ~ as.double(dose))) %>%
  unite(col = "drug_dose", sep = "_", drug, dose, remove = FALSE) %>%
  mutate(dose = as.factor(dose))

plot1 <- ggplot(counter_dist_df %>% mutate(counter_factual_distance = log2(counter_factual_distance + 1)) %>%
                  left_join(TDC_df) %>%
                  mutate(drug = fct_reorder(drug, `Drug Cluster`)),        
                aes(x = dose, y = counter_factual_distance)) +
  facet_wrap(~drug) +
  geom_violin(aes(fill = dose), linewidth = 0.1) +
  stat_summary(fun = "mean", size = 0.01) +
  xlab("Log10 (Dose)") +
  ylab("Log2 (Counterfactual Distance + 1)") + 
  viridis::scale_fill_viridis(option = "magma", discrete = TRUE) +
  guides(fill = guide_legend(title = "Log10 Dose", override.aes = list(size = 2))) +
  theme(text =element_text(size = 6),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.5, "cm"),
        axis.line = element_line(linewidth = 0.5),
        axis.ticks = element_line(linewidth = 0.1)) +
  monocle3:::monocle_theme_opts()
ggsave(plot = plot1, filename = "single_cell_distance_violin_plots_BT228.png",
       dpi = 900, width = 7.5, height = 11)

counter_dist_df <- read_csv("BT333_SingleCellCounterfactualSampleDistance.csv") %>%
  select(X1 = 1, counter_factual_distance = 2) %>%
  separate(col = X1, into = c("drug", "dose"), sep = "_") %>%
  mutate(dose = case_when(
    !drug %in% c("DMSO", "Media", "PBS", "Panitumumab") ~ as.double(dose) - 1,
    TRUE ~ as.double(dose))) %>%
  unite(col = "drug_dose", sep = "_", drug, dose, remove = FALSE) %>%
  mutate(dose = as.factor(dose))

plot1 <- ggplot(counter_dist_df %>% mutate(counter_factual_distance = log2(counter_factual_distance + 1)) %>%
                  left_join(TDC_df) %>%
                  mutate(drug = fct_reorder(drug, `Drug Cluster`)),        
                aes(x = dose, y = counter_factual_distance)) +
  facet_wrap(~drug) +
  geom_violin(aes(fill = dose), linewidth = 0.1) +
  stat_summary(fun = "mean", size = 0.01) +
  xlab("Log10 (Dose)") +
  ylab("Log2 (Counterfactual Distance + 1)") + 
  viridis::scale_fill_viridis(option = "magma", discrete = TRUE) +
  guides(fill = guide_legend(title = "Log10 Dose", override.aes = list(size = 2))) +
  theme(text =element_text(size = 6),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.5, "cm"),
        axis.line = element_line(linewidth = 0.5),
        axis.ticks = element_line(linewidth = 0.1)) +
  monocle3:::monocle_theme_opts()
ggsave(plot = plot1, filename = "single_cell_distance_violin_plots_BT333.png",
       dpi = 900, width = 7.5, height = 11)

