library(tidyverse)
library(monocle3)
library(ComplexHeatmap)

this.dir <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(this.dir))
rm(this.dir)

source("cell_cycle.R")
cc.genes <- readRDS("cc.genes.RDS")
source("calculate_aggreg_expression.R")

cds_7d <- readRDS("sci-Plex-EGFRi-7d-largescreen-cds.RDS")
colData(cds_7d)$time_point <- "7d"
cds_24hr <- readRDS("sci-Plex-EGFRi-largescreen-cds.RDS")
colData(cds_24hr)$time_point <- "24hr"

cds <- combine_cds(list(cds_24hr, cds_7d), cell_names_unique = TRUE)

cds <- estimate_size_factors(cds)
cds <- detect_genes(cds)
cds <- estimate_cell_cycle(cds, g1s_markers = cc.genes$s.genes, g2m_markers = cc.genes$g2m.genes)

adaptive_up_genes <- read_csv("adaptive_program.csv",
                              col_names = TRUE) %>%
  filter(grepl(pattern = "pregulated", gene_module))

colData(cds)$adaptive_up <- calculate_aggregate_expression_score(cds, 
                                                                 signature_genes = adaptive_up_genes$gene_short_name,
                                                                 from_id = FALSE)

mean_adaptive_up <- colData(cds) %>% 
  as.data.frame() %>% 
  group_by(cell_line, drug, dose, time_point) %>%
  dplyr::summarise(mean_adaptive = mean(adaptive_up), 
            count = n(),
            .groups = "drop") %>%
  select(-count) %>%
  group_by(cell_line) %>% # scale by PDCL
  mutate(mean_adaptive = scale(mean_adaptive)) %>%
  mutate(dose = case_when(
    drug == "Panitumumab" ~ as.character(log10(dose * 10)),
    drug %in% c("DMSO", "PBS", "Media") ~ drug,
    TRUE ~ as.character(log10(dose))
  )) %>%
  unite(col = "cell_line_drug_time_point",cell_line, drug, time_point, sep = "_") %>%
  pivot_wider(names_from = "cell_line_drug_time_point", values_from = mean_adaptive, id_cols = dose) %>%
  column_to_rownames("dose") %>%
  t() %>%
  as.data.frame() %>%
  select(-PBS, -Media) %>%
  rownames_to_column(var = "cell_line_drug_time_point") %>%
  separate(cell_line_drug_time_point, into = c("cell_line", "drug", "time_point"), remove = FALSE, sep = "_")

mean_adaptive_up_controls <- mean_adaptive_up %>%
  filter(drug %in% c("DMSO"))

mean_adaptive_up <- mean_adaptive_up %>%
  mutate(DMSO = case_when(
    cell_line == "BT112" & time_point == "24hr" ~ mean_adaptive_up_controls %>% filter(cell_line_drug_time_point == "BT112_DMSO_24hr") %>% select(DMSO) %>% pull(),
    cell_line == "BT228" & time_point == "24hr" ~ mean_adaptive_up_controls %>% filter(cell_line_drug_time_point == "BT228_DMSO_24hr") %>% select(DMSO) %>% pull(),
    cell_line == "BT333" & time_point == "24hr" ~ mean_adaptive_up_controls %>% filter(cell_line_drug_time_point == "BT333_DMSO_24hr") %>% select(DMSO) %>% pull(),
    cell_line == "BT112" & time_point == "7d" ~ mean_adaptive_up_controls %>% filter(cell_line_drug_time_point == "BT112_DMSO_7d") %>% select(DMSO) %>% pull(),
    cell_line == "BT228" & time_point == "7d" ~ mean_adaptive_up_controls %>% filter(cell_line_drug_time_point == "BT228_DMSO_7d") %>% select(DMSO) %>% pull(),
    cell_line == "BT333" & time_point == "7d" ~ mean_adaptive_up_controls %>% filter(cell_line_drug_time_point == "BT333_DMSO_7d") %>% select(DMSO) %>% pull()
  )) %>%
  filter(!drug %in% c("DMSO", "PBS", "Media")) %>%
  select(-cell_line, -drug) %>%
  pivot_longer(cols = c("1","2","3","4"), 
               names_to = "dose", 
               values_to = "agg_score") %>%
  mutate(norm_agg_score = case_when(
    is.na(agg_score) ~ NA,
    TRUE ~ agg_score
  )) %>%
  select(-agg_score) %>%
  pivot_wider(id_cols = c(cell_line_drug_time_point, DMSO), 
              names_from = dose, 
              values_from = norm_agg_score) %>%
  column_to_rownames(var = "cell_line_drug_time_point") %>%
  t()

# ================================================================================
# Comparing 24hr and 7d adaptive resistance
# ================================================================================
# arranging by cell line

mean_adaptive_up_byline <- mean_adaptive_up %>% as.data.frame() %>%
  rownames_to_column(var = "dose") %>%
  pivot_longer(2:last_col(), names_to = "cell_line_drug_time_point", values_to = "mean_adaptive") %>%
  separate(col = cell_line_drug_time_point, sep = "_", into = c("cell_line", "drug", "time_point"))

mean_adaptive_DMSO <- mean_adaptive_up_byline %>% filter(dose == "DMSO") %>%
  distinct(cell_line, dose, time_point, mean_adaptive) %>%
  mutate(drug = "DMSO")

mean_adaptive_up_10uM <- mean_adaptive_up_byline %>% filter(dose == 4) %>%
  bind_rows(mean_adaptive_DMSO)

ggplot(data = mean_adaptive_up_10uM %>%
         mutate(color_plot = case_when(
           drug == "Osimertinib" ~ "drug1",
           # drug == "MTX-211" ~ "drug3",
           TRUE ~ NA
         )),
       aes(x = time_point, y = mean_adaptive)) +
  facet_wrap(~cell_line) +
  geom_point(size = 0.1) +
  geom_line(aes(color = color_plot, size = color_plot, alpha = color_plot, group = drug),
             show.legend = F, lineend = "round") +
  scale_size_manual(na.value = 0.2, values = c(0.75,0.5,0.5)) +
  scale_alpha_manual(na.value = 0.5, values = c(1,1,1)) +
  geom_hline(aes(yintercept = 0), linetype = 2, linewidth = 0.2) +
  scale_color_manual(breaks = c("drug1", "drug3"),values = c("#F76F8E","#54F2F2")) +
  ylab("Scaled Mean Adaptive Score") +
  xlab("Drug Exposure Time") +
  theme(text = element_text(size = 6),
        axis.ticks = element_line(linewidth = 0.3),
        axis.ticks.length = unit(0.7, "mm")) +
  monocle3:::monocle_theme_opts()
ggsave("scaled_separately/adaptive_resistance_osimertinib_scaled_PDCL_timepoint.png",
       dpi = 600, height = 1.15, width = 1.4)

adaptive_delta_10uM <- mean_adaptive_up_10uM %>%
  pivot_wider(id_cols = c("cell_line", "dose", "drug"), names_from = "time_point", values_from = "mean_adaptive") %>%
  dplyr::rename(mean_adaptive_7d = `7d`, mean_adaptive_24hr = `24hr`) %>%
  mutate(delta_adaptive = mean_adaptive_7d - mean_adaptive_24hr,
         percent_change_adaptive = (mean_adaptive_7d - mean_adaptive_24hr)/abs(mean_adaptive_24hr))

proportion_delta_10uM <- readRDS("sci-Plex-EGFRi-7d-largescreen-delta-proportion-10uM.RDS") %>%
  select(proportion_24hr = `24hr`, 
         proportion_7d = `7d`,
         delta_proportion = delta, 
         percent_change_proportion = percent_change) %>%
  ungroup()

test <- full_join(adaptive_delta_10uM %>% mutate, 
                  proportion_delta_10uM, 
                  by = c("cell_line" = "cell_line",
                         "drug" = "drug"))

# Annotation of high persisting or subsiding 
test_high <- test %>%
  group_by(cell_line) %>%
  mutate(status = case_when(
    percent_change_adaptive > 0.5 ~ "Increase",
    TRUE ~ "Maintain/\nSubside"
  )) %>%
  filter(!is.na(status)) %>%
  select(cell_line, drug, status)

status_plot <- test %>% 
  filter(drug %in% c(test_high$drug)) %>%
  left_join(test_high, by = c("cell_line", "drug")) %>%
  filter(!is.na(status))

# plotting percent change in proportion
ggplot(data = status_plot,
       aes(x = status, y = percent_change_proportion,fill = status)) +
  geom_violin(aes(), show.legend = F,
              linewidth = 0.25) +
  stat_summary(fun.y=median, geom="point", color="black", stroke = 0, size = 1, show.legend = F) +
  ggpubr::stat_compare_means(method = "wilcox", label = "p.format", size = 1.75, label.x.npc = 0.4, label.y.npc = 0.9) +
  ylab("% Change Proportion") +
  scale_fill_manual(values = c("firebrick", "royalblue"), 
                    breaks = c("Increase", "Maintain/\nSubside")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 6),
        axis.ticks.length = unit(0.5, "mm"),
        axis.ticks = element_line(linewidth = 0.25),
        axis.title.x = element_blank()) +
  monocle3:::monocle_theme_opts()
ggsave("adaptive_resistance_persist_vs_subside_delta_violin_allPDCL_noP.png",
       dpi = 600, height = 1.25, width = 1.2)

# visualizing inhibitors that fall under increase or maintain/subside
ggplot(data = test %>% 
         # filter(drug %in% c(test_high$drug)) %>%
         left_join(test_high, by = c("cell_line", "drug")) %>%
         select(cell_line, drug, mean_adaptive_24hr, mean_adaptive_7d, status) %>%
         pivot_longer(cols = "mean_adaptive_24hr":"mean_adaptive_7d", 
                      names_to = "time", values_to = "mean_adaptive") %>%
         mutate(time = case_when(
           grepl(pattern = "24hr", time) ~ "24hr",
           TRUE ~ "7d"
         )) %>%
         group_by(cell_line, drug) %>%
         mutate(color_plot = case_when(
           status == "Increase" ~ "drug3",
           status == "Maintain/\nSubside" ~ "drug4",
           TRUE ~ NA
         )),
       aes(x = time, y = mean_adaptive, group = drug)) +
  facet_wrap(~cell_line) +
  geom_point(size = 0.05) +
  geom_line( aes(color = color_plot, size = color_plot, alpha = color_plot),
             show.legend = F, lineend = "round") +
  scale_size_manual(na.value = 0.2, values = c(0.25,0.4)) +
  scale_alpha_manual(na.value = 0.5, values = c(0.7,0.9,1)) +
  # geom_line(linewidth = 0.2, aes(color = col)) +
  geom_hline(aes(yintercept = 0), linetype = 2, linewidth = 0.2) +
  scale_color_manual(values = c("firebrick", "royalblue")) +
  ylab("Scaled Mean Adaptive Score") +
  xlab("Drug Exposure Time") +
  theme(text = element_text(size = 6),
        axis.ticks = element_line(linewidth = 0.3),
        axis.ticks.length = unit(0.7, "mm")) +
  monocle3:::monocle_theme_opts()
ggsave("resistance_proportion_comparison/test_adaptive_resistance_increase_newlines.png",
       dpi = 600, height = 1.15, width = 1.4)

  
# ========================================================================================
# Highlighting specific inhibitors on the adaptive comparison plot
# ========================================================================================

ggplot(data = mean_adaptive_up_10uM %>%
         # group_by(cell_line, drug) %>%
         mutate(color_plot = case_when(
           # drug == "MTX-211" ~ "drug1",
           drug == "Osimertinib" ~ "drug2",
           # drug == "DMSO" ~ "drug3",
           TRUE ~ NA
         )),
       aes(x = time_point, y = mean_adaptive)) +
  facet_wrap(~cell_line) +
  geom_point(size = 0.1) +
  geom_line(aes(color = color_plot, size = color_plot, alpha = color_plot, group = drug),
            show.legend = F, lineend = "round") +
  scale_size_manual(na.value = 0.2, values = c(0.75,0.5,0.5)) +
  scale_alpha_manual(na.value = 0.5, values = c(1,1,1)) +
  geom_hline(aes(yintercept = 0), linetype = 2, linewidth = 0.2) +
  scale_color_manual(values = c("magenta1")) +
  ylab("Scaled Mean Adaptive Score") +
  xlab("Drug Exposure Time") +
  theme(text = element_text(size = 6)) +
  monocle3:::monocle_theme_opts()
ggsave("adaptive_resistance_osimertinib_scaled_PDCL_timepoint.png",
       dpi = 600, height = 1.5, width = 2.25)


