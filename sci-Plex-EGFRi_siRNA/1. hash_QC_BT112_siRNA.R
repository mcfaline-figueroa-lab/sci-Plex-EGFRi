library(tidyverse)
library(monocle3)
this.dir <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(this.dir))
rm(this.dir)

cds <- readRDS("raw_data/cds_precell_prehash.RDS")
hash <- read_tsv("raw_data/hashTable.out", col_names = c("sample", "cell_ID", "oligo", "hash_umis"))

cds_col <- colData(cds) %>% as.data.frame() %>%
  mutate(cell_ID=Cell) %>%
  tidyr::separate(col = Cell, into = c('P7', 'P5', NA, NA, 'n1', NA, NA, 'n2'), sep = "_") %>%
  mutate(RT = paste('RT_BC_', n1, sep = '')) %>%
  mutate(Lig = paste('Lig_BC_', n2, sep = '')) %>%
  select(-n1, -n2)

RT_map <- read_csv("raw_data/RT_Map.csv", col_names = TRUE)
cds_col <- left_join(cds_col, RT_map, by = c("RT" = "RT"))

hashTable_summary <- hash %>%
  group_by(cell_ID) %>%
  mutate(proportion = hash_umis/sum(hash_umis),total_hash_umis_per_cell = sum(hash_umis)) %>%
  arrange(desc(proportion)) %>%
  mutate(rank = row_number()) %>%
  dplyr::filter(rank %in% c(1,2)) %>%
  mutate(top_to_second_best_ratio = ifelse(sum(rank) > 1, proportion[rank == 1]/proportion[rank == 2], 10)) %>%
  dplyr::filter(rank == 1)

cds_col <- left_join(cds_col, hashTable_summary, by = c("sample" = "sample", "cell_ID" = "cell_ID"))

cds_col_temp <- cds_col %>%
  separate(oligo, into = c("siRNA_condition",
                           "replicate",
                           "hash_id"), 
           sep = "_",
           remove = F
           ) %>%
  separate(hash_id,
           into = c("hash_plate", "hash_well"),
           sep = 2,
           remove = F) %>%
  dplyr::rename(n_umi = n.umi)

colData(cds) <- DataFrame(cds_col_temp)
colData(cds)$Cell <- colData(cds)$cell_ID
row.names(colData(cds)) <- colData(cds)$cell_ID

colData(cds)$log10_umi <- colData(cds)$n_umi %>% log10()
mt_genes <- rowData(cds) %>% as.data.frame() %>%
  filter(grepl("MT-",gene_short_name) == TRUE)
mt_genes_id <- rownames(mt_genes)
colData(cds)$percent_mito <- 100 * (colSums(exprs(cds)[mt_genes_id,])/colSums(exprs(cds)))

kp <- ggplot(colData(cds) %>% as.data.frame() %>% arrange(desc(n_umi)) %>% mutate(cell_rank = dplyr::row_number()), aes(x= log10(cell_rank), y = log10(n_umi))) +
  geom_point(size = 0.5, stroke = 0) +
  geom_hline(linewidth = 0.2, yintercept = log10(300)) +
  monocle3:::monocle_theme_opts() +
  scale_y_continuous(breaks = c(1,2,3,4,5)) +
  scale_x_continuous(breaks = c(0,1,2,3,4,5)) +
  xlab("Log10 CellRank") +
  ylab("Log10 UMI") +
  guides(guides(color = guide_legend(override.aes = list(size=2)))) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom", 
        legend.key.width = unit(0.5,"line"),
        legend.key.height = unit(0.2,"line"),
        text = element_text(size = 6))
ggsave(plot = kp, "QC/Kneeplot_300umi.png", dpi = 600, height = 1.5, width = 1.5)

kp <- ggplot(colData(cds) %>% as.data.frame() %>% filter(n_umi >= 10)) +
  ggridges::geom_density_ridges(aes(x = log10(n_umi), y = sample, fill = sample), 
                                alpha = 0.5, 
                                rel_min_height = 0.001,
                                position = position_nudge(y = -0.5),
                                show.legend = FALSE) +
  geom_vline(xintercept = log10(300)) +
  monocle3:::monocle_theme_opts() +
  xlab("Log10 UMI") +
  ylab("Density") +
  scale_x_continuous(breaks = seq(0,5,1)) +
  viridis::scale_fill_viridis(option = "B", discrete = TRUE)
ggsave(plot = kp,"QC/umi_density_300umi.png", dpi = 600, height = 2, width = 2.5)

ggplot(colData(cds) %>% as.data.frame() %>% filter(n_umi >= 100)) +
  ggridges::geom_density_ridges(aes(x = log10(hash_umis), y = sample, fill = sample), 
                                alpha = 0.5, 
                                rel_min_height = 0.001,
                                position = position_nudge(y = -0.5),
                                show.legend = FALSE) +
  geom_vline(xintercept = log10(5)) +
  monocle3:::monocle_theme_opts() +
  xlab("Log10 Hash UMI") +
  ylab("Density") +
  viridis::scale_fill_viridis(option = "B", discrete = TRUE)
ggsave("QC/hash_density_100umi.png", dpi = 600, height = 2, width = 2.5)

ggplot(colData(cds) %>% as.data.frame() %>% filter(n_umi >= 100)) +
  ggridges::geom_density_ridges(aes(x = log10(top_to_second_best_ratio), y = sample, fill = sample), 
                                alpha = 0.5, 
                                rel_min_height = 0.001,
                                position = position_nudge(y = -0.5),
                                show.legend = FALSE) +
  geom_vline(xintercept = log10(2.5)) +
  monocle3:::monocle_theme_opts() +
  xlab("Log10 Top to Second Best Ratio") +
  ylab("Density") +
  viridis::scale_fill_viridis(option = "B", discrete = TRUE)
ggsave("QC/hash_ratio_density_100umi.png", dpi = 600, height = 2, width = 3)

ggplot(colData(cds) %>% as.data.frame() %>% filter(n_umi >= 100) %>% arrange(desc(hash_umis)) %>% dplyr::mutate(rank = row_number()), 
       aes(x = log10(rank), y = log10(hash_umis))) +
  geom_point(size = 0.5, stroke = 0) +
  geom_hline(linewidth = 0.2, yintercept = log10(10)) +
  monocle3:::monocle_theme_opts() +
  scale_y_continuous(breaks = c(1,2,3,4,5)) +
  scale_x_continuous(breaks = c(0,1,2,3,4,5)) +
  xlab("Log10 CellRank") +
  ylab("Log10 Hash UMI") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom", 
        legend.key.width = unit(0.5,"line"),
        legend.key.height = unit(0.2,"line"),
        text = element_text(size = 6))
ggsave("QC/hash_kneeplot.png", dpi = 600, height = 1.5, width = 2)

ggplot(colData(cds) %>% as.data.frame() %>% filter(n_umi >= 300 & hash_umis >= 5),
       aes(x = sample, y = percent_mito)) +
  geom_violin() +
  monocle3:::monocle_theme_opts() +
  viridis::scale_fill_viridis(option = "B", discrete = TRUE)

# Final QC cutoffs
cds <- cds[,colData(cds)$n_umi >= 500 & 
             !is.na(colData(cds)$hash_umis) &
             colData(cds)$hash_umis >= 4 & 
             colData(cds)$top_to_second_best_ratio >= 2]

cds <- estimate_size_factors(cds)
cds <- detect_genes(cds)

source("cell_cycle.R")
cc.genes <- readRDS("cc.genes.RDS")
cds <- estimate_cell_cycle(cds, g1s_markers = cc.genes$s.genes, g2m_markers = cc.genes$g2m.genes)

saveRDS(cds, "sci-Plex-EGFRi-siRNA-cds.RDS")

# Final cell counts
n_umi_per_PCR_tot <- colData(cds) %>% as.data.frame() %>%
  summarise(count = n(), mean_umi_per_cell = round(mean(n_umi), 3), 
            med_umi_per_cell = median(n_umi), total_umi = sum(n_umi), median_hash_umi = median(hash_umis))
write_tsv(n_umi_per_PCR_tot, "QC/counts_umi_post_QC_filtering.txt")

temp_count <- colData(cds) %>% as.data.frame() %>% 
  group_by(timepoint,siRNA_condition, replicate) %>%
  dplyr::summarise(count = n(), .groups = "drop") %>%
  summarise(med_count = median(count))
write_tsv(temp_count, "QC/counts_per_siRNA_time_post_QC_filtering.txt")

temp_count <- colData(cds) %>% as.data.frame() %>% 
  group_by(siRNA_condition, drug, dose) %>%
  dplyr::summarise(count = n(), .groups = "drop") %>%
  dplyr::summarise(median_count_per_drug_dose = median(count),
                   mean_count_per_drug_dose = mean(count),
                   total_count = sum(count))
write_tsv(temp_count, "QC/counts_post_QC_filtering.txt")

# ==============================================================================================
# BT112 siRNA KD
# ==============================================================================================
test <- cds
colData(test)$condition <- paste0(colData(test)$timepoint, "_", colData(test)$siRNA_condition)
p1 <- plot_percent_cells_positive(cds_subset = test[rowData(test)$gene_short_name %in% c("EGFR", "MKI67"),], 
                                  group_cells_by = "condition",
                                  min_expr = 0) +
  theme()
temp_percent_pos <- p1$data


test <- temp_percent_pos %>% 
  separate(condition, into = c("timepoint", "siRNA_condition")) %>%
  filter(feature_label == "EGFR") %>%
  mutate(siRNA_condition = factor(siRNA_condition, levels = c("siNTC", "siEGFR1", "siEGFR2"))) %>%
  filter(timepoint == "24hr") %>%
  mutate(target_fraction_low = target_fraction_low/target_fraction_mean[siRNA_condition == "siNTC"],
         target_fraction_high = target_fraction_high/target_fraction_mean[siRNA_condition == "siNTC"]) %>%
  mutate(target_fraction_mean = target_fraction_mean/target_fraction_mean[siRNA_condition == "siNTC"])

ggplot(test,
       aes(x = siRNA_condition, y = target_fraction_mean)) +
  geom_bar(stat = "identity", aes(fill = siRNA_condition), show.legend = F) +
  geom_linerange(inherit.aes = F, aes(x = siRNA_condition,
                                      ymin = target_fraction_low,
                                      ymax = target_fraction_high)) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  ylab("% Cells Pos (Relative to siNTC)") +
  xlab("siRNA Condition") +
  theme(plot.title = element_text(size = 6, hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        text = element_text(size = 7, color = "black"),
        axis.ticks.length = unit(0.75, "mm"),
        axis.ticks = element_line(linewidth = unit(0.35, "mm"))) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(n = 3, name = "Dark2")) +
  monocle3:::monocle_theme_opts()
ggsave(filename = "EGFR_expression/BT112_siRNA_EGFR_percent_cells_pos_min0_bar.png",
       dpi = 600, height = 1.25, width = 1.25)
