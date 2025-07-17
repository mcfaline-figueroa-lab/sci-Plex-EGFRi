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

RT_map <- read_csv("raw_data/CellLine_RT_Map.csv", col_names = TRUE)

cds_col <- left_join(cds_col, RT_map, by = c("RT" = "RT_barcode"))

hashTable_summary <- hash %>%
  group_by(cell_ID) %>%
  mutate(proportion = hash_umis/sum(hash_umis),total_hash_umis_per_cell = sum(hash_umis)) %>%
  arrange(desc(proportion)) %>%
  mutate(rank = row_number()) %>%
  dplyr::filter(rank %in% c(1,2)) %>%
  mutate(top_to_second_best_ratio = ifelse(sum(rank) > 1, proportion[rank == 1]/proportion[rank == 2], 1)) %>%
  dplyr::filter(rank == 1)

cds_col <- left_join(cds_col, hashTable_summary, by = c("sample" = "sample", "cell_ID" = "cell_ID"))

cds_col_temp <- cds_col %>%
  separate(oligo, into = c("drug1", "dose1", "drug2", "dose2","hash"), sep = "_") %>%
  separate(hash, sep = 2, into = c("hash_plate", "hash_well")) %>%
  dplyr::rename(n_umi = n.umi)

test <- cds_col_temp %>% group_by(cell_line, drug1, dose1, drug2, dose2) %>% summarise(count = n())

colData(cds) <- DataFrame(cds_col_temp)
colData(cds)$Cell <- colData(cds)$cell_ID
row.names(colData(cds)) <- colData(cds)$cell_ID

colData(cds)$log10_umi <- colData(cds)$n_umi %>% log10()
mt_genes <- rowData(cds) %>% as.data.frame() %>%
  filter(grepl("MT-",gene_short_name) == TRUE)
mt_genes_id <- rownames(mt_genes)
colData(cds)$percent_mito <- 100 * (colSums(exprs(cds)[mt_genes_id,])/colSums(exprs(cds)))

kp <- ggplot(colData(cds) %>% as.data.frame() %>% arrange(desc(n_umi)) %>% mutate(cell_rank = dplyr::row_number()), aes(x= log10(cell_rank), y = log10(n_umi))) +
  geom_point(size = 0.5, stroke = 0, aes(color = cell_line)) +
  geom_hline(linewidth = 0.2, yintercept = log10(300)) +
  monocle3:::monocle_theme_opts() +
  scale_color_manual(name = "Cell Line",
                     breaks = c("BT112", "BT228", "BT333"),
                     values = c("darkblue", "darkred", "darkgreen")) +
  facet_wrap(~cell_line, nrow = 1) +
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
ggsave(plot = kp, "QC/Kneeplot_300umi.png", dpi = 600, height = 1.5, width = 3)

kp <- ggplot(colData(cds) %>% as.data.frame() %>% arrange(desc(n_umi)) %>% mutate(cell_rank = dplyr::row_number()), aes(x= log10(cell_rank), y = log10(n_umi))) +
  geom_point(size = 0.5, stroke = 0) +
  geom_hline(linewidth = 0.2, yintercept = log10(300)) +
  monocle3:::monocle_theme_opts() +
  scale_y_continuous(breaks = c(1,2,3,4,5)) +
  scale_x_continuous(breaks = c(0,1,2,3,4,5)) +
  xlab("Log10 CellRank") +
  ylab("Log10 UMI") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom", 
        legend.key.width = unit(0.5,"line"),
        legend.key.height = unit(0.2,"line"),
        text = element_text(size = 6))
ggsave(plot = kp, "QC/Kneeplot_300umi_combined.png", dpi = 600, height = 1, width = 1.5)

kp <- ggplot(colData(cds) %>% as.data.frame() %>% filter(n_umi >= 10)) +
  ggridges::geom_density_ridges(aes(x = log10(n_umi), y = cell_line, fill = cell_line), 
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

dim(cds[,colData(cds)$n_umi >= 300 & !is.na(colData(cds)$hash_umis) & colData(cds)$hash_umis >= 5 & colData(cds)$top_to_second_best_ratio >= 2.5])
dim(cds[,colData(cds)$n_umi >= 500])

ggplot(colData(cds) %>% as.data.frame() %>% filter(n_umi >= 100)) +
  ggridges::geom_density_ridges(aes(x = log10(hash_umis), y = cell_line, fill = cell_line), 
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
  ggridges::geom_density_ridges(aes(x = log10(top_to_second_best_ratio), y = cell_line, fill = cell_line), 
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
       aes(x = cell_line, y = percent_mito)) +
  geom_violin(aes(fill = cell_line)) +
  monocle3:::monocle_theme_opts() +
  viridis::scale_fill_viridis(option = "B", discrete = TRUE)

# ==============================================================================================

# Final QC cutoffs
cds <- cds[,colData(cds)$n_umi >= 100 & 
             !is.na(colData(cds)$hash_umis) &
             colData(cds)$hash_umis >= 5 & 
             colData(cds)$top_to_second_best_ratio >= 2.5]

cds <- estimate_size_factors(cds)
cds <- detect_genes(cds)

source("/Users/rossgiglio/Documents/McFaline-Figueroa_Lab/Experiment/20220919_combined_3run/cell_cycle.R")
cc.genes <- readRDS("/Users/rossgiglio/Documents/McFaline-Figueroa_Lab/Experiment/20220919_combined_3run/cc.genes.RDS")
cds <- estimate_cell_cycle(cds, g1s_markers = cc.genes$s.genes, g2m_markers = cc.genes$g2m.genes)

saveRDS(cds, file = "sci-Plex-EGFRi-combination-PI3Ki-cds.RDS")

cds <- readRDS("sci-Plex-EGFRi-combination-PI3Ki-cds.RDS")

# Final cell counts
n_umi_per_PCR_tot <- colData(cds) %>% as.data.frame() %>%
  summarise(count = n(), mean_umi_per_cell = mean(n_umi), med_umi_per_cell = median(n_umi), total_umi = sum(n_umi)) %>%
  mutate(P7_row = "Total")

n_umi_per_PCR <- colData(cds) %>% as.data.frame() %>%
  separate(P7, into = c(NA, "P7_row"), sep = 2) %>%
  group_by(P7_row) %>%
  summarise(count = n(), mean_umi_per_cell = mean(n_umi), med_umi_per_cell = median(n_umi), total_umi = sum(n_umi)) %>%
  bind_rows(n_umi_per_PCR_tot)

mean_counts <- colData(cds) %>% as.data.frame() %>%
  group_by(cell_line, drug1, dose1, drug2, dose2) %>%
  summarise(count = n(), .groups = "drop") %>%
  # group_by(cell_line) %>%
  summarise(average = mean(count), median = median(count)) 
