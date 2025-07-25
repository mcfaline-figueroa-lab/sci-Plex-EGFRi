library(tidyverse)
library(monocle3)
this.dir <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(this.dir))
rm(this.dir)

# ======================================================================================================

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
  separate(oligo, into = c("drug", "dose", "hash"), sep = "_") %>%
  separate(dose, sep = -2, into = c("dose", NA), convert = TRUE) %>%
  separate(hash, sep = 2, into = c("hash_plate", "hash_well")) %>%
  mutate(log10_dose = log10(dose)) %>%
  dplyr::rename(n_umi = n.umi)

test <- cds_col_temp %>% group_by(cell_line, drug, dose) %>% summarise(count = n())

colData(cds) <- DataFrame(cds_col_temp)
colData(cds)$Cell <- colData(cds)$cell_ID
row.names(colData(cds)) <- colData(cds)$cell_ID

colData(cds)$log10_umi <- colData(cds)$n_umi %>% log10()
mt_genes <- rowData(cds) %>% as.data.frame() %>%
  filter(grepl("MT-",gene_short_name) == TRUE)
mt_genes_id <- rownames(mt_genes)
colData(cds)$percent_mito <- 100 * (colSums(exprs(cds)[mt_genes_id,])/colSums(exprs(cds)))

ggplot(colData(cds) %>% as.data.frame() %>% arrange(desc(n_umi)) %>% mutate(cell_rank = dplyr::row_number()), aes(x= log10(cell_rank), y = log10(n_umi))) +
  geom_point(size = 0.5, stroke = 0, aes(color = cell_line)) +
  geom_hline(linewidth = 0.2, yintercept = log10(500)) +
  monocle3:::monocle_theme_opts() +
  scale_color_manual(name = "Cell Line",
                     breaks = c("BT112", "BT228", "BT333"),
                     values = c("darkblue", "darkred", "darkgreen")) +
  facet_wrap(~cell_line, nrow = 1) +
  scale_y_continuous(breaks = c(1,2,3,4,5)) +
  xlab("Log10 CellRank") +
  ylab("Log10 UMI") +
  guides(guides(color = guide_legend(override.aes = list(size=2)))) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom", 
        legend.key.width = unit(0.5,"line"),
        legend.key.height = unit(0.2,"line"),
        text = element_text(size = 6))
ggsave("QC/Kneeplot_500umi.png", dpi = 600, height = 1.5, width = 3)

ggplot(colData(cds) %>% as.data.frame() %>% separate(col = P5, into = c(NA, "PCR_plate"), sep = -1) %>% arrange(desc(n_umi)) %>% mutate(cell_rank = dplyr::row_number()), aes(x= log10(cell_rank), y = log10(n_umi))) +
  geom_point(size = 0.5, stroke = 0, aes(color = PCR_plate)) +
  geom_hline(size = 0.2, yintercept = log10(500)) +
  monocle3:::monocle_theme_opts() +
  # viridis::scale_color_viridis(name = "PCR Plate",
  #                              discrete = TRUE,
  #                              option = "D", direction = -1) +
  scale_color_manual(values = viridis::magma(n = 4)[2:3], name = "PCR Plate") +
  scale_y_continuous(breaks = c(1,2,3,4,5)) +
  scale_x_continuous(breaks = c(1,2,3,4,5,6)) +
  xlab("Log10 CellRank") +
  ylab("Log10 UMI") +
  guides(guides(color = guide_legend(override.aes = list(size=2)))) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom", 
        legend.key.width = unit(1,"line"),
        legend.key.height = unit(0.2,"line"),
        text = element_text(size = 6))
ggsave("QC/Kneeplot_PCRplate_500umi.png", dpi = 600, height = 2, width = 2)

ggplot(colData(cds) %>% as.data.frame() %>% filter(hash_umis >= 1)) +
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
ggsave("QC/umi_density_300umi.png", dpi = 600, height = 2, width = 2.5)

ggplot(colData(cds) %>% as.data.frame() %>% filter(n_umi >= 300)) +
  ggridges::geom_density_ridges(aes(x = log10(hash_umis), y = cell_line, fill = cell_line), 
                                alpha = 0.5, 
                                rel_min_height = 0.001,
                                position = position_nudge(y = -0.5),
                                show.legend = FALSE) +
  geom_vline(xintercept = log10(5)) +
  monocle3:::monocle_theme_opts() +
  xlab("Log10 UMI") +
  ylab("Density") +
  scale_x_continuous(breaks = seq(0,5,1)) +
  viridis::scale_fill_viridis(option = "B", discrete = TRUE)
ggsave("QC/hash_umi_density_300umifilter.png", dpi = 600, height = 2, width = 2.5)

ggplot(colData(cds) %>% as.data.frame() %>% filter(n_umi >= 100)) +
  ggridges::geom_density_ridges(aes(x = log10(hash_umis), y = cell_line, fill = cell_line), 
                                alpha = 0.5, 
                                rel_min_height = 0.001,
                                position = position_nudge(y = -0.5),
                                show.legend = FALSE) +
  geom_vline(xintercept = log10(10)) +
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
  geom_vline(xintercept = log10(3)) +
  monocle3:::monocle_theme_opts() +
  xlab("Log10 Top to Second Best Ratio") +
  ylab("Density") +
  viridis::scale_fill_viridis(option = "B", discrete = TRUE)
ggsave("QC/hash_ratio_density_100umi.png", dpi = 600, height = 2, width = 3)

ggplot(colData(cds) %>% as.data.frame() %>% arrange(desc(hash_umis)) %>% dplyr::mutate(rank = row_number()), 
       aes(x = log10(rank), y = log10(hash_umis))) +
  geom_point() 

ggplot(colData(cds) %>% as.data.frame() %>% filter(n_umi >= 300 & hash_umis >= 5),
       aes(x = cell_line, y = percent_mito)) +
  geom_violin(aes(fill = cell_line)) +
  monocle3:::monocle_theme_opts() +
  viridis::scale_fill_viridis(option = "B", discrete = TRUE)

ggplot(colData(cds) %>% as.data.frame() %>% filter(n_umi >= 300 & hash_umis >= 10 & top_to_second_best_ratio >= 3),
       aes(x = drug, y = log10_umi)) +
  geom_violin() +
  monocle3:::monocle_theme_opts() +
  viridis::scale_fill_viridis(option = "B", discrete = TRUE) +
  theme(axis.text.x = element_blank())

# ==============================================================================================
cds <- cds[,colData(cds)$n_umi >= 300]

# predicting and filtering doublets, print matrix file and run scrublet
Matrix::writeMM(t(exprs(cds)), file = "QC/UMI_count_filt.matrix")

# after running scrublet
scrublet_scores <- read_tsv("QC/scrublet/doublet_scores_EGFRi.txt",
                            col_names = F)

colData(cds)$doublet_score <- scrublet_scores$X1
ggplot(colData(cds) %>% as.data.frame, aes(x=doublet_score)) +
  geom_histogram(color = 'lightpink2', fill = 'lightpink2', binwidth = 0.02) +
  geom_vline(aes(xintercept=0.40), color = 'black') +
  annotate("text", label = "Filter = 0.40", x = .52, y = 9000) +
  monocle3:::monocle_theme_opts() +
  xlab("Doublet Score") +
  ylab("Count")
ggsave('QC/scrublet/doublet_histogram.png',
       height = 4, width = 6)

# Final QC cutoffs
cds <- cds[,colData(cds)$n_umi >= 300 & 
             colData(cds)$doublet_score <= 0.40 &
             !is.na(colData(cds)$hash_umis) &
             colData(cds)$hash_umis >= 5 & 
             colData(cds)$top_to_second_best_ratio >= 2.5]

cds <- estimate_size_factors(cds)
cds <- detect_genes(cds)

source("cell_cycle.R")
cc.genes <- readRDS("cc.genes.RDS")
cds <- estimate_cell_cycle(cds, g1s_markers = cc.genes$s.genes, g2m_markers = cc.genes$g2m.genes)

saveRDS(cds, file = "sci-Plex-EGFRi-7d-largescreen-cds.RDS")

# Final cell counts
n_umi_per_PCR_tot <- colData(cds) %>% as.data.frame() %>%
  dplyr::summarise(count = n(), mean_umi_per_cell = mean(n_umi), med_umi_per_cell = median(n_umi), total_umi = sum(n_umi)) %>%
  mutate(P7_row = "Total")

n_umi_per_PCR <- colData(cds) %>% as.data.frame() %>%
  separate(P7, into = c(NA, "P7_row"), sep = 2) %>%
  group_by(P7_row) %>%
  summarise(count = n(), mean_umi_per_cell = mean(n_umi), med_umi_per_cell = median(n_umi), total_umi = sum(n_umi)) %>%
  bind_rows(n_umi_per_PCR_tot)

test <- colData(cds) %>% as.data.frame() %>%
  group_by(cell_line, drug, dose) %>%
  summarise(count = n(), .groups = "drop")
  # group_by(cell_line) %>%
  summarise(average = mean(count), median = median(count)) 

n_umi_per_PDCL <- colData(cds) %>% as.data.frame() %>%
  group_by(cell_line) %>%
  summarise(count = n(), mean_umi_per_cell = mean(n_umi), med_umi_per_cell = median(n_umi), total_umi = sum(n_umi))
