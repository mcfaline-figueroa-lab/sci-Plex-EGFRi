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
  separate(oligo, into = c("drug", "dose", "hash"), sep = "_") %>%
  separate(dose, sep = -2, into = c("dose", NA), convert = TRUE) %>%
  separate(hash, sep = 2, into = c("hash_plate", "hash_well")) %>%
  mutate(dose = case_when(
    cell_line == "BT333" & replicate == "1" & hash_plate == "02" & hash_well == "H01" ~ 10,
    cell_line == "BT333" & replicate == "1" & hash_plate == "03" & hash_well == "H01" ~ 100,
    TRUE ~ dose
  )) %>%
  mutate(drug = case_when(
    cell_line == "BT228" & replicate == "2" & hash_plate == "01" & hash_well == "C01" ~ "AG-490",
    cell_line == "BT228" & replicate == "2" & hash_plate == "01" & hash_well == "C02" ~ "Pelitinib",
    cell_line == "BT228" & replicate == "2" & hash_plate == "01" & hash_well == "C03" ~ "ChrysophanicAcid",
    cell_line == "BT228" & replicate == "2" & hash_plate == "01" & hash_well == "C04" ~ "TAK-285",
    cell_line == "BT228" & replicate == "2" & hash_plate == "01" & hash_well == "C05" ~ "AG555",
    cell_line == "BT228" & replicate == "2" & hash_plate == "01" & hash_well == "C06" ~ "AZ5104",
    cell_line == "BT228" & replicate == "2" & hash_plate == "01" & hash_well == "C07" ~ "EAI045",
    cell_line == "BT228" & replicate == "2" & hash_plate == "01" & hash_well == "C08" ~ "Cyasterone",
    cell_line == "BT228" & replicate == "2" & hash_plate == "01" & hash_well == "C09" ~ "EGFRInhibitor",
    cell_line == "BT228" & replicate == "2" & hash_plate == "01" & hash_well == "C10" ~ "JBJ-04-125-02",
    cell_line == "BT228" & replicate == "2" & hash_plate == "01" & hash_well == "D01" ~ "AC480",
    cell_line == "BT228" & replicate == "2" & hash_plate == "01" & hash_well == "D02" ~ "Genistein",
    cell_line == "BT228" & replicate == "2" & hash_plate == "01" & hash_well == "D03" ~ "Epigallocatechin-Gallate",
    cell_line == "BT228" & replicate == "2" & hash_plate == "01" & hash_well == "D04" ~ "Varlitinib",
    cell_line == "BT228" & replicate == "2" & hash_plate == "01" & hash_well == "D05" ~ "AG494",
    cell_line == "BT228" & replicate == "2" & hash_plate == "01" & hash_well == "D06" ~ "Osimertinib",
    cell_line == "BT228" & replicate == "2" & hash_plate == "01" & hash_well == "D07" ~ "Butein",
    cell_line == "BT228" & replicate == "2" & hash_plate == "01" & hash_well == "D08" ~ "HS-10296",
    cell_line == "BT228" & replicate == "2" & hash_plate == "01" & hash_well == "D09" ~ "RG13022",
    cell_line == "BT228" & replicate == "2" & hash_plate == "01" & hash_well == "D10" ~ "AG-1557",
    TRUE ~ drug
  )) %>%
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
  xlab("Log10 CellRank") +
  ylab("Log10 UMI") +
  guides(guides(color = guide_legend(override.aes = list(size=2)))) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom", 
        legend.key.width = unit(1,"line"),
        legend.key.height = unit(0.2,"line"),
        text = element_text(size = 6))
ggsave("QC/Kneeplot_PCRplate_500umi.png", dpi = 600, height = 2, width = 2)

saveRDS(cds, "BTline_EGFRi_24hr_unfiltered.RDS")

dim(cds[,colData(cds)$n_umi >= 500 & !is.na(colData(cds)$hash_umis) & colData(cds)$hash_umis >= 5 & colData(cds)$top_to_second_best_ratio >= 2.5])
dim(cds[,colData(cds)$n_umi >= 500])

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
cds <- readRDS("BTline_EGFRi_24hr_unfiltered.RDS")

cds <- cds[,colData(cds)$n_umi >= 300]

# predicting and filtering doublets
Matrix::writeMM(t(exprs(cds)), file = "QC/UMI_count_filt.matrix") # input matrix into scrublet

scrublet_scores <- read_tsv("QC/scrublet/doublet_scores_EGFRi.txt",
                            col_names = F)

colData(cds)$doublet_score <- scrublet_scores$X1
ggplot(colData(cds) %>% as.data.frame, aes(x=doublet_score)) +
  geom_histogram(color = 'lightpink2', fill = 'lightpink2', binwidth = 0.02) +
  geom_vline(aes(xintercept=0.50), color = 'black') +
  annotate("text", label = "Filter = 0.50", x = .62, y = 9000) +
  monocle3:::monocle_theme_opts() +
  xlab("Doublet Score") +
  ylab("Count")
ggsave('QC/scrublet/doublet_histogram.png',
       height = 4, width = 6)

# Final QC cutoffs
cds <- cds[,colData(cds)$n_umi >= 500 & 
             colData(cds)$doublet_score <= 0.50 &
             !is.na(colData(cds)$hash_umis) &
             colData(cds)$hash_umis >= 5 & 
             colData(cds)$top_to_second_best_ratio >= 2.5]

cds <- estimate_size_factors(cds)
cds <- detect_genes(cds)

source("cell_cycle.R")
cc.genes <- readRDS("cc.genes.RDS")
cds <- estimate_cell_cycle(cds, g1s_markers = cc.genes$s.genes, g2m_markers = cc.genes$g2m.genes)

saveRDS(cds, file = "cds_large_screen.RDS")

# ==========================================================================================
# Getting cell counts, median umi, median hash umi, unique genes per cell line
# ==========================================================================================

cell_count_metrics_tot <- colData(cds) %>% as.data.frame() %>%
  summarise(cell_count = n(), 
            median_umi = median(n_umi),
            median_hash_umi = median(hash_umis),
            mean_express_genes = mean(num_genes_expressed)) %>%
  mutate(cell_line = "Total") %>%
  select(cell_line, everything())

cell_count_metrics <- colData(cds) %>% as.data.frame() %>%
  group_by(cell_line) %>%
  summarise(cell_count = n(), 
            median_umi = median(n_umi),
            median_hash_umi = median(hash_umis),
            mean_express_genes = mean(num_genes_expressed))

cell_count_metrics_per_cell_drug_dose <- colData(cds) %>% as.data.frame() %>%
  group_by(cell_line, drug, dose) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(cell_line) %>%
  summarise(mean_drug_dose_count = mean(cell_count))

cell_count_metrics_per_drug_dose_tot <- colData(cds) %>% as.data.frame() %>%
  group_by(cell_line, drug, dose) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  summarise(mean_drug_dose_count = mean(cell_count)) %>%
  mutate(cell_line = "Total") %>%
  select(cell_line, everything())

#========================================================================================
# Correlating replicates for QC
#========================================================================================

expressed_genes <- row.names(rowData(cds)[Matrix::rowSums(exprs(cds) > 0) > 
                                            dim(cds)[2]*0.05 ,])

saveRDS(expressed_genes, file = "large_screen_expressed_genes_5percent_allcells.RDS")

# by line
for (line in c("BT112", "BT228", "BT333")) {
  for (drug in c("Osimertinib")) { #input drug of interest
    cds_temp <- cds[expressed_genes,
                    colData(cds)$cell_line == line & 
                      colData(cds)$drug == drug &
                      colData(cds)$log10_dose == 4 &
                      colData(cds)$replicate == 1]
    rep1 <- data.frame("rep_1" = Matrix::rowMeans(normalized_counts(cds_temp))) %>%
      rownames_to_column(var = "gene")
    cds_temp <- cds[expressed_genes,
                    colData(cds)$cell_line == line & 
                      colData(cds)$drug == drug & 
                      colData(cds)$log10_dose == 4 &
                      colData(cds)$replicate == 2]
    rep2 <- data.frame("rep_2" = Matrix::rowMeans(normalized_counts(cds_temp))) %>%
      rownames_to_column(var = "gene")
    
    corr_df <- left_join(rep1, rep2, by = c("gene" = "gene")) %>%
      column_to_rownames(var = "gene")
    corr <- cor(x = corr_df, method = "pearson")
    
    temp <- ggplot(data = corr_df, aes(x = rep_1, y = rep_2)) +
      geom_point(size = 0.2) +
      geom_smooth(method = "lm", 
                  linewidth = 0.3) +
      annotate(geom = "text",
               x = (min(corr_df$rep_1) + 3 * sd(corr_df$rep_1)),
               y = max(corr_df$rep_2),
               label = paste("r =", round(corr[1,2], 3)),
               size = 2) +
      xlab("Replicate 1") +
      ylab("Replicate 2") +
      theme(text = element_text(size = 6)) +
      monocle3:::monocle_theme_opts()
    ggsave(plot = temp, filename = paste("QC/replicate_correlation/normalized_counts_replicate_corr_", line, "_", drug, ".png", sep = ""),
           dpi = 600, height = 2, width = 2)
    
  }
}


# for all PDCLs
for (drug in c("Osimertinib")) {
  cds_temp <- cds[expressed_genes,
                  colData(cds)$drug == drug &
                    colData(cds)$log10_dose == 4 &
                    colData(cds)$replicate == 1]
  rep1 <- data.frame("rep_1" = Matrix::rowMeans(normalized_counts(cds_temp))) %>%
    rownames_to_column(var = "gene")
  cds_temp <- cds[expressed_genes,
                    colData(cds)$drug == drug &
                    colData(cds)$replicate == 2]
  rep2 <- data.frame("rep_2" = Matrix::rowMeans(normalized_counts(cds_temp))) %>%
    rownames_to_column(var = "gene")
  
  corr_df <- left_join(rep1, rep2, by = c("gene" = "gene")) %>%
    column_to_rownames(var = "gene")
  corr <- cor(x = corr_df, method = "pearson")
  
  temp <- ggplot(data = corr_df, aes(x = rep_1, y = rep_2)) +
    geom_point(size = 0.2) +
    geom_smooth(method = "lm", 
                linewidth = 0.3) +
    annotate(geom = "text",
             x = (min(corr_df$rep_1) + 4.5 * sd(corr_df$rep_1)),
             y = max(corr_df$rep_2),
             label = paste("r =", round(corr[1,2], 3)),
             size = 2) +
    xlab("Replicate 1") +
    ylab("Replicate 2") +
    ggtitle(drug) +
    theme(text = element_text(size = 6),
          plot.title = element_text(hjust = 0.5)) +
    monocle3:::monocle_theme_opts()
  ggsave(plot = temp, filename = paste("QC/replicate_correlation/normalized_counts_allPDCL_highdose_replicate_corr_", drug, ".png", sep = ""),
         dpi = 600, height = 1.2, width = 1.2)
  
}

# for all drugs
for (line in c("BT112", "BT228", "BT333")) {
  cds_temp <- cds[expressed_genes,
                  colData(cds)$cell_line == line &
                    colData(cds)$replicate == 1]
  rep1 <- data.frame("rep_1" = Matrix::rowMeans(normalized_counts(cds_temp))) %>%
    rownames_to_column(var = "gene")
  cds_temp <- cds[expressed_genes,
                  colData(cds)$drug == drug &
                    colData(cds)$replicate == 2]
  rep2 <- data.frame("rep_2" = Matrix::rowMeans(normalized_counts(cds_temp))) %>%
    rownames_to_column(var = "gene")
  
  corr_df <- left_join(rep1, rep2, by = c("gene" = "gene")) %>%
    column_to_rownames(var = "gene")
  corr <- cor(x = corr_df, method = "pearson")
  
  temp <- ggplot(data = corr_df, aes(x = rep_1, y = rep_2)) +
    geom_point(size = 0.2) +
    geom_smooth(method = "lm", 
                linewidth = 0.3) +
    annotate(geom = "text",
             x = (min(corr_df$rep_1) + 3 * sd(corr_df$rep_1)),
             y = max(corr_df$rep_2),
             label = paste("r =", round(corr[1,2], 3)),
             size = 2) +
    xlab("Replicate 1") +
    ylab("Replicate 2") +
    theme(text = element_text(size = 6)) +
    monocle3:::monocle_theme_opts()
  ggsave(plot = temp, filename = paste("QC/replicate_correlation/normalized_counts_alldrug_replicate_corr_", line, ".png", sep = ""),
         dpi = 600, height = 2, width = 2)
  
}



