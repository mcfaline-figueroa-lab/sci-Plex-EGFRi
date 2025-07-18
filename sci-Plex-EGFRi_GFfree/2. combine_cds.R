library(tidyverse)
library(monocle3)
this.dir <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(this.dir))
rm(this.dir)

cds_1 <- readRDS("sci-Plex-EGFRi-GFfree-cds-P7-A-P5-10.RDS")
cds_2 <- readRDS("sci-Plex-EGFRi-GFfree-cds-P7-C-P5-7.RDS")

cds <- combine_cds(cds_list = list(cds_1, cds_2), cell_names_unique = T)

cds <- detect_genes(cds)

summary <- colData(cds) %>% as.data.frame() %>%
  dplyr::group_by(cell_line, GF_condition, drug) %>%
  summarise(count = n(), mean_umi = mean(n_umi))

# printing metrics
n_umi_per_PCR_tot <- colData(cds) %>% as.data.frame() %>%
  group_by(cell_line) %>%
  summarise(count = n(), mean_umi_per_cell = mean(n_umi), med_umi_per_cell = median(n_umi), 
            total_umi = sum(n_umi), median_hash_umi = median(hash_umis))
write_tsv(n_umi_per_PCR_tot, "QC/combined/umis_post_QC_filtering.txt")

test <- colData(cds) %>% as.data.frame() %>%
  group_by(cell_line, GF_condition, drug, dose) %>%
  dplyr::summarise(count = n(), .groups = "drop") %>%
  group_by(cell_line) %>%
  summarise(mean_cell_count = round(mean(count),3), median_cell_count = round(median(count), 3)) 
write_tsv(test, "QC/combined/cell_counts_post_QC_filtering.txt")

test <- colData(cds) %>% as.data.frame() %>%
  group_by(cell_line, GF_condition, drug) %>%
  dplyr::summarise(count = n(), .groups = "drop") %>%
  # group_by(cell_line) %>%
  summarise(total_cell_count = sum(count),
            mean_cell_count = round(mean(count),3), 
            median_cell_count = round(median(count), 3)) 
write_tsv(test, "QC/combined/cell_counts_per_PDCL_GF_drug_post_QC_filtering.txt")

saveRDS(cds, "sci-Plex-EGFRi-GFfree-cds.RDS")