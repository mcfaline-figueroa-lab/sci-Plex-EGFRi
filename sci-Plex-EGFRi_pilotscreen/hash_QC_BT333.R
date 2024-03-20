# ==============================================================================
# Adding hash information to cds
# ==============================================================================
library(tidyverse)
library(monocle3)

cds.pre <- readRDS("cds_precell_prehash.RDS")

# Assigning cell line to colData

cols <- colData(cds.pre) %>% as.data.frame() %>%
  mutate(cell_ID=Cell) %>%
  tidyr::separate(col = Cell, into = c('P5', 'P7', NA, NA, 'n1', NA, NA, 'n2'), sep = "_") %>%
  mutate(RT = paste('RT_BC_', n1, sep = '')) %>%
  mutate(Lig = paste('Lig_BC_', n2, sep = '')) %>%
  select(-n1, -n2)

cell.line.map <- read_csv(file = cell.line.map.path)
RT.barcodes <- read.table(file = RT.barcodes.path, 
                          col.names = c("RT", "Barcode"))

cell.line.map <- left_join(cell.line.map, RT.barcodes, by = c("Barcode" = "Barcode"))

cols <- cols %>%
  left_join(cell.line.map, by = c("RT" = "RT")) %>%
  select(-`Well Position`, -Barcode) %>%
  dplyr::rename(Cell.Line = `Cell Line`)

# Assigning hash to colData
hashTable <- read.table(hash.path, 
                        col.names = c("sample", "cell_ID", "oligo", "axis", "hash_umis"))

hashTable_summary <- hashTable %>%
  group_by(cell_ID) %>%
  mutate(proportion = hash_umis/sum(hash_umis),total_hash_umis_per_cell = sum(hash_umis)) %>%
  arrange(desc(proportion)) %>%
  mutate(rank = row_number()) %>%
  dplyr::filter(rank %in% c(1,2)) %>%
  mutate(top_to_second_best_ratio = ifelse(sum(rank) > 1, proportion[rank == 1]/proportion[rank == 2], 1)) %>%
  dplyr::filter(rank == 1)

hashTable_summary$treatment <- sapply(hashTable_summary$oligo,function(x){stringr::str_split(x,pattern = "_")[[1]][1]})
hashTable_summary$dose <- sapply(hashTable_summary$oligo,function(x){
  y <- stringr::str_split(x,pattern = "_")[[1]][2]
  z <- stringr::str_split(y,pattern = "nM")[[1]][1]
  z <- as.numeric(z)
  return(z)})
hashTable_summary$well_ID <- sapply(hashTable_summary$oligo,function(x){stringr::str_split(x,pattern = "_")[[1]][3]})

hashTable_summary <- hashTable_summary %>% 
  mutate(well_ID_temp = well_ID) %>%
  tidyr::separate(well_ID_temp, into = c("plate", "plate_position"), sep = -3) %>%
  mutate(plate_coor = plate_position) %>%
  tidyr::separate(plate_coor, into = c("plate_row", "plate_col"), sep = -2)

cols <- cols %>% 
  left_join(hashTable_summary, by = c("cell_ID" = "cell_ID", "sample" = "sample")) %>%
  select(-Size_Factor, -axis, -proportion, -rank)

rownames(cols) <- cols$cell_ID

# assiging our cols dataframe as colData(cds)
colData(cds.pre) <- DataFrame(cols)

# adding columns for log10(n.umi) and % mito
colData(cds.pre)$log10.umi <- colData(cds.pre)$n.umi %>% log10()

mt_genes <- rowData(cds.pre) %>% as.data.frame() %>%
  filter(grepl("^MT-",gene_short_name) == TRUE)
mt_genes_id <- rownames(mt_genes)
colData(cds.pre)$percent_mito <- 100 * (colSums(exprs(cds.pre)[mt_genes_id,])/colSums(exprs(cds.pre)))

saveRDS(cds.pre, 
        "cds_hash.RDS")

# ==============================================================================
# QC 
# ==============================================================================
library(tidyverse)
library(monocle3)

# Load in 3 sequencing run cds
cds <- readRDS("cds_hash.RDS")

# Plot ridgeplot (density plot) of n.umi for each cell line
ggplot(data = as.data.frame(colData(cds)) %>% filter(is.na(Cell.Line) == FALSE)) +
  ggridges::geom_density_ridges(aes(x = log10(n.umi), y = Cell.Line, fill = Cell.Line), 
                                alpha = 0.5, 
                                rel_min_height = 0.001,
                                position = position_nudge(y = -0.5),
                                show.legend = FALSE) +
  geom_vline(xintercept = log10(100)) +
  annotate(geom = "text", label = "Cutoff = 100", x = log10(330), y = 1.4) + 
  monocle3:::monocle_theme_opts() +
  xlab("Log10(UMI)") +
  ylab("Density") +
  viridis::scale_fill_viridis(option = "C", discrete = TRUE)
  #scale_x_continuous(limits = c(0,3.5), breaks = c(0,1,2,3))
ggsave("Log10_UMI_density_plot_100umi.png", dpi = 600, height = 2, width = 5)

ggplot(colData(cds) %>% as.data.frame() %>% arrange(desc(n.umi)) %>% mutate(cell_rank = dplyr::row_number()), aes(x= log10(cell_rank), y = log10(n.umi))) +
  geom_point(size = 0.5, stroke = 0, aes(color = "BT333")) +
  geom_hline(size = 0.2, yintercept = log10(100)) +
  annotate(geom = "text", label = "Cutoff = 100", x = 0.5, y = log10(130)) +
  monocle3:::monocle_theme_opts() +
  viridis::scale_color_viridis(name = "Cell Line", 
                              breaks = c("BT333"),
                              discrete = TRUE,
                              option = "C") +
  xlab("Log10 CellRank") +
  ylab("Log10 UMI") +
  guides(guides(color = guide_legend(override.aes = list(size=2))))
ggsave("Kneeplot_cellline_100umi.png", dpi = 600, height = 4, width = 6)


# Preparing cds for scrublet doublet detection
cds <- cds[,!(is.na(colData(cds)$hash_umis)) & 
                    !(is.na(colData(cds)$Cell.Line)) &
                    colData(cds)$n.umi >= 100]
test <- exprs(cds)
test <- t(test)
Matrix::writeMM(test, "UMI_count_filt.matrix") # import this file into scrublet .ipynb to get doublet_scores

scores <- read.table(file = "QC/doublet_scores_EGFRi.txt",
                     header = F)
colData(cds)$doublet_score <- scores$V1

ggplot(colData(cds) %>% as.data.frame, aes(x=doublet_score)) +
  geom_histogram(color = 'lightpink2', fill = 'lightpink2', binwidth = 0.02) +
  geom_vline(aes(xintercept=0.50), color = 'black') +
  annotate("text", label = "Filter = 0.50", x = .62, y = 9000) +
  monocle3:::monocle_theme_opts() +
  xlab("Doublet Score") +
  ylab("Count")
ggsave('doublet_histogram.png',
       height = 4, width = 6)

cds <- cds[, colData(cds)$doublet_score <= 0.50]
dim(cds)

ggplot(data = as.data.frame(colData(cds)) %>% filter(is.na(Cell.Line) == FALSE)) +
  ggridges::geom_density_ridges(aes(x = top_to_second_best_ratio, y = Cell.Line, fill = Cell.Line), 
                                alpha = 0.5, 
                                rel_min_height = 0.001,
                                position = position_nudge(y = -0.5),
                                show.legend = FALSE) +
  geom_vline(xintercept = 2) +
  annotate(geom = "text", label = "Cutoff = 2", x = 3, y = 1.4) + 
  monocle3:::monocle_theme_opts() +
  xlab("Top to second best ratio") +
  ylab("Density") +
  viridis::scale_fill_viridis(option = "C", discrete = TRUE) +
  scale_x_continuous(limits = c(0,10))
ggsave("Log10_UMI_density_plot_t2sbr2.png", dpi = 600, height = 2, width = 5)

cds <- cds[, colData(cds)$hash_umis >= 10 &
             colData(cds)$top_to_second_best_ratio >= 5]


# Saving CDS as RDS that has been filtered by umi, hash, and doublet
cds <- estimate_size_factors(cds)
cds <- detect_genes(cds)
colData(cds)$log10_dose <- log10(colData(cds)$dose)
colData(cds)$cell <- colData(cds)$cell_ID
saveRDS(object = cds, file = "cds_QC_BT333.RDS")
