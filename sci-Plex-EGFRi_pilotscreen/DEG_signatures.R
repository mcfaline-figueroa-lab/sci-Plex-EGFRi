library(tidyverse)
library(monocle3)
library(ComplexHeatmap)
this.dir <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(this.dir))
rm(this.dir)

diff_test_result <- read_csv("supplementary_table_S1_pilot_screen_dose_DEG_fit_models_result.csv")

DEG_filtered <- diff_test_result %>%
  filter(grepl("log", term) & q_value < 0.001 & abs(normalized_effect)>0.05) %>%
  pull(id) %>%
  unique()

upset.list <- list()
for (drug in c("Afatinib", "Brigatinib", "CUDC-101", "Neratinib", "Osimertinib")) {
  temp <- diff_test_result %>%
    filter(term == "log(dose + 0.1)", 
           treatment == drug, 
           id %in% DEG_filtered, 
           normalized_effect > 0, 
           q_value < 0.001, 
           !is.na(normalized_effect)) %>%
    select(gene_short_name) %>%
    distinct() %>%
    pull()
  upset.list[[drug]] <- temp
}

intersection <- ComplexHeatmap::make_comb_mat(upset.list, mode = "distinct")

intersect_plot <- ComplexHeatmap::UpSet(intersection[comb_size(intersection) >= 30], 
                                   comb_order = order(comb_size(intersection[comb_size(intersection) >= 30])),
                                   top_annotation = upset_top_annotation(intersection[comb_size(intersection) >= 30],  
                                                                         # ylim = c(0, 500), 
                                                                         add_numbers = TRUE, 
                                                                         height = unit(1, "cm"), 
                                                                         annotation_name_gp = gpar(fontsize = 3), show_annotation_name = F),
                                   right_annotation = upset_right_annotation(intersection[comb_size(intersection) >= 30], 
                                                                             add_numbers = FALSE, 
                                                                             width = unit(1, "cm"), 
                                                                             annotation_name_gp = gpar(fontsize = 3), show_annotation_name = F
                                                                             )
)
intersect_plot

intersection <- intersection[comb_size(intersection) >= 30]
gene_sets <- sapply(comb_name(intersection), function(nm) extract_comb(intersection, nm))
test <- set_name(intersection)
test_names <- data.frame(x = names(gene_sets)) %>%
  separate(x, into = c("1","2","3","4","5"), sep = c(1,2,3,4)) %>%
  mutate(`1` = ifelse(`1` == 1, "Afat",""),
         `2` = ifelse(`2` == 1, "Brigat",""),
         `3` = ifelse(`3` == 1, "CUDC",""),
         `4` = ifelse(`4` == 1, "Nerat",""),
         `5` = ifelse(`5` == 1, "Osimert","")) %>%
  unite(col = "name", sep = "/", 1:5) %>%
  mutate(name = paste("_", name, "_", sep = "")) %>%
  mutate(name = str_replace_all(name, c("////" = "/", "///" = "/", "//" = "/"))) %>%
  mutate(name = str_replace_all(name, c("_/" = "", "/_" = "", "_" = ""))) %>%
  pull()
names(gene_sets) <- test_names

saveRDS(gene_sets, "pilot_screen_EGFRi_intersection_sigantures.RDS")


# ==============================================================================
# Correlation coefficient 
# ==============================================================================

deg_list <- diff_test_result %>% 
  filter(grepl(pattern = "log", term), 
         q_value < 0.01) %>%
  pull(id) %>%
  unique()

deg_precorr <- list()
deg_corr <- list()
for (cell_type in c("A172", "T98G", "U87MG", "BT333")){
  
  deg_precorr[[cell_type]] <- diff_test_result %>%
    filter(Cell.Line == cell_type) %>%
    filter(grepl("log", term),id %in% deg_list) %>%
    mutate(condition = paste0(Cell.Line,"_",treatment)) %>%
    dplyr::select(condition, id, normalized_effect) %>%
    tidyr::spread(key = condition, value = normalized_effect) %>%
    as.data.frame()
  
  row.names(deg_precorr[[cell_type]]) <- deg_precorr[[cell_type]]$id
  deg_precorr[[cell_type]]$id <- NULL
  
  deg_corr[[cell_type]] <- cor(deg_precorr[[cell_type]], method = "pearson")
 
}

deg_corr_df <- do.call("cbind", deg_corr)
row.names(deg_corr_df) <- sapply(row.names(deg_corr_df),function(x){stringr::str_split(x, pattern = "_")[[1]][2]})
saveRDS(deg_corr_df, "DEG_correlation_matrtix.RDS")

gene_meta <- data.frame(id = c("Afatinib", "Brigatinib", "CUDC-101", "EAI045", "Neratinib", "Osimertinib"), gene_short_name = NA)
row.names(gene_meta) <- gene_meta$id
cell_meta <- data.frame(id = colnames(deg_corr_df)) %>%
  separate(id, into = c("Cell.Line", "treatment"), sep = "_", remove = F)
row.names(cell_meta) <- cell_meta$id

cds_corr <- new_cell_data_set(deg_corr_df,
                              gene_metadata = gene_meta,
                              cell_metadata = cell_meta)

reducedDims(cds_corr)[["PCA"]] <- t(exprs(cds_corr))

cds_corr <- reduce_dimension(cds_corr, umap.min_dist = 0.1, umap.n_neighbors = 10)
colData(cds_corr)$UMAP1 <- reducedDims(cds_corr)[["UMAP"]][,1]
colData(cds_corr)$UMAP2 <- reducedDims(cds_corr)[["UMAP"]][,2]

colData(cds_corr) %>% 
  as.data.frame %>% 
  mutate(Cell.Line = factor(Cell.Line, levels = c("A172", "T98G", "U87MG", "BT333"))) %>%
  mutate(treatment = factor(treatment, levels = c("EAI045", "Afatinib", "Neratinib", "Osimertinib", "Brigatinib", "CUDC-101"))) %>%
  ggplot() +
  geom_point(aes(x = UMAP1, y = UMAP2, shape = Cell.Line, fill = treatment,color = treatment), size = 2.5, stroke = 0.5, alpha = 0.9) +
  monocle3:::monocle_theme_opts() +
  guides(guides(fill = guide_legend(override.aes = list(size=2)),
                shape ="none")) +
  theme(
    legend.position = "bottom",
    text = element_text(size = 10),
    legend.key.width = unit(0.2,"line"),
    legend.key.height = unit(0.75,"line"),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 10)) +
  scale_shape_manual(name = "Cell Line", values = c(21,22,24,25)) +
  scale_fill_manual(name = "Drug", values = c(viridis::magma(n = 6)[1:6]), 
                    aesthetics = c("fill", "color"))

ggsave(filename = "covariate/covariate_UMAPs/UMAP_correlation_coeff_legend.png", non, width = 3, height = 2, dpi = 600)

