library(tidyverse)
library(monocle3)
library(ComplexHeatmap)
this.dir <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(this.dir))
rm(this.dir)

diff_test_result <- read_csv("supplementary_table_S1_pilot_screen_dose_DEG_fit_models_signif_result.csv")

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
