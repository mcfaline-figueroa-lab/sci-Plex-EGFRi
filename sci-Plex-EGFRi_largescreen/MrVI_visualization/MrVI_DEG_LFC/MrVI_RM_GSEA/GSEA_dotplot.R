library(monocle3)
library(tidyverse)
library(ComplexHeatmap)

this.dir <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(this.dir))
rm(this.dir)

source("GSA_helper_functions.r")
source("loadGSCSafe.R")

hallmarksGSC <- loadGSCSafe(file = "h.all.v6.0.symbols.gmt",
                            type = "gmt")
oncogenicGSC <- loadGSCSafe(file = "c6.all.v2023.2.Hs.symbols.gmt",
                            type = "gmt")
kinaseupGSC <- loadGSCSafe(file = "gene_set_library_up_crisp.gmt",
                           type = "gmt")
kinasedownGSC <- loadGSCSafe(file = "gene_set_library_dn_crisp.gmt",
                             type = "gmt")

gsc_list <- list("hallmarksGSC" = hallmarksGSC, 
                 "oncogenicGSC" = oncogenicGSC, 
                 "kinaseupGSC" = kinaseupGSC, 
                 "kinasedownGSC" = kinasedownGSC)
rm(hallmarksGSC, oncogenicGSC, kinaseupGSC, kinasedownGSC)

cds <- readRDS("cds_large_screen.RDS")
gene_universe_ensembl <- row.names(rowData(cds)[Matrix::rowSums(exprs(cds) > 0) > 
                                                                     dim(cds)[2]*0.05 ,])

gene_universe <- rowData(cds) %>%
  as.data.frame() %>%
  filter(id %in% gene_universe_ensembl) %>%
  select(gene_short_name) %>%
  distinct() %>%
  pull()

GSA_full_up <- data.frame()
suppressMessages(
  for(line in c("bt228")) {
    for (DRC_temp in 1:20) {
      top_genes_temp <- read_csv(paste0("DRC_up_signature_LFCandQVal_Filtered/",line,"_DRC", DRC_temp,"_up_signature_02062024.csv")) 
      
      if (length(top_genes_temp) > 1) {
        top_genes <- top_genes_temp %>%
          select(2) %>%
          pull()
        
        for(set in c("hallmarksGSC", "oncogenicGSC", "kinaseupGSC", "kinasedownGSC")){
          test <- piano::runGSAhyper(genes = top_genes, 
                                     gsc = gsc_list[[set]], 
                                     adjMethod = "none", 
                                     universe = gene_universe)
          GSAhyper_df <- as.data.frame(test$p.adj) %>%
            rownames_to_column(var = "gene_set")
          colnames(GSAhyper_df) <- c("gene_set","p_value")
          GSAhyper_df <- GSAhyper_df %>%
            mutate(cell_line = toupper(line),
                   DRC = DRC_temp,
                   gsc_set = set,
                   size = length(top_genes))
          GSA_full_up <- bind_rows(GSA_full_up, GSAhyper_df)
        }
      }
    }
  }
)
GSA_full_up$q_value <- p.adjust(GSA_full_up$p_value, method = "BH")
saveRDS(GSA_full_up, "GSEA_results/BT228_up_df_LFCandQval_filt.RDS")

GSA_full_up <- readRDS("GSEA_results/BT228_up_df_LFCandQval_filt.RDS")

GSA_full_up_plot <- GSA_full_up %>%
  group_by(DRC) %>%
  arrange(q_value) %>%
  mutate(rank = row_number()) %>%
  filter(rank <= 3) %>%
  filter(q_value <= 0.05)

ggplot(GSA_full_up_plot %>% arrange(gsc_set) %>% 
         filter(gene_set %in% c("MAP2K4_knockdown_62_GSE19091",
                                "HALLMARK_MTORC1_SIGNALING",
                                "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                                "ALK_DN.V1_UP",
                                "RAF_UP.V1_DN",
                                "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                                "IL2_UP.V1_DN",
                                "E2F1_UP.V1_DN",
                                "RAF_UP.V1_DN",
                                "KRAS.KIDNEY_UP.V1_UP",
                                "BRAF_overexpression_180_GSE46801",
                                "EGFR_drugactivation_19_GDS2146")) %>%
         mutate(gene_set = factor(gene_set, levels = gene_set)), 
       aes(x = as.factor(DRC), y = gene_set)) +
  geom_point(aes(color = -log10(q_value), size = -log10(q_value)), stroke = 0.1) +
  scale_color_gradient2(low = "white", high = "red", midpoint = 0, 
                        name = "-Log10 FDR",
                        guide = "legend") +
  scale_size_continuous(range = c(1,3)) +
  xlab("Response Module") +
  guides(size = guide_legend(title = "-Log10 FDR", nrow = 1)) +
  theme(text = element_text(size = 6),
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.2,"cm"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.box.margin = unit(c(-3,50,-5,0), "mm"),
        # axis.text.y = element_text(angle = 45, hjust = 1),
        axis.title.y = element_blank(),
        axis.ticks = element_line(linewidth = unit(0.2, "mm"))) +
  monocle3:::monocle_theme_opts()
ggsave("GSEA_results/figures/BT228_DRC_up_GSEA_top2_FDR005_chosen_display_long.png",
       dpi = 600, height = 2, width = 4)

ggplot(GSA_full_up_plot %>% arrange(gsc_set) %>%
         mutate(gene_set = factor(gene_set, levels = gene_set)), 
       aes(x = as.factor(DRC), y = gene_set)) +
  geom_point(aes(color = -log10(q_value), size = -log10(q_value)), stroke = 0.1) +
  scale_color_gradient2(low = "white", high = "red", midpoint = 0, 
                        name = "-Log10 FDR",
                        guide = "legend") +
  scale_size_continuous(range = c(1,3)) +
  xlab("Response Module") +
  guides(size = guide_legend(title = "-Log10 FDR", nrow = 1)) +
  theme(text = element_text(size = 6),
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.2,"cm"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.box.margin = unit(c(-3,50,-3,0), "mm"),
        # axis.text.y = element_text(angle = 45, hjust = 1),
        axis.title.y = element_blank(),
        axis.ticks = element_line(linewidth = unit(0.2, "mm"))) +
  monocle3:::monocle_theme_opts()

ggsave("GSEA_results/figures/supplemental/BT228_DRC_up_GSEA_top3_FDR005_chosen_display_long.png",
       dpi = 600, height = 3, width = 3.5)

GSA_full_down <- data.frame()
suppressMessages(
  for(line in c("bt228")) {
    for (DRC_temp in 1:20) {
      top_genes_temp <- read_csv(paste0("DRC_down_signature_LFCandQVal_Filtered/",line,"_DRC", DRC_temp,"_down_signature_02062024.csv")) 
      
      if (length(top_genes_temp) > 1) {
        top_genes <- top_genes_temp %>%
          select(2) %>%
          pull()
        
        for(set in c("hallmarksGSC", "oncogenicGSC", "kinaseupGSC", "kinasedownGSC")){
          test <- piano::runGSAhyper(genes = top_genes, 
                                     gsc = gsc_list[[set]], 
                                     adjMethod = "none", 
                                     universe = gene_universe)
          GSAhyper_df <- as.data.frame(test$p.adj) %>%
            rownames_to_column(var = "gene_set")
          colnames(GSAhyper_df) <- c("gene_set","p_value")
          GSAhyper_df <- GSAhyper_df %>%
            mutate(cell_line = toupper(line),
                   DRC = DRC_temp,
                   gsc_set = set,
                   size = length(top_genes))
          GSA_full_down <- bind_rows(GSA_full_down, GSAhyper_df)
        }
      }
    }
  }
)
GSA_full_down$q_value <- p.adjust(GSA_full_down$p_value, method = "BH")
saveRDS(GSA_full_down, "GSEA_results/BT228_down_df_LFCandQval_filt.RDS")

GSA_full_down <- readRDS("GSEA_results/BT228_down_df_LFCandQval_filt.RDS")

GSA_full_down_plot <- GSA_full_down %>%
  group_by(DRC) %>%
  arrange(q_value) %>%
  mutate(rank = row_number()) %>%
  filter(rank <= 3) %>%
  filter(q_value <= 0.05)

ggplot(GSA_full_down_plot %>% arrange(gsc_set) %>% 
         filter(gene_set %in% c("EGFR_druginhibition_82_GSE27638",
                                "MAP2K1_druginhibition_172_GSE39984",
                                "IL15_UP.V1_UP",
                                "KRAS.600_UP.V1_UP",
                                "CDK8_knockdown_63_GSE19199",
                                "EGFR_UP.V1_UP")) %>%
         mutate(gene_set = factor(gene_set, levels = gene_set)), 
       aes(x = as.factor(DRC), y = gene_set)) +
  geom_point(aes(color = -log10(q_value), size = -log10(q_value)), stroke = 0.1) +
  scale_color_gradient2(low = "white", high = "darkblue", midpoint = 0, 
                        name = "-Log10 FDR",
                        guide = "legend") +
  scale_size_continuous(range = c(1,3),
                        breaks = c(2,3,4,5,6),
                        labels = c("2", "3", "4", "5", "6+")) +
  xlab("Drug Response Class") +
  guides(size = guide_legend(title = "-Log10 FDR", nrow = 1, 
                             override.aes = list(shape = 21, stroke = 0.8)),
         color = "none") +
  theme(text = element_text(size = 6),
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.2,"cm"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.box.margin = unit(c(-3,50,-5,0), "mm"),
        # axis.text.y = element_text(angle = 45, hjust = 1),
        axis.title.y = element_blank(),
        axis.ticks = element_line(linewidth = unit(0.2, "mm"))) +
  monocle3:::monocle_theme_opts()
ggsave("GSEA_results/figures/BT228_DRC_down_GSEA_top2_FDR005_chosen_display_legend.png",
       dpi = 600, height = 1, width = 3)

ggplot(GSA_full_down_plot %>% arrange(gsc_set) %>%
         mutate(gene_set = factor(gene_set, levels = gene_set)), 
       aes(x = as.factor(DRC), y = gene_set)) +
  geom_point(aes(color = -log10(q_value), size = -log10(q_value)), stroke = 0.1) +
  scale_color_gradient2(low = "white", high = "darkblue", midpoint = 0, 
                        name = "-Log10 FDR",
                        guide = "legend") +
  scale_size_continuous(range = c(1,3)) +
  xlab("Response Module") +
  guides(size = guide_legend(title = "-Log10 FDR", nrow = 1)) +
  theme(text = element_text(size = 6),
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.2,"cm"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.box.margin = unit(c(-3,50,-5,0), "mm"),
        # axis.text.y = element_text(angle = 45, hjust = 1),
        axis.title.y = element_blank(),
        axis.ticks = element_line(linewidth = unit(0.2, "mm"))) +
  monocle3:::monocle_theme_opts()
ggsave("GSEA_results/figures/supplemental/BT228_DRC_down_GSEA_top2_FDR005_chosen_display_legend.png",
       dpi = 600, height = 3, width = 3)

# ===============================================================================
# BT112
# ===============================================================================

GSA_full_up <- data.frame()
suppressMessages(
  for(line in c("bt112")) {
    for (DRC_temp in 1:23) {
      top_genes_temp <- read_csv(paste0("DRC_up_signature_LFCandQVal_Filtered/",line,"_DRC", DRC_temp,"_up_signature_02062024.csv")) 
      
      if (length(top_genes_temp) > 1) {
        top_genes <- top_genes_temp %>%
          select(2) %>%
          pull()
        
        for(set in c("hallmarksGSC", "oncogenicGSC", "kinaseupGSC", "kinasedownGSC")){
          test <- piano::runGSAhyper(genes = top_genes, 
                                     gsc = gsc_list[[set]], 
                                     adjMethod = "none", 
                                     universe = gene_universe)
          GSAhyper_df <- as.data.frame(test$p.adj) %>%
            rownames_to_column(var = "gene_set")
          colnames(GSAhyper_df) <- c("gene_set","p_value")
          GSAhyper_df <- GSAhyper_df %>%
            mutate(cell_line = toupper(line),
                   DRC = DRC_temp,
                   gsc_set = set,
                   size = length(top_genes))
          GSA_full_up <- bind_rows(GSA_full_up, GSAhyper_df)
        }
      }
    }
  }
)
GSA_full_up$q_value <- p.adjust(GSA_full_up$p_value, method = "BH")
saveRDS(GSA_full_up, "GSEA_results/BT112_up_df_LFCandQval_filt.RDS")

GSA_full_up<- readRDS("GSEA_results/BT112_up_df_LFCandQval_filt.RDS")

GSA_full_up_plot <- GSA_full_up %>%
  group_by(DRC) %>%
  arrange(q_value) %>%
  mutate(rank = row_number()) %>%
  filter(rank <= 1)
  # filter(q_value <= 0.05)

ggplot(GSA_full_up_plot %>% arrange(gsc_set) %>%
         mutate(gene_set = factor(gene_set, levels = gene_set)), 
       aes(x = as.factor(DRC), y = gene_set)) +
  geom_point(aes(color = -log10(q_value), size = -log10(q_value)), stroke = 0.1) +
  scale_color_gradient2(low = "white", high = "red", midpoint = 0, 
                        name = "-Log10 FDR",
                        guide = "legend") +
  scale_size_continuous(range = c(1,3)) +
  xlab("Response Module") +
  guides(size = guide_legend(title = "-Log10 FDR", nrow = 1)) +
  theme(text = element_text(size = 6),
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.2,"cm"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.box.margin = unit(c(-3,50,-3,0), "mm"),
        # axis.text.y = element_text(angle = 45, hjust = 1),
        axis.title.y = element_blank(),
        axis.ticks = element_line(linewidth = unit(0.2, "mm"))) +
  monocle3:::monocle_theme_opts()

ggsave("GSEA_results/figures/supplemental/BT112_DRC_up_GSEA_top3_chosen_display_long.png",
       dpi = 600, height = 3.5, width = 3.1)


GSA_full_down <- data.frame()
suppressMessages(
  for(line in c("bt112")) {
    for (DRC_temp in 1:23) {
      top_genes_temp <- read_csv(paste0("DRC_down_signature_LFCandQVal_Filtered/",line,"_DRC", DRC_temp,"_down_signature_02062024.csv")) 
      
      if (length(top_genes_temp) > 1) {
        top_genes <- top_genes_temp %>%
          select(2) %>%
          pull()
        
        for(set in c("hallmarksGSC", "oncogenicGSC", "kinaseupGSC", "kinasedownGSC")){
          test <- piano::runGSAhyper(genes = top_genes, 
                                     gsc = gsc_list[[set]], 
                                     adjMethod = "none", 
                                     universe = gene_universe)
          GSAhyper_df <- as.data.frame(test$p.adj) %>%
            rownames_to_column(var = "gene_set")
          colnames(GSAhyper_df) <- c("gene_set","p_value")
          GSAhyper_df <- GSAhyper_df %>%
            mutate(cell_line = toupper(line),
                   DRC = DRC_temp,
                   gsc_set = set,
                   size = length(top_genes))
          GSA_full_down <- bind_rows(GSA_full_down, GSAhyper_df)
        }
      }
    }
  }
)
GSA_full_down$q_value <- p.adjust(GSA_full_down$p_value, method = "BH")
saveRDS(GSA_full_down, "GSEA_results/BT112_down_df_LFCandQval_filt.RDS")

GSA_full_down <- readRDS("GSEA_results/BT112_down_df_LFCandQval_filt.RDS")

GSA_full_down_plot <- GSA_full_down %>%
  group_by(DRC) %>%
  arrange(q_value) %>%
  mutate(rank = row_number()) %>%
  filter(rank <= 1)
  # filter(q_value <= 0.05)

ggplot(GSA_full_down_plot %>% arrange(gsc_set) %>%
         mutate(gene_set = factor(gene_set, levels = gene_set)), 
       aes(x = as.factor(DRC), y = gene_set)) +
  geom_point(aes(color = -log10(q_value), size = -log10(q_value)), stroke = 0.1) +
  scale_color_gradient2(low = "white", high = "darkblue", midpoint = 0, 
                        name = "-Log10 FDR",
                        guide = "legend") +
  scale_size_continuous(range = c(1,3)) +
  xlab("Response Module") +
  guides(size = guide_legend(title = "-Log10 FDR", nrow = 1)) +
  theme(text = element_text(size = 6),
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.2,"cm"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.box.margin = unit(c(-3,50,-3,0), "mm"),
        # axis.text.y = element_text(angle = 30, hjust = 1),
        axis.title.y = element_blank(),
        axis.ticks = element_line(linewidth = unit(0.2, "mm"))) +
  monocle3:::monocle_theme_opts()
ggsave("GSEA_results/figures/supplemental/BT112_DRC_down_GSEA_top2_FDR005_chosen_display_legend.png",
       dpi = 600, height = 3.5, width = 3.1)


# ===============================================================================
# BT333
# ===============================================================================

GSA_full_up <- data.frame()
suppressMessages(
  for(line in c("bt333")) {
    for (DRC_temp in 1:25) {
      top_genes_temp <- read_csv(paste0("DRC_up_signature_LFCandQVal_Filtered/",line,"_DRC", DRC_temp,"_up_signature_02062024.csv")) 
      
      if (length(top_genes_temp) > 1) {
        top_genes <- top_genes_temp %>%
          select(2) %>%
          pull()
        
        for(set in c("hallmarksGSC", "oncogenicGSC", "kinaseupGSC", "kinasedownGSC")){
          test <- piano::runGSAhyper(genes = top_genes, 
                                     gsc = gsc_list[[set]], 
                                     adjMethod = "none", 
                                     universe = gene_universe)
          GSAhyper_df <- as.data.frame(test$p.adj) %>%
            rownames_to_column(var = "gene_set")
          colnames(GSAhyper_df) <- c("gene_set","p_value")
          GSAhyper_df <- GSAhyper_df %>%
            mutate(cell_line = toupper(line),
                   DRC = DRC_temp,
                   gsc_set = set,
                   size = length(top_genes))
          GSA_full_up <- bind_rows(GSA_full_up, GSAhyper_df)
        }
      }
    }
  }
)
GSA_full_up$q_value <- p.adjust(GSA_full_up$p_value, method = "BH")
saveRDS(GSA_full_up, "GSEA_results/BT333_up_df_LFCandQval_filt.RDS")

GSA_full_up <- readRDS("GSEA_results/BT333_up_df_LFCandQval_filt.RDS")

GSA_full_up_plot <- GSA_full_up %>%
  group_by(DRC) %>%
  arrange(q_value) %>%
  mutate(rank = row_number()) %>%
  filter(rank <= 3) %>%
  filter(q_value <= 0.3)

ggplot(GSA_full_up_plot %>% arrange(gsc_set) %>%
         mutate(gene_set = factor(gene_set, levels = gene_set)), 
       aes(x = as.factor(DRC), y = gene_set)) +
  geom_point(aes(color = -log10(q_value), size = -log10(q_value)), stroke = 0.1) +
  scale_color_gradient2(low = "white", high = "lightpink", midpoint = 0,  
                        name = "-Log10 FDR",
                        guide = "legend") +
  scale_size_continuous(range = c(1,3)) +
  xlab("Response Module") +
  guides(size = guide_legend(title = "-Log10 FDR", nrow = 1)) +
  theme(text = element_text(size = 6),
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.2,"cm"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.box.margin = unit(c(-3,50,-3,0), "mm"),
        # axis.text.y = element_text(angle = 45, hjust = 1),
        axis.title.y = element_blank(),
        axis.ticks = element_line(linewidth = unit(0.2, "mm"))) +
  monocle3:::monocle_theme_opts()

ggsave("GSEA_results/figures/supplemental/BT333_DRC_up_GSEA_top3_chosen_display_long.png",
       dpi = 600, height = 2.25, width = 3)


GSA_full_down <- data.frame()
suppressMessages(
  for(line in c("bt333")) {
    for (DRC_temp in 1:25) {
      top_genes_temp <- read_csv(paste0("DRC_down_signature_LFCandQVal_Filtered/",line,"_DRC", DRC_temp,"_down_signature_02062024.csv")) 
      
      if (length(top_genes_temp) > 1) {
        top_genes <- top_genes_temp %>%
          select(2) %>%
          pull()
        
        for(set in c("hallmarksGSC", "oncogenicGSC", "kinaseupGSC", "kinasedownGSC")){
          test <- piano::runGSAhyper(genes = top_genes, 
                                     gsc = gsc_list[[set]], 
                                     adjMethod = "none", 
                                     universe = gene_universe)
          GSAhyper_df <- as.data.frame(test$p.adj) %>%
            rownames_to_column(var = "gene_set")
          colnames(GSAhyper_df) <- c("gene_set","p_value")
          GSAhyper_df <- GSAhyper_df %>%
            mutate(cell_line = toupper(line),
                   DRC = DRC_temp,
                   gsc_set = set,
                   size = length(top_genes))
          GSA_full_down <- bind_rows(GSA_full_down, GSAhyper_df)
        }
      }
    }
  }
)
GSA_full_down$q_value <- p.adjust(GSA_full_down$p_value, method = "BH")
saveRDS(GSA_full_down, "GSEA_results/BT333_down_df_LFCandQval_filt.RDS")

GSA_full_down <- readRDS("GSEA_results/BT333_down_df_LFCandQval_filt.RDS")

GSA_full_down_plot <- GSA_full_down %>%
  group_by(DRC) %>%
  arrange(q_value) %>%
  mutate(rank = row_number()) %>%
  filter(rank <= 3) %>%
  filter(q_value <= 0.2)

ggplot(GSA_full_down_plot %>% arrange(gsc_set) %>%
         mutate(gene_set = factor(gene_set, levels = gene_set)), 
       aes(x = as.factor(DRC), y = gene_set)) +
  geom_point(aes(color = -log10(q_value), size = -log10(q_value)), stroke = 0.1) +
  scale_color_gradient2(low = "white", high = "dodgerblue", midpoint = 0, 
                        name = "-Log10 FDR",
                        guide = "legend") +
  scale_size_continuous(range = c(1,3)) +
  xlab("Response Module") +
  guides(size = guide_legend(title = "-Log10 FDR", nrow = 1)) +
  theme(text = element_text(size = 6),
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.2,"cm"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.box.margin = unit(c(-3,50,-3,0), "mm"),
        # axis.text.y = element_text(angle = 45, hjust = 1),
        axis.title.y = element_blank(),
        axis.ticks = element_line(linewidth = unit(0.2, "mm"))) +
  monocle3:::monocle_theme_opts()
ggsave("GSEA_results/figures/supplemental/BT333_DRC_down_GSEA_top2_FDR005_chosen_display_legend.png",
       dpi = 600, height = 2, width = 3)

