library(devtools)
library(parallel)
library(ggplot2)
library(ggridges)
library(gridExtra)
library(plyr)
library(dplyr)
library(tidyr)
library(pheatmap)
library(piano)
library(UpSetR)
library(monocle3)

# Read in sci-plex repo functions
source("~/Documents/github_repos/sci-plex/bin/cell_cycle.R")
cc.genes = readRDS("~/Documents/github_repos/sci-plex/bin/cc.genes.RDS")

### Re-work data frames to join Adaptation, proliferation and MKI67 data
Adaptive_signature_df <- readRDS("adaptive_resistance_upgenes_zscored.RDS")

Adaptive_signature_df_BT112 <- Adaptive_signature_df[c("DMSO","1","2","3","4"),grepl("BT112",colnames(Adaptive_signature_df))]
colnames(Adaptive_signature_df_BT112) <- sapply(colnames(Adaptive_signature_df_BT112), function(x){stringr::str_split(x,pattern = "BT112_")[[1]][2]})
row.names(Adaptive_signature_df_BT112) <- paste0(row.names(Adaptive_signature_df_BT112),"_BT112")

Adaptive_signature_df_BT228 <- Adaptive_signature_df[c("DMSO","1","2","3","4"),grepl("BT228",colnames(Adaptive_signature_df))]
colnames(Adaptive_signature_df_BT228) <- sapply(colnames(Adaptive_signature_df_BT228), function(x){stringr::str_split(x,pattern = "BT228_")[[1]][2]})
row.names(Adaptive_signature_df_BT228) <- paste0(row.names(Adaptive_signature_df_BT228),"_BT228")

Adaptive_signature_df_BT333 <- Adaptive_signature_df[c("DMSO","1","2","3","4"),grepl("BT333",colnames(Adaptive_signature_df))]
colnames(Adaptive_signature_df_BT333) <- sapply(colnames(Adaptive_signature_df_BT333), function(x){stringr::str_split(x,pattern = "BT333_")[[1]][2]})
row.names(Adaptive_signature_df_BT333) <- paste0(row.names(Adaptive_signature_df_BT333),"_BT333")

Adaptive_signature_df_joint <- rbind(Adaptive_signature_df_BT112,Adaptive_signature_df_BT228,Adaptive_signature_df_BT333)

hmcols <- colorRampPalette(c("blue","white","red"))(35)

paletteLength <- 35

myBreaks <- c(seq(min(Adaptive_signature_df_joint[!is.na(Adaptive_signature_df_joint)]), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(Adaptive_signature_df_joint[!is.na(Adaptive_signature_df_joint)])/paletteLength, max(Adaptive_signature_df_joint[!is.na(Adaptive_signature_df_joint)]), length.out=floor(paletteLength/2)))

drug_annotations <- read.csv("Annotations_final_RG.csv")
drug_annotations <- replace(drug_annotations, drug_annotations=='', "Unknown")
drug_annotations$Reversibility <- sapply(drug_annotations$reversible, function(x){
  if(x == "Yes")return("Reversible")
  if(x == "No")return("Covalent")
  return("Unknown")
})

ann_col <- list(Reversibility = c("Covalent" = "#FF0022","Reversible" = "#011627","Unknown" = "white"),
                Class = c("acetamide" = RColorBrewer::brewer.pal(9,"Set1")[1], 
                          "aminobenzimidazole" = RColorBrewer::brewer.pal(9,"Set1")[2],
                          "aminoethylamide" = RColorBrewer::brewer.pal(9,"Set1")[3], 
                          "aminonucleoside" = RColorBrewer::brewer.pal(9,"Set1")[4],      
                          "aminopyrazine"  = RColorBrewer::brewer.pal(9,"Set1")[5], 
                          "benzoimidazole" = RColorBrewer::brewer.pal(9,"Set1")[6], 
                          "dicarboxylic acid" = RColorBrewer::brewer.pal(9,"Set1")[7], 
                          "EGFR-activator"  = RColorBrewer::brewer.pal(9,"Set1")[8], 
                          "indole" = RColorBrewer::brewer.pal(9,"Set1")[9],
                          "monoclonal antibody" = RColorBrewer::brewer.pal(8,"Set2")[1], 
                          "naturally-derived" = RColorBrewer::brewer.pal(8,"Set2")[2],
                          "PROTAC" = RColorBrewer::brewer.pal(8,"Set2")[3],             
                          "pyridine" = RColorBrewer::brewer.pal(8,"Set2")[4],           
                          "pyrimidine" = RColorBrewer::brewer.pal(8,"Set2")[5],          
                          "quinazolinamine" = RColorBrewer::brewer.pal(8,"Set2")[6], 
                          "quinazoline" = RColorBrewer::brewer.pal(8,"Set2")[7],
                          "quinoline" = RColorBrewer::brewer.pal(8,"Set2")[8], 
                          "sulfanomide" = RColorBrewer::brewer.pal(3,"Paired")[1],
                          "tyrphostin" = RColorBrewer::brewer.pal(3,"Paired")[2]))

pheatmap::pheatmap(Adaptive_signature_df_joint,
                   clustering_method = "ward.D2",
                   cluster_rows = FALSE,
                   gaps_row = c(5,10),
                   treeheight_col = 10,
                   color = hmcols,
                   breaks = myBreaks,
                   border_color = NA,
                   fontsize = 12,
                   cutree_cols = 4,
                   annotation_col = data.frame(row.names = drug_annotations$drug,
                                               Reversibility = drug_annotations$Reversibility,
                                               Class = drug_annotations$class_broad),
                   annotation_colors = ann_col,
                   file = "Heatmaps/Adaptive_signature_heatmap.png",
                   width = 14,
                   height = 6)

cds <- readRDS("cds_large_screen.RDS")

colData(cds)$MKI67_expression <- t(t(exprs(cds[rowData(cds)$gene_short_name == "MKI67",]))/colData(cds)$Size_Factor) %>% as.numeric()

MKI67_expression_df <- colData(cds) %>%
  as.data.frame() %>%
  filter(!(drug %in% c("PBS","Panitumumab","Mobocertinib"))) %>%
  group_by(cell_line, drug, dose, replicate) %>%
  summarize(MKI67_expression = log(mean(MKI67_expression) + 1)) %>%
  group_by(cell_line, drug, dose) %>%
  summarize(MKI67_expression = mean(MKI67_expression)) %>%
  group_by(cell_line) %>%
  mutate(MKI67_expression = MKI67_expression/MKI67_expression[drug == "DMSO"]) %>%
  mutate(condition_id = paste0(cell_line,"_",dose)) %>%
  ungroup() %>%
  dplyr::select(condition_id, drug, MKI67_expression) %>%
  tidyr::spread(key = drug, value = MKI67_expression) %>%
  as.data.frame()

BT112_MKI67_expression_DMSO_df <- rep(MKI67_expression_df %>% filter(condition_id == "BT112_0") %>% dplyr::select("DMSO") %>% as.numeric(), 72)
BT112_MKI67_expression_Media_df <- rep(MKI67_expression_df %>% filter(condition_id == "BT112_0") %>% dplyr::select("Media") %>% as.numeric(), 72)
BT228_MKI67_expression_DMSO_df <- rep(MKI67_expression_df %>% filter(condition_id == "BT228_0") %>% dplyr::select("DMSO") %>% as.numeric(), 72)
BT228_MKI67_expression_Media_df <- rep(MKI67_expression_df %>% filter(condition_id == "BT228_0") %>% dplyr::select("Media") %>% as.numeric(), 72)
BT333_MKI67_expression_DMSO_df <- rep(MKI67_expression_df %>% filter(condition_id == "BT333_0") %>% dplyr::select("DMSO") %>% as.numeric(), 72)
BT333_MKI67_expression_Media_df <- rep(MKI67_expression_df %>% filter(condition_id == "BT333_0") %>% dplyr::select("Media") %>% as.numeric(), 72)

MKI67_expression_vehicle_df <- rbind(BT112_MKI67_expression_DMSO_df,BT112_MKI67_expression_Media_df,
                                     BT228_MKI67_expression_DMSO_df,BT228_MKI67_expression_Media_df,
                                     BT333_MKI67_expression_DMSO_df,BT333_MKI67_expression_Media_df) %>%
  as.data.frame()

MKI67_expression_vehicle_df$condition_id <- c("DMSO_BT112","Media_BT112","DMSO_BT228","Media_BT228","DMSO_BT333","Media_BT333")
colnames_vehicle_df <- colnames(MKI67_expression_vehicle_df)
colnames_vehicle_df <- colnames_vehicle_df[colnames_vehicle_df != "condition_id"]

MKI67_expression_vehicle_df <- MKI67_expression_vehicle_df[,c("condition_id",colnames_vehicle_df)]

MKI67_expression_heatmap_df <- MKI67_expression_df
colnames_MKI67_expression_df <- colnames(MKI67_expression_heatmap_df)
colnames_MKI67_expression_df <- colnames_MKI67_expression_df[!(colnames_MKI67_expression_df %in% c("DMSO","Media"))]

colnames(MKI67_expression_vehicle_df) <- colnames_MKI67_expression_df
MKI67_expression_heatmap_df <- MKI67_expression_heatmap_df[,colnames_MKI67_expression_df]
MKI67_expression_heatmap_df <- rbind(MKI67_expression_vehicle_df,MKI67_expression_heatmap_df)

MKI67_expression_heatmap_df <- MKI67_expression_heatmap_df %>%
  filter(!(condition_id %in% c("BT112_0","BT228_0","BT333_0")))

row.names(MKI67_expression_heatmap_df) <- MKI67_expression_heatmap_df$condition_id
MKI67_expression_heatmap_df$condition_id <- NULL

row_order <- c("DMSO_BT112","Media_BT112","BT112_10","BT112_100","BT112_1000","BT112_10000",
               "DMSO_BT228","Media_BT228","BT228_10","BT228_100","BT228_1000","BT228_10000",
               "DMSO_BT333","Media_BT333","BT333_10","BT333_100","BT333_1000","BT333_10000")

MKI67_expression_ph <- pheatmap::pheatmap(MKI67_expression_heatmap_df[row_order,],
                                          clustering_method = "ward.D2",
                                          cluster_rows = FALSE,
                                          cutree_cols = 6,
                                          gaps_row = c(6,12),
                                          treeheight_col = 10,
                                          color = viridis::magma(35),
                                          fontsize = 12,
                                          border_color = NA,
                                          height = 5,
                                          width = 13,
                                          file = "Heatmaps/Mean_MIK67_Joint.png")

###
cds.list <- list()
expressed_genes.list <- list()

for(PDCL in unique(colData(cds)$cell_line) %>% sort()){
  
  cds.list[[PDCL]] <- cds[,colData(cds)$cell_line == PDCL]
  cds.list[[PDCL]] <- detect_genes(cds.list[[PDCL]])
  cds.list[[PDCL]] <- estimate_size_factors(cds.list[[PDCL]])
  
  expressed_genes.list[[PDCL]] <- rowData(cds.list[[PDCL]])[rowData(cds.list[[PDCL]])$num_cells_expressed >= nrow(cds.list[[PDCL]])*0.01,]$id
  print(length(expressed_genes.list[[PDCL]]))
  
}

# Re-do MKI67 heatmap with DEG information
treatment_diff_test.list <- readRDS("full_treatment_diff_test.list.rds")

for(PDCL in sort(unique(colData(cds)$cell_line))){
  
  treatment_diff_test.list[[PDCL]] <- do.call("rbind",treatment_diff_test.list[[PDCL]])
  treatment_diff_test.list[[PDCL]]$cell_line <- rep(PDCL, nrow(treatment_diff_test.list[[PDCL]]))
  treatment_diff_test.list[[PDCL]] <- treatment_diff_test.list[[PDCL]] %>% filter(id %in% expressed_genes.list[[PDCL]])
  treatment_diff_test.list[[PDCL]]$q_value <- p.adjust(treatment_diff_test.list[[PDCL]]$p_value, method = "BH")
  
}

treatment_diff_test_results <- do.call("rbind",treatment_diff_test.list)

BT112_degs_df <- treatment_diff_test_results %>% 
  filter(term == "log(dose + 0.1)", q_value < 0.01, abs(normalized_effect) > 0.05, cell_line == "BT112", drug != "Vehicle") %>%
  group_by(cell_line,drug) %>%
  summarize(degs = n())

BT112_max <- BT112_degs_df %>% filter(drug != "CUDC-101") %>% pull(degs) %>% max()
BT112_degs_df <- BT112_degs_df %>% mutate(degs = ifelse(degs > BT112_max, BT112_max, degs))

BT228_degs_df <- treatment_diff_test_results %>% 
  filter(term == "log(dose + 0.1)", q_value < 0.01, abs(normalized_effect) > 0.05, cell_line == "BT228", drug != "Vehicle") %>%
  group_by(cell_line,drug) %>%
  summarize(degs = n())

BT228_max <- BT228_degs_df %>% filter(drug != "CUDC-101") %>% pull(degs) %>% max()
BT228_degs_df <- BT228_degs_df %>% mutate(degs = ifelse(degs > BT228_max, BT228_max, degs))

BT333_degs_df <- treatment_diff_test_results %>% 
  filter(term == "log(dose + 0.1)", q_value < 0.01, abs(normalized_effect) > 0.05, cell_line == "BT333", drug != "Vehicle") %>%
  group_by(cell_line,drug) %>%
  summarize(degs = n()) 

BT333_max <- BT333_degs_df %>% filter(drug != "CUDC-101") %>% pull(degs) %>% max()
BT333_degs_df <- BT333_degs_df %>% mutate(degs = ifelse(degs > BT333_max, BT333_max, degs))

pheatmap::pheatmap(MKI67_expression_heatmap_df[row_order[grepl("BT112",row_order)],],
                   clustering_method = "ward.D2",
                   cluster_rows = FALSE,
                   cutree_cols = 6,
                   cluster_cols = MKI67_expression_ph$tree_col,
                   treeheight_col = 10,
                   color = viridis::plasma(35),
                   annotation_col = data.frame(row.names = BT112_degs_df$drug, log10_DEGs = log10(BT112_degs_df$degs)),
                   annotation_colors = list(log10_DEGs = viridis::magma(35)),
                   fontsize = 12,
                   border_color = NA,
                   height = 3.5,
                   width = 14,
                   file = "Heatmaps/Mean_MIK67_BT112.png")

pheatmap::pheatmap(MKI67_expression_heatmap_df[row_order[grepl("BT228",row_order)],],
                   clustering_method = "ward.D2",
                   cluster_rows = FALSE,
                   cutree_cols = 6,
                   cluster_cols = MKI67_expression_ph$tree_col,
                   treeheight_col = 10,
                   color = viridis::plasma(35),
                   annotation_col = data.frame(row.names = BT228_degs_df$drug, log10_DEGs = log10(BT228_degs_df$degs)),
                   annotation_colors = list(log10_DEGs = viridis::magma(35)),
                   fontsize = 12,
                   border_color = NA,
                   height = 3.5,
                   width = 14,
                   file = "Heatmaps/Mean_MIK67_BT228.png")

pheatmap::pheatmap(MKI67_expression_heatmap_df[row_order[grepl("BT333",row_order)],],
                   clustering_method = "ward.D2",
                   cluster_rows = FALSE,
                   cutree_cols = 6,
                   cluster_cols = MKI67_expression_ph$tree_col,
                   treeheight_col = 10,
                   color = viridis::plasma(35),
                   annotation_col = data.frame(row.names = BT333_degs_df$drug, log10_DEGs = log10(BT333_degs_df$degs)),
                   annotation_colors = list(log10_DEGs = viridis::magma(35)),
                   fontsize = 12,
                   border_color = NA,
                   height = 3.5,
                   width = 14,
                   file = "Heatmaps/Mean_MIK67_BT333.png")

### Merge adaptation and MKI67 data frames
Adaptive_signature_df_long <- reshape2::melt(Adaptive_signature_df)
colnames(Adaptive_signature_df_long) <- c("dose", "condition_id", "adaptive_score")

Adaptive_signature_df_long$new_dose <- sapply(Adaptive_signature_df_long$dose, function(x){
  if(x == "DMSO")return("DMSO")
  if(x == "PBS")return("PBS")
  if(x == "Media")return("Media")
  if(x == "1") return("10")
  if(x == "2") return("100")
  if(x == "3") return("1000")
  if(x == "4") return("10000")
  return(NA)
})

Adaptive_signature_df_long$new_drug <- sapply(Adaptive_signature_df_long$condition_id, function(x){
  stringr::str_split(x, pattern = "_")[[1]][2]
})

Adaptive_signature_df_long$new_PDCL <- sapply(Adaptive_signature_df_long$condition_id, function(x){
  stringr::str_split(x, pattern = "_")[[1]][1]
})

Adaptive_signature_df_long$new_condition_id <- paste0(Adaptive_signature_df_long$new_drug,"_",
                                                      Adaptive_signature_df_long$new_dose,"_",
                                                      Adaptive_signature_df_long$new_PDCL)

Adaptive_signature_df_long <- Adaptive_signature_df_long %>% filter(new_dose %in% c("10","100","1000","10000"))

MKI67_expression_heatmap_df_long <- reshape2::melt(t(MKI67_expression_heatmap_df))
colnames(MKI67_expression_heatmap_df_long) <- c("drug", "condition_id", "MKI67_score")


MKI67_expression_heatmap_df_long$new_dose <- sapply(MKI67_expression_heatmap_df_long$condition_id, function(x){
  stringr::str_split(x, pattern = "_")[[1]][2]
})
MKI67_expression_heatmap_df_long <- MKI67_expression_heatmap_df_long %>% filter(new_dose %in% c("10","100","1000","10000"))
MKI67_expression_heatmap_df_long$new_PDCL <- sapply(MKI67_expression_heatmap_df_long$condition_id, function(x){
  stringr::str_split(x, pattern = "_")[[1]][1]
})
MKI67_expression_heatmap_df_long$new_drug <- MKI67_expression_heatmap_df_long$drug
MKI67_expression_heatmap_df_long$new_condition_id <- paste0(MKI67_expression_heatmap_df_long$new_drug,"_",
                                                            MKI67_expression_heatmap_df_long$new_dose,"_",
                                                            MKI67_expression_heatmap_df_long$new_PDCL)

joint_MKI67_Adaptive_df <- left_join(MKI67_expression_heatmap_df_long %>% dplyr::select(new_condition_id, MKI67_score), 
                                     Adaptive_signature_df_long %>% dplyr::select(new_condition_id, adaptive_score),
                                     by = "new_condition_id")

joint_MKI67_Adaptive_df$PDCL <- sapply(joint_MKI67_Adaptive_df$new_condition_id, function(x){
  stringr::str_split(x, pattern = "_")[[1]][3]})
joint_MKI67_Adaptive_df$drug <- sapply(joint_MKI67_Adaptive_df$new_condition_id, function(x){
  stringr::str_split(x, pattern = "_")[[1]][1]})
joint_MKI67_Adaptive_df$dose <- sapply(joint_MKI67_Adaptive_df$new_condition_id, function(x){
  stringr::str_split(x, pattern = "_")[[1]][2]})

joint_MKI67_Adaptive_df <- joint_MKI67_Adaptive_df %>% filter(!is.na(MKI67_score), !is.na(adaptive_score))

mean_adaptive_score <- mean(joint_MKI67_Adaptive_df$adaptive_score)
sd_adaptive_score <- sd(joint_MKI67_Adaptive_df$adaptive_score)
mean_MKI67_score <- mean(joint_MKI67_Adaptive_df$MKI67_score)
sd_MKI67_score <- sd(joint_MKI67_Adaptive_df$MKI67_score)

joint_MKI67_Adaptive_df <- joint_MKI67_Adaptive_df %>%
  mutate(phenotype = ifelse(MKI67_score < mean_MKI67_score - sd_MKI67_score & adaptive_score < mean_adaptive_score + sd_adaptive_score, "Viability low / Adaptive low",
                            ifelse(MKI67_score < mean_MKI67_score - sd_MKI67_score & adaptive_score > mean_adaptive_score + (2*sd_adaptive_score), "Viability low / Adaptive high",
                                   "Other")))

joint_MKI67_Adaptive_df$phenotype <- factor(joint_MKI67_Adaptive_df$phenotype, levels = c("Viability low / Adaptive low","Viability low / Adaptive high","Other"))

ggplot(joint_MKI67_Adaptive_df, aes(x = MKI67_score, y = adaptive_score, 
                                    fill = phenotype, color = log10(as.numeric(dose))),) +
  geom_point(size = 1, stroke = 0.3, shape = 21) +
  facet_wrap(~PDCL, ncol = 3) +
  geom_hline(yintercept = mean_adaptive_score + sd_adaptive_score, size = 0.25) +
  geom_hline(yintercept = mean_adaptive_score + (2*sd_adaptive_score), linetype = "dashed", size = 0.25) +
  geom_vline(xintercept = mean_MKI67_score - sd_MKI67_score, size = 0.25) +
  ggrepel::geom_text_repel(data = joint_MKI67_Adaptive_df %>% filter(dose %in% c("1000","10000"), phenotype %in% c("Viability low / Adaptive low","Viability low / Adaptive high")), 
                           aes(x=MKI67_score,y=adaptive_score,label=drug),
                           size = 1.75,
                           force = 5,
                           segment.size = 0.1,
                           #label.size = 0.1,
                           box.padding = 0.1,
                           #label.padding = 0.1,
                           max.overlaps = 100,
                           color = "black") +
  theme(text = element_text(size = 9),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.position = "none",
        legend.key.width = unit(1,"line"), 
        legend.key.height = unit(0.5,"line")) +
  monocle3:::monocle_theme_opts() +
  scale_fill_manual("Class", values = c("Viability low / Adaptive low" = "darkcyan", "Viability low / Adaptive high" = "firebrick2", "Other" = "grey90")) +
  viridis::scale_color_viridis("log10(Dose [nM])", option = "magma") +
  guides(fill = guide_legend(override.aes = list(size=4))) +
  ylab("Adaptive score") +
  xlab("MKI67 score")
ggsave("Marker_plots/Proliferation_vs_adaptation_supplemental.png", dpi = 900, height = 3.5, width = 6.5)

ggplot(joint_MKI67_Adaptive_df, aes(x = MKI67_score, y = adaptive_score, 
                                    fill = phenotype, color = log10(as.numeric(dose))),) +
  geom_point(size = 1, stroke = 0.3, shape = 21) +
  facet_wrap(~PDCL, ncol = 3) +
  geom_hline(yintercept = mean_adaptive_score + sd_adaptive_score, size = 0.25) +
  geom_hline(yintercept = mean_adaptive_score + (2*sd_adaptive_score), linetype = "dashed", size = 0.25) +
  geom_vline(xintercept = mean_MKI67_score - sd_MKI67_score, size = 0.25) +
  ggrepel::geom_text_repel(data = joint_MKI67_Adaptive_df %>% filter(dose %in% c("1000","10000"), phenotype %in% c("Viability low / Adaptive low","Viability low / Adaptive high")), 
                  aes(x=MKI67_score,y=adaptive_score,label=drug),
                  size = 1.75,
                  force = 5,
                  segment.size = 0.1,
                  #label.size = 0.1,
                  box.padding = 0.1,
                  #label.padding = 0.1,
                  max.overlaps = 100,
                  color = "black") +
  theme(text = element_text(size = 9),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.position = "right",
        legend.key.width = unit(1,"line"), 
        legend.key.height = unit(0.5,"line")) +
  monocle3:::monocle_theme_opts() +
  scale_fill_manual("Class", values = c("Viability low / Adaptive low" = "darkcyan", "Viability low / Adaptive high" = "firebrick2", "Other" = "grey90")) +
  viridis::scale_color_viridis("log10(Dose [nM])", option = "magma") +
  guides(fill = guide_legend(override.aes = list(size=4))) +
  ylab("Adaptive score") +
  xlab("MKI67 score")
ggsave("Marker_plots/Proliferation_vs_adaptation_for_legend.png", dpi = 900, height = 4, width = 7.5)

cor.test(joint_MKI67_Adaptive_df$adaptive_score, joint_MKI67_Adaptive_df$MKI67_score)
cor.test(joint_MKI67_Adaptive_df[joint_MKI67_Adaptive_df$drug == "Osimertinib",]$adaptive_score, joint_MKI67_Adaptive_df[joint_MKI67_Adaptive_df$drug == "Osimertinib",]$MKI67_score)
cor.test(joint_MKI67_Adaptive_df[joint_MKI67_Adaptive_df$drug == "MTX-211",]$adaptive_score, joint_MKI67_Adaptive_df[joint_MKI67_Adaptive_df$drug == "MTX-211",]$MKI67_score)

