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
library(monocle3)

source("~/Documents/github_repos/sci-plex/bin/cell_cycle.R")
cc.genes = readRDS("~/Documents/github_repos/sci-plex/bin/cc.genes.RDS")

calculate_aggregate_expression_score <- function(cds, signature_genes, from_id = FALSE){
  
  if(from_id == TRUE){
    cds_subset = cds[rowData(cds)$id %in% signature_genes,]
  }
  else(cds_subset = cds[rowData(cds)$gene_short_name %in% signature_genes,])
  aggregate_signature_expression = exprs(cds_subset)
  aggregate_signature_expression = t(t(aggregate_signature_expression) / pData(cds_subset)$Size_Factor)
  aggregate_signature_expression = Matrix::colSums(aggregate_signature_expression)
  signature_score = log(aggregate_signature_expression+1)
  return(aggregate_signature_expression)
}

cds <- readRDS("cds_large_screen.RDS")
colData(cds)$cell <- colData(cds)$cell_ID

cds <- estimate_cell_cycle(cds, 
                           g1s_markers = cc.genes$s.genes, 
                           g2m_markers = cc.genes$g2m.genes)

proliferation_index_heatmap_summary <- colData(cds) %>%
  as.data.frame() %>%
  ungroup() %>%
  mutate(dose = as.character(dose)) %>%
  group_by(cell_line, drug, dose) %>%
  summarise(mean_proliferation_index = mean(proliferation_index)) %>%
  ungroup() %>%
  mutate(mean_proliferation_index = mean_proliferation_index/mean_proliferation_index[drug == "DMSO" & dose == "0"]) 

cds <- detect_genes(cds)
expressed_genes <- rowData(cds)[rowData(cds)$num_cells_expressed >= ncol(cds)*.01,]$id

length(expressed_genes)

# genesets can be downloaded online
hallmarksGSC <- loadGSCSafe(file="h.all.v6.0.symbols.gmt")
OncogenicSignaturesGSC <- loadGSCSafe(file="c6.all.v6.0.OncogenicSignatures.symbols.gmt")

colData(cds)$KRAS_up_score <- calculate_aggregate_expression_score(cds, signature_genes = hallmarksGSC$gsc$HALLMARK_KRAS_SIGNALING_UP)
colData(cds)$KRAS_dn_score <- calculate_aggregate_expression_score(cds, signature_genes = hallmarksGSC$gsc$HALLMARK_KRAS_SIGNALING_DN)
colData(cds)$PI3K_AKT_MTOR_score <- calculate_aggregate_expression_score(cds, signature_genes = hallmarksGSC$gsc$HALLMARK_PI3K_AKT_MTOR_SIGNALING)
colData(cds)$MTORC1_score <- calculate_aggregate_expression_score(cds, signature_genes = hallmarksGSC$gsc$HALLMARK_MTORC1_SIGNALING)
colData(cds)$E2F_score <- calculate_aggregate_expression_score(cds, signature_genes = hallmarksGSC$gsc$HALLMARK_E2F_TARGETS)
colData(cds)$G2M_checkpoint_score <- calculate_aggregate_expression_score(cds, signature_genes = hallmarksGSC$gsc$HALLMARK_G2M_CHECKPOINT)
colData(cds)$P53_score <- calculate_aggregate_expression_score(cds, signature_genes = hallmarksGSC$gsc$HALLMARK_P53_PATHWAY)
colData(cds)$EGFR_score <- calculate_aggregate_expression_score(cds, signature_genes = OncogenicSignaturesGSC$gsc$EGFR_UP.V1_UP)
colData(cds)$MEK_score <- calculate_aggregate_expression_score(cds, signature_genes = OncogenicSignaturesGSC$gsc$MEK_UP.V1_UP)
colData(cds)$AKT_score <- calculate_aggregate_expression_score(cds, signature_genes = OncogenicSignaturesGSC$gsc$AKT_UP.V1_UP)

APM_MHC_CI_genes <- c("PSMB5",
               "PSMB6",
               "PSMB7",
               "PSMB8",
               "PSMB9",
               "PSMB10",
               "TAP1",
               "TAP2",
               "ERAP1",
               "ERAP2",
               "CANX",
               "CALR",
               "PDIA3",
               "TAPBP",
               "B2M",
               "HLA-A",
               "HLA-B",
               "HLA-C")

APM_MHC_CII_genes <- c("CIITA", 
                       "HLA-DRA", 
                       "HLA-DRB1", 
                       "HLA-DRB5", 
                       "HLA-DRB6", 
                       "HLA-DQA1", 
                       "HLA-DQB1", 
                       "HLA-DPA1", 
                       "HLA-DMA", 
                       "HLA-DMB", 
                       "HLA-DOA")

colData(cds)$APM_MHC_CI_score <- calculate_aggregate_expression_score(cds, signature_genes = APM_MHC_CI_genes)
colData(cds)$APM_MHC_CII_score <- calculate_aggregate_expression_score(cds, signature_genes = APM_MHC_CII_genes)

drugs <- unique(colData(cds)$drug)
drugs <- drugs[!(drugs %in% c("DMSO","PBS","Media", "Panitumumab"))] %>% sort()

gene_set_score_diff_test.list <- list()

for(PDCL in c("BT112","BT228","BT333")){
  
  gene_set_score_diff_test.list[[PDCL]] <- list()
  
  for(treatment in drugs){
    
    gene_set_score_diff_test.list[[PDCL]][[treatment]] <- list()
    
    for(gene_set in c("EGFR_score","MEK_score","AKT_score","KRAS_up_score","KRAS_dn_score","PI3K_AKT_MTOR_score","MTORC1_score","E2F_score","G2M_checkpoint_score","P53_score","APM_MHC_CI_score","APM_MHC_CII_score")){
    
      score_df_subset <- colData(cds) %>% as.data.frame() %>% filter(cell_line == PDCL, drug %in% c("DMSO",treatment))
      
      model_to_test <- paste0(gene_set," ~ log(dose + 0.1) + replicate")
        
      diff_test_results <- speedglm::speedglm(data = score_df_subset, formula = model_to_test)
      diff_test_results <- broom::tidy(diff_test_results)
      
      gene_set_score_diff_test.list[[PDCL]][[treatment]][[gene_set]] <- data.frame(cell_line = PDCL,
                                                                                   signature = gene_set,
                                                                                   drug = treatment, 
                                                                                   term = diff_test_results$term, 
                                                                                   beta = diff_test_results$estimate, 
                                                                                   p_value = diff_test_results$p.value)
    message("Finished ",gene_set)
      }
    
  }
  
  message("Finished ",PDCL)
  
}

for(PDCL in c("BT112","BT228","BT333")){
  
  for(treatment in "Panitumumab"){
    
    for(gene_set in c("EGFR_score","MEK_score","AKT_score","KRAS_up_score","KRAS_dn_score","PI3K_AKT_MTOR_score","MTORC1_score","E2F_score","G2M_checkpoint_score","P53_score","APM_MHC_CI_score","APM_MHC_CII_score")){
      
      score_df_subset <- colData(cds) %>% as.data.frame() %>% filter(cell_line == PDCL, drug %in% c("PBS",treatment)) %>%
        mutate(dose = dose*10)
    
    model_to_test <- paste0(gene_set," ~ log(dose + 0.1) + replicate")
    
    diff_test_results <- speedglm::speedglm(data = score_df_subset, formula = model_to_test)
    diff_test_results <- broom::tidy(diff_test_results)
    
    gene_set_score_diff_test.list[[PDCL]][[treatment]][[gene_set]] <- data.frame(cell_line = PDCL,
                                                                                 signature = gene_set,
                                                                                 drug = treatment, 
                                                                                 term = diff_test_results$term, 
                                                                                 beta = diff_test_results$estimate, 
                                                                                 p_value = diff_test_results$p.value)
    
    message("Finished ",gene_set)
    
    }
    
  }
  
  message("Finished ",PDCL)
  
}

gene_set_score_diff_test_results.list <- gene_set_score_diff_test.list

for(PDCL in c("BT112","BT228","BT333")){
  
  for(treatment in names(gene_set_score_diff_test_results.list[[PDCL]])){
    
  gene_set_score_diff_test_results.list[[PDCL]][[treatment]] <- do.call("rbind",gene_set_score_diff_test.list[[PDCL]][[treatment]])
  
  }

}

for(PDCL in c("BT112","BT228","BT333")){
    
  gene_set_score_diff_test_results.list[[PDCL]] <- do.call("rbind",gene_set_score_diff_test_results.list[[PDCL]])
  gene_set_score_diff_test_results.list[[PDCL]]$q_value <- p.adjust(gene_set_score_diff_test_results.list[[PDCL]]$p_value, method = "BH")
  
}

gene_set_score_diff_test_results <- do.call("rbind",gene_set_score_diff_test_results.list)
# saveRDS(gene_set_score_diff_test_results, "gene_set_score_diff_test_results.rds")
# gene_set_score_diff_test_results <- readRDS("gene_set_score_diff_test_results.rds")

APM_MHC_CI_score_results <- gene_set_score_diff_test_results %>%
  filter(term == "log(dose + 0.1)", signature == "APM_MHC_CI_score") %>%
  group_by(cell_line, drug) %>%
  arrange(desc(beta))

### Focus in on APM scores. Signaling scores are difficult to interpret given adaptation and other phenotypes, can be revisited
BT112_APM_MHCI_perturbing_drugs <- gene_set_score_diff_test_results %>%
  filter(term == "log(dose + 0.1)", q_value < 0.05, cell_line == "BT112",signature %in% c("APM_MHC_CI_score")) %>%
  pull(drug) %>%
  unique() %>%
  sort()

BT228_APM_MHCI_perturbing_drugs <- gene_set_score_diff_test_results %>%
  filter(term == "log(dose + 0.1)", q_value < 0.05, cell_line == "BT228",signature %in% c("APM_MHC_CI_score")) %>%
  pull(drug) %>%
  unique() %>%
  sort()

BT333_APM_MHCI_perturbing_drugs <- gene_set_score_diff_test_results %>%
  filter(term == "log(dose + 0.1)", q_value < 0.05, cell_line == "BT333",signature %in% c("APM_MHC_CI_score")) %>%
  pull(drug) %>%
  unique() %>%
  sort()

union_APM_MHCI_perturbing_drugs <- Reduce("union", list(BT112_APM_MHCI_perturbing_drugs,BT228_APM_MHCI_perturbing_drugs,BT333_APM_MHCI_perturbing_drugs))

APM_MHCI_score_sig_results <- gene_set_score_diff_test_results %>%
  filter(term == "log(dose + 0.1)", drug %in% union_APM_MHCI_perturbing_drugs, signature == "APM_MHC_CI_score") %>%
  mutate(beta = ifelse(q_value <= 0.05, beta, 0)) %>%
  dplyr::select(cell_line, drug, beta) %>%
  tidyr::spread(key = cell_line, value = beta) %>%
  as.data.frame()

row.names(APM_MHCI_score_sig_results) <- APM_MHCI_score_sig_results$drug
APM_MHCI_score_sig_results$drug <- NULL

col_fun1 = circlize::colorRamp2(c(min(APM_MHCI_score_sig_results), 0, max(APM_MHCI_score_sig_results)), c("blue", "white", "red"))

png("Marker_plots/APM_MHCI_scores.png", res = 600, pointsize = 9, units = "in", width = 2, height = 2)
circlize::circos.heatmap(APM_MHCI_score_sig_results %>% as.matrix(), 
                         col = circlize::colorRamp2(c(min(APM_MHCI_score_sig_results), 0, max(APM_MHCI_score_sig_results)), c("blue", "white", "red")), 
                         dend.side = "inside",
                         rownames.side = "outside",
                         dend.track.height = 0.5)
circlize::circos.clear()
dev.off()

png("Marker_plots/APM_MHCI_scores_legend.png", res = 600, pointsize = 9, units = "in", width = 2, height = 2)
circlize::circos.heatmap(APM_MHCI_score_sig_results %>% as.matrix(), 
                         col = circlize::colorRamp2(c(min(APM_MHCI_score_sig_results), 0, max(APM_MHCI_score_sig_results)), c("blue", "white", "red")), 
                         dend.side = "inside",
                         rownames.side = "outside",
                         dend.track.height = 0.5)
circlize::circos.clear()
lgd_1 = ComplexHeatmap::Legend(title = "APM_CI", col_fun = col_fun1)
grid::grid.draw(lgd_1)
dev.off()

BT112_APM_MHCII_perturbing_drugs <- gene_set_score_diff_test_results %>%
  filter(term == "log(dose + 0.1)", q_value < 0.05, cell_line == "BT112",signature %in% c("APM_MHC_CII_score")) %>%
  pull(drug) %>%
  unique() %>%
  sort()

BT228_APM_MHCII_perturbing_drugs <- gene_set_score_diff_test_results %>%
  filter(term == "log(dose + 0.1)", q_value < 0.05, cell_line == "BT228",signature %in% c("APM_MHC_CII_score")) %>%
  pull(drug) %>%
  unique() %>%
  sort()

BT333_APM_MHCII_perturbing_drugs <- gene_set_score_diff_test_results %>%
  filter(term == "log(dose + 0.1)", q_value < 0.05, cell_line == "BT333",signature %in% c("APM_MHC_CII_score")) %>%
  pull(drug) %>%
  unique() %>%
  sort()

union_APM_MHCII_perturbing_drugs <- Reduce("union", list(BT112_APM_MHCII_perturbing_drugs,BT228_APM_MHCII_perturbing_drugs,BT333_APM_MHCII_perturbing_drugs))

APM_MHCII_score_sig_results <- gene_set_score_diff_test_results %>%
  filter(term == "log(dose + 0.1)", drug %in% union_APM_MHCII_perturbing_drugs, signature == "APM_MHC_CII_score") %>%
  mutate(beta = ifelse(q_value <= 0.05, beta, 0)) %>%
  dplyr::select(cell_line, drug, beta) %>%
  tidyr::spread(key = cell_line, value = beta) %>%
  as.data.frame()

row.names(APM_MHCII_score_sig_results) <- APM_MHCII_score_sig_results$drug
APM_MHCII_score_sig_results$drug <- NULL

png("Marker_plots/APM_MHCII_scores.png", res = 600, pointsize = 9, units = "in", width = 2, height = 2)
circlize::circos.heatmap(APM_MHCII_score_sig_results %>% as.matrix(), 
                         col = circlize::colorRamp2(c(min(APM_MHCII_score_sig_results), 0, max(APM_MHCII_score_sig_results)), c("blue", "white", "red")), 
                         dend.side = "inside",
                         rownames.side = "outside",
                         dend.track.height = 0.5)
circlize::circos.clear()
dev.off()

col_fun2 = circlize::colorRamp2(c(min(APM_MHCII_score_sig_results), 0, max(APM_MHCII_score_sig_results)), c("blue", "white", "red"))

png("Marker_plots/APM_MHCII_scores_legend.png", res = 600, pointsize = 9, units = "in", width = 2, height = 2)
circlize::circos.heatmap(APM_MHCII_score_sig_results %>% as.matrix(), 
                         col = circlize::colorRamp2(c(min(APM_MHCII_score_sig_results), 0, max(APM_MHCII_score_sig_results)), c("blue", "white", "red")), 
                         dend.side = "inside",
                         rownames.side = "outside",
                         dend.track.height = 0.5)
circlize::circos.clear()
lgd_2 = ComplexHeatmap::Legend(title = "APM_CII", col_fun = col_fun2)
grid::grid.draw(lgd_2)
dev.off()


#### Read in DEG test results to examine individual APM genes
cds.list <- list()
expressed_genes.list <- list()

for(PDCL in unique(colData(cds)$cell_line) %>% sort()){
  
  cds.list[[PDCL]] <- cds[,colData(cds)$cell_line == PDCL]
  cds.list[[PDCL]] <- detect_genes(cds.list[[PDCL]])
  cds.list[[PDCL]] <- estimate_size_factors(cds.list[[PDCL]])
  
  expressed_genes.list[[PDCL]] <- rowData(cds.list[[PDCL]])[rowData(cds.list[[PDCL]])$num_cells_expressed >= nrow(cds.list[[PDCL]])*0.01,]$id
  print(length(expressed_genes.list[[PDCL]]))
  
}

treatment_diff_test.list <- readRDS("full_treatment_diff_test.list.rds")

for(PDCL in sort(unique(colData(cds)$cell_line))){
  
  treatment_diff_test.list[[PDCL]] <- do.call("rbind",treatment_diff_test.list[[PDCL]])
  treatment_diff_test.list[[PDCL]]$cell_line <- rep(PDCL, nrow(treatment_diff_test.list[[PDCL]]))
  treatment_diff_test.list[[PDCL]] <- treatment_diff_test.list[[PDCL]] %>% filter(id %in% expressed_genes.list[[PDCL]])
  treatment_diff_test.list[[PDCL]]$q_value <- p.adjust(treatment_diff_test.list[[PDCL]]$p_value, method = "BH")
  
}

treatment_diff_test_results <- do.call("rbind",treatment_diff_test.list)

### Focus on APM alterations in 2 or more lines
top_APM_MHC_CI_genes <- treatment_diff_test_results %>%
  filter(term == "log(dose + 0.1)", 
         gene_short_name %in% APM_MHC_CI_genes,
         q_value < 0.05,
         drug %in% c("AG494","AG-18","AG555","AG-490","TyrphostinAG-528","Tyrphostin9","CUDC-101","Puromycin")) %>%
  pull(gene_short_name) %>%
  unique()

top_APM_MHC_CI_genes_matrix <- treatment_diff_test_results %>%
  filter(term == "log(dose + 0.1)", 
         # gene_short_name %in% top_APM_MHC_CI_genes,
         gene_short_name %in% APM_MHC_CI_genes,
         drug %in% c("AG494","AG-18","AG555","AG-490","TyrphostinAG-528","Tyrphostin9","CUDC-101","Puromycin")) %>%
  mutate(condition_id = paste0(cell_line,"_",drug)) %>%
  dplyr::select(condition_id, gene_short_name, normalized_effect) %>%
  tidyr::spread(key = condition_id, value = normalized_effect) %>%
  as.data.frame()

row.names(top_APM_MHC_CI_genes_matrix) <- top_APM_MHC_CI_genes_matrix$gene_short_name
top_APM_MHC_CI_genes_matrix$gene_short_name <- NULL

top_APM_MHC_CI_genes_matrix[is.na(top_APM_MHC_CI_genes_matrix)] <- 0

hmcols <- colorRampPalette(c("blue","white","red"))(35)

paletteLength <- 35
myBreaks <- c(seq(min(top_APM_MHC_CI_genes_matrix), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(top_APM_MHC_CI_genes_matrix)/paletteLength, max(top_APM_MHC_CI_genes_matrix), length.out=floor(paletteLength/2)))

annotation_df <- data.frame(row.names = colnames(top_APM_MHC_CI_genes_matrix), 
                            condition_id = colnames(top_APM_MHC_CI_genes_matrix))

annotation_df$cell_line <- sapply(annotation_df$condition_id, function(x){stringr::str_split(x, pattern = "_")[[1]][1]})
annotation_df$drug <- sapply(annotation_df$condition_id, function(x){stringr::str_split(x, pattern = "_")[[1]][2]})

drug_order_1 <- c("BT112_AG-18",
                  "BT228_AG-18",
                  "BT333_AG-18",
                  "BT112_AG-490",
                  "BT228_AG-490",
                  "BT333_AG-490",
                  "BT112_AG494",
                  "BT228_AG494",
                  "BT333_AG494",
                  "BT112_AG555",
                  "BT228_AG555",
                  "BT333_AG555",
                  "BT112_TyrphostinAG-528",
                  "BT228_TyrphostinAG-528",
                  "BT333_TyrphostinAG-528",
                  "BT112_Tyrphostin9",
                  "BT228_Tyrphostin9",
                  "BT333_Tyrphostin9",
                  "BT112_CUDC-101",
                  "BT228_CUDC-101",
                  "BT333_CUDC-101")

pheatmap::pheatmap(top_APM_MHC_CI_genes_matrix[,drug_order_1],
                   clustering_method = "ward.D2",
                   cluster_cols = FALSE,
                   color = hmcols,
                   breaks = myBreaks,
                   gaps_col = c(3,6,9,12,15,18,21),
                   cutree_rows = 6,
                   treeheight_row = 10,
                   fontsize = 6,
                   fontsize_col = 4,
                   height = 2.75,
                   width = 3,
                   border_color = NA,
                   file = "Marker_plots/top_APM_MHC_CI_genes.png")

### CII heatmap
top_APM_MHC_CII_genes <- treatment_diff_test_results %>%
  filter(term == "log(dose + 0.1)", 
         gene_short_name %in% APM_MHC_CII_genes,
         q_value < 0.05,
         drug %in% union_APM_MHCII_perturbing_drugs) %>%
  pull(gene_short_name) %>%
  unique()

top_APM_MHC_CII_genes_matrix <- treatment_diff_test_results %>%
  filter(term == "log(dose + 0.1)", 
         gene_short_name %in% APM_MHC_CII_genes,
         drug %in% union_APM_MHCII_perturbing_drugs,
         cell_line == "BT333") %>%
  mutate(condition_id = drug) %>%
  dplyr::select(condition_id, gene_short_name, normalized_effect) %>%
  tidyr::spread(key = condition_id, value = normalized_effect) %>%
  as.data.frame()

row.names(top_APM_MHC_CII_genes_matrix) <- top_APM_MHC_CII_genes_matrix$gene_short_name
top_APM_MHC_CII_genes_matrix$gene_short_name <- NULL

top_APM_MHC_CII_genes_matrix[is.na(top_APM_MHC_CII_genes_matrix)] <- 0

hmcols_2 <- colorRampPalette(c("blue","white","red"))(35)

paletteLength <- 35
myBreaks_2 <- c(seq(min(top_APM_MHC_CII_genes_matrix), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(top_APM_MHC_CII_genes_matrix)/paletteLength, max(top_APM_MHC_CII_genes_matrix), length.out=floor(paletteLength/2)))

annotation_df_2 <- data.frame(row.names = colnames(top_APM_MHC_CII_genes_matrix), 
                            condition_id = colnames(top_APM_MHC_CII_genes_matrix))

annotation_df_2$cell_line <- sapply(annotation_df_2$condition_id, function(x){stringr::str_split(x, pattern = "_")[[1]][1]})
annotation_df_2$drug <- sapply(annotation_df_2$condition_id, function(x){stringr::str_split(x, pattern = "_")[[1]][2]})

drug_order_CII_1 <- c("CNX-2006",
                  "MTX-211",
                  "Avitinib",
                  "AZ5104",
                  "Poziotinib",
                  "Rociletinib",
                  "CL-387785",
                  "Gefitinib",
                  "AEE788",
                  "TAS6417",
                  "Alflutinib",
                  "AZD3759",
                  "Panitumumab",
                  "AG-1478",
                  "AP26113-analog",
                  "WZ3146",
                  "AG-1557",
                  "PD153035",
                  "Sapitinib",
                  "Naquotinib",
                  "NSC228155")

pheatmap::pheatmap(top_APM_MHC_CII_genes_matrix[,drug_order_CII_1],
                   clustering_method = "ward.D2",
                   cluster_cols = FALSE,
                   color = hmcols_2,
                   breaks = myBreaks_2,
                   gaps_col = c(20),
                   cutree_rows = 4,
                   treeheight_row = 10,
                   fontsize = 6,
                   fontsize_col = 6,
                   height = 2,
                   width = 3,
                   border_color = NA,
                   file = "Marker_plots/top_APM_MHC_CII_genes.png")

pheatmap::pheatmap(top_APM_MHC_CII_genes_matrix[,drug_order_CII_1],
                   clustering_method = "ward.D2",
                   cluster_cols = FALSE,
                   color = hmcols_2,
                   breaks = myBreaks_2,
                   gaps_col = c(20),
                   cutree_rows = 4,
                   treeheight_row = 10,
                   fontsize = 6,
                   fontsize_col = 6,
                   height = 3,
                   width = 4,
                   border_color = NA,
                   file = "Marker_plots/top_APM_MHC_CII_genes_for_legend.png")

### Focus on a couple of interesting genes associated with tumor immunogenicity
treatment_diff_test_results %>%
  filter(term == "log(dose + 0.1)", 
         gene_short_name %in% c("JAK1","STAT1","STAT3"), 
         q_value < 0.05,
         drug == "CUDC-101") %>%
  dplyr::select(cell_line, drug, gene_short_name, q_value) %>%
  print(n = 30)

STAT_CUDC_cds <- cds[rowData(cds)$gene_short_name %in% c("STAT1","JAK1","STAT3"), 
                     (colData(cds)$drug %in% c("CUDC-101") & colData(cds)$dose == 10000) | colData(cds)$drug == "DMSO"]

colData(STAT_CUDC_cds)$treatment_id <- paste0(colData(STAT_CUDC_cds)$cell_line,"_",colData(STAT_CUDC_cds)$dose)

plot_percent_cells_positive(STAT_CUDC_cds, group_cells_by = "treatment_id", ncol = 1) +
  scale_fill_manual("PDCL", values = c("BT112_0" = "navy", "BT112_10000" = "navy","BT228_0" = "brown4", "BT228_10000" = "brown4","BT333_0" = "darkgreen", "BT333_10000" = "darkgreen")) +
  theme(text = element_text(size = 6),
        legend.position = "none",
        legend.key.width = unit(0.6,"line"), 
        legend.key.height = unit(0.6,"line")) +
  guides(guides(fill = guide_legend(override.aes = list(size=3))))
ggsave("Marker_plots/CUDC_effects_JAK_STAT_expression.png", dpi = 900, width = 1.1, height = 3.5)

treatment_diff_test_results %>%
  filter(term == "log(dose + 0.1)", 
         gene_short_name %in% c("HLA-DMA","HLA-DMB"), 
         q_value < 0.05,
         drug %in% c("MTX-211","Avitinib","BDTX-189","Panitumumab","Poziotinib")) %>%
  dplyr::select(cell_line, drug, gene_short_name, q_value) %>%
  print(n = 30)

HLADMA_cds <- cds[rowData(cds)$gene_short_name %in% c("HLA-DMA"), 
                  (colData(cds)$drug %in% c("MTX-211","Avitinib","BDTX-189") & colData(cds)$dose == 10000) | 
                    (colData(cds)$drug == "DMSO") | (colData(cds)$drug == "Panitumumab" & colData(cds)$dose == 1000)]

colData(HLADMA_cds)$treatment_id <- paste0(colData(HLADMA_cds)$cell_line,"_",colData(HLADMA_cds)$drug,"_",colData(HLADMA_cds)$dose)

colData(HLADMA_cds)$treatment_id <- factor(colData(HLADMA_cds)$treatment_id,
                                           levels = c("BT112_DMSO_0","BT112_BDTX-189_10000","BT112_Avitinib_10000", "BT112_Panitumumab_1000","BT112_MTX-211_10000",
                                                      "BT228_DMSO_0","BT228_BDTX-189_10000","BT228_Avitinib_10000", "BT228_Panitumumab_1000","BT228_MTX-211_10000",
                                                      "BT333_DMSO_0","BT333_BDTX-189_10000","BT333_Avitinib_10000", "BT333_Panitumumab_1000","BT333_MTX-211_10000"))

plot_percent_cells_positive(HLADMA_cds, group_cells_by = "treatment_id") +
  scale_fill_manual("PDCL", values = c("BT112_DMSO_0" = "navy","BT112_BDTX-189_10000" = "navy","BT112_Avitinib_10000" = "navy", "BT112_Panitumumab_1000" = "navy","BT112_MTX-211_10000" = "navy",
                                       "BT228_DMSO_0" = "brown4","BT228_BDTX-189_10000" = "brown4","BT228_Avitinib_10000" = "brown4", "BT228_Panitumumab_1000" = "brown4","BT228_MTX-211_10000" = "brown4",
                                       "BT333_DMSO_0" = "darkgreen","BT333_BDTX-189_10000" = "darkgreen","BT333_Avitinib_10000" = "darkgreen", "BT333_Panitumumab_1000" = "darkgreen","BT333_MTX-211_10000" = "darkgreen")) +
  theme(text = element_text(size = 6),
        legend.position = "none",
        legend.key.width = unit(0.6,"line"), 
        legend.key.height = unit(0.6,"line")) +
  guides(guides(fill = guide_legend(override.aes = list(size=3))))
ggsave("Marker_plots/EGFRi_effects_HLA-DMA_expression.png", dpi = 900, width = 2, height = 1.65)

HLADMB_cds <- cds[rowData(cds)$gene_short_name %in% c("HLA-DMB"), 
                  (colData(cds)$drug %in% c("Poziotinib","Avitinib") & colData(cds)$dose == 10000) | 
                    (colData(cds)$drug == "DMSO")]

colData(HLADMB_cds)$treatment_id <- paste0(colData(HLADMB_cds)$cell_line,"_",colData(HLADMB_cds)$drug,"_",colData(HLADMB_cds)$dose)

colData(HLADMB_cds)$treatment_id <- factor(colData(HLADMB_cds)$treatment_id,
                                           levels = c("BT112_DMSO_0","BT112_Avitinib_10000","BT112_Poziotinib_10000",
                                                      "BT228_DMSO_0","BT228_Avitinib_10000","BT228_Poziotinib_10000",
                                                      "BT333_DMSO_0","BT333_Avitinib_10000","BT333_Poziotinib_10000"))

plot_percent_cells_positive(HLADMB_cds, group_cells_by = "treatment_id", min_expr = 0) +
  scale_fill_manual("PDCL", values = c("BT112_DMSO_0" = "navy","BT112_Avitinib_10000" = "navy","BT112_Poziotinib_10000" = "navy",
                                       "BT228_DMSO_0" = "brown4","BT228_Avitinib_10000" = "brown4","BT228_Poziotinib_10000" = "brown4",
                                       "BT333_DMSO_0" = "darkgreen","BT333_Avitinib_10000" = "darkgreen","BT333_Poziotinib_10000" = "darkgreen")) +
  theme(text = element_text(size = 6),
        legend.position = "none",
        legend.key.width = unit(0.6,"line"), 
        legend.key.height = unit(0.6,"line")) +
  guides(guides(fill = guide_legend(override.aes = list(size=3))))
ggsave("Marker_plots/EGFRi_effects_HLA-DMB_expression.png", dpi = 900, width = 1.5, height = 1.65)
