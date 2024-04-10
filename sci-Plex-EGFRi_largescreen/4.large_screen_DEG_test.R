library(devtools)
library(parallel)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(monocle3)

cds <- readRDS("cds_large_screen.RDS")

cds.list <- list()
expressed_genes.list <- list()

for(PDCL in unique(colData(cds)$cell_line) %>% sort()){
  
  cds.list[[PDCL]] <- cds[,colData(cds)$cell_line == PDCL]
  cds.list[[PDCL]] <- detect_genes(cds.list[[PDCL]])
  cds.list[[PDCL]] <- estimate_size_factors(cds.list[[PDCL]])
  
  expressed_genes.list[[PDCL]] <- rowData(cds.list[[PDCL]])[rowData(cds.list[[PDCL]])$num_cells_expressed >= nrow(cds.list[[PDCL]])*0.01,]$id
  print(length(expressed_genes.list[[PDCL]]))
  
}

expressed_genes <- Reduce("intersect", expressed_genes.list) %>% unique()

drugs <- unique(colData(cds)$drug)
drugs <- drugs[!(drugs %in% c("DMSO","PBS","Media", "Panitumumab"))] %>% sort()

treatment_diff_test.list <- list()

for(PDCL in sort(unique(colData(cds)$cell_line))){
  
  treatment_diff_test.list[[PDCL]] <- list()
  
  for(treatment in drugs){
    
    cds_subset <- cds.list[[PDCL]][,colData(cds.list[[PDCL]])$drug %in% c("DMSO", treatment)]
    
    message("Running DEG analysis for ",unique(colData(cds_subset)$drug[colData(cds_subset)$drug != "DMSO"]), " treated ",unique(colData(cds_subset)$cell_line))
    
    treatment_diff_test.list[[PDCL]][[treatment]] <- fit_models(cds_subset[expressed_genes.list[[PDCL]],],
                                                                model_formula_str = "~log(dose + 0.1) + replicate",
                                                                cores = 1)
    treatment_diff_test.list[[PDCL]][[treatment]] <- coefficient_table(treatment_diff_test.list[[PDCL]][[treatment]]) %>% dplyr::select(-model,-model_summary)
    treatment_diff_test.list[[PDCL]][[treatment]]$drug <- rep(treatment, nrow(treatment_diff_test.list[[PDCL]][[treatment]]))
    
    message("Finished ",treatment)
    
  }
  
  message("Finished ",PDCL)
  
}

# saveRDS(treatment_diff_test.list, "treatment_diff_test.list.rds")
# treatment_diff_test.list <- readRDS("treatment_diff_test.list.rds")

### Test Panitumumab separately since doses are different
for(PDCL in sort(unique(colData(cds)$cell_line))){
  
  for(treatment in "Panitumumab"){
    
    cds_subset <- cds.list[[PDCL]][,colData(cds.list[[PDCL]])$drug %in% c("PBS", treatment) & colData(cds.list[[PDCL]])$cell_line == PDCL]
    treatment_diff_test.list[[PDCL]][[treatment]] <- fit_models(cds_subset[expressed_genes.list[[PDCL]],], model_formula_str = "~log(dose + 0.1) + replicate",
                                                                cores = 1)
    treatment_diff_test.list[[PDCL]][[treatment]] <- coefficient_table(treatment_diff_test.list[[PDCL]][[treatment]]) %>% dplyr::select(-model,-model_summary)
    treatment_diff_test.list[[PDCL]][[treatment]]$drug <- rep(treatment, nrow(treatment_diff_test.list[[PDCL]][[treatment]]))
    
    message("Finished ",treatment)
    
  }
  
}

### Randomly assign doses to vehicle cells as negative control tests
for(PDCL in sort(unique(colData(cds)$cell_line))){
  
  for(treatment in "Vehicle"){
    
    cds_subset <- cds.list[[PDCL]][,colData(cds.list[[PDCL]])$drug %in% c("DMSO","Media","PBS") & colData(cds.list[[PDCL]])$cell_line == PDCL]
    print(unique(colData(cds_subset)$dose))
    colData(cds_subset)$dose <- sapply(colData(cds_subset)$drug, function(x){
      if(x == "DMSO")return(0)
      if(x %in% c("Media","PBS"))return(sample(c(1,10,100,1000,10000),1))
      return(NA)
    })
    
    print(unique(colData(cds_subset)$dose))
    treatment_diff_test.list[[PDCL]][[treatment]] <- fit_models(cds_subset[expressed_genes.list[[PDCL]],], model_formula_str = "~log(dose + 0.1) + replicate",
                                                                cores = 1)
    treatment_diff_test.list[[PDCL]][[treatment]] <- coefficient_table(treatment_diff_test.list[[PDCL]][[treatment]]) %>% dplyr::select(-model,-model_summary)
    treatment_diff_test.list[[PDCL]][[treatment]]$drug <- rep(treatment, nrow(treatment_diff_test.list[[PDCL]][[treatment]]))
    
    message("Finished ",treatment)
    
  }
  
}

rm(cds_subset)

# saveRDS(treatment_diff_test.list, "full_treatment_diff_test.list.rds")
treatment_diff_test.list <- readRDS("full_treatment_diff_test.list.rds")

