library(tidyverse)
library(monocle3)
this.dir <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(this.dir))
rm(this.dir)

cds <- readRDS("sci-Plex-EGFRi-combination-PI3Ki-cds.RDS")

# adaptive gene list
adaptive_up_genes <- read_csv("adaptive_program.csv",
                              col_names = TRUE) %>%
  filter(grepl(pattern = "pregulated", gene_module))

adaptive_up_genes_id <- rowData(cds) %>% as.data.frame() %>%
  filter(gene_short_name %in% adaptive_up_genes$gene_short_name) %>%
  select(id)

write_csv(x = adaptive_up_genes_id, file = "adaptive_up_genes_id.csv", col_names = F)

colData(cds)$replicate <- case_when(
  colData(cds)$hash_plate == "01" ~ "1",
  colData(cds)$hash_plate == "02" ~ "2"
)

colData(cds)$sample_key <- paste0(colData(cds)$drug2, "_", 
                                  colData(cds)$dose2, "_",
                                  colData(cds)$drug1,"_",
                                  colData(cds)$dose1)

# printing matrices for each cell line
Matrix::writeMM(t(exprs(cds[adaptive_up_genes_id %>% pull(),])), 
                file = "EGFRi_PI3Ki_UMI_count_filt_transpose.matrix")
write_csv(colData(cds) %>% as.data.frame(), 
          file = "EGFRi_PI3Ki_colData.csv")
write_csv(rowData(cds[adaptive_up_genes_id %>% pull(),]) %>% as.data.frame(), 
          file = "EGFRi_PI3Ki_rowData.csv")

dim(cds[adaptive_up_genes_id %>% pull(),])


