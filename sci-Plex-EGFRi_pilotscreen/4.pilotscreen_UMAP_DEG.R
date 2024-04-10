library(devtools)
library(ggplot2)
library(ggridges)
library(gridExtra)
library(plyr)
library(dplyr)
library(tidyr)
library(monocle3)

### Point to sci-plex repo below
source("~/Documents/github_repos/sci-plex/bin/cell_cycle.R")
source("~/Documents/github_repos/sci-plex/bin/loadGSCSafe.R")
source("~/Documents/github_repos/sci-plex/bin/GSA_helper_functions.R")
cc.genes = readRDS("~/Documents/github_repos/sci-plex/bin/cc.genes.RDS")

#### Functions ####
### Define Helper function for downstream Gene Set Enrichment Analysis
replace_gene_names_vec <- function(input_vec, name_vec, retain_inds = c(-1,-2)) {
  temp <- merge(name_vec, input_vec, by="row.names")
  temp2 <- temp[,retain_inds]
  names(temp2) <- temp[,2]
  return(temp2)
}

loadGSCSafe <- function (file, type = "auto", addInfo, sep="\t", encoding="latin1") 
{
  if (missing(addInfo)) {
    addUserInfo <- "skip"
    addInfo <- "none"
  }
  else {
    addUserInfo <- "yes"
  }
  tmp <- try(type <- match.arg(type, c("auto", "gmt", "sbml", 
                                       "sif", "data.frame"), several.ok = FALSE), silent = TRUE)
  if (class(tmp) == "try-error") {
    stop("argument type set to unknown value")
  }
  if (type == "auto") {
    if (class(file) == "character") {
      tmp <- unlist(strsplit(file, "\\."))
      type <- tolower(tmp[length(tmp)])
      if (!type %in% c("gmt", "sif", "sbml", "xml")) 
        stop(paste("can not handle .", type, " file extension, read manually using e.g. read.delim() and load as data.frame", 
                   sep = ""))
    }
    else {
      type <- "data.frame"
    }
  }
  if (type == "gmt") {
    con <- file(file, encoding=encoding)
    tmp <- try(suppressWarnings(open(con)), silent = TRUE)
    if (class(tmp) == "try-error") 
      stop("file could not be read")
    if (addUserInfo == "skip") 
      addInfo <- vector()
    gscList <- list()
    i <- 1
    tmp <- try(suppressWarnings(while (length(l <- scan(con, 
                                                        nlines = 1, what = "character", quiet = T, sep=sep)) > 0) {
      if (addUserInfo == "skip") 
        addInfo <- rbind(addInfo, l[1:2])
      tmp <- l[3:length(l)]
      gscList[[l[1]]] <- unique(tmp[tmp != "" & tmp != 
                                      " " & !is.na(tmp)])
      i <- i + 1
    }), silent = TRUE)
    if (class(tmp) == "try-error") 
      stop("file could not be read")
    close(con)
    gsc <- gscList[!duplicated(names(gscList))]
    if (addUserInfo == "skip") 
      addInfo <- unique(addInfo)
  }
  else if (type %in% c("sbml", "xml")) {
    require(rsbml)
    tmp <- try(sbml <- rsbml_read(file))
    if (class(tmp) == "try-error") {
      stop("file could not be read by rsbml_read()")
    }
    gsc <- list()
    for (iReaction in 1:length(reactions(model(sbml)))) {
      metIDs <- names(c(reactants(reactions(model(sbml))[[iReaction]]), 
                        products(reactions(model(sbml))[[iReaction]])))
      geneIDs <- names(modifiers(reactions(model(sbml))[[iReaction]]))
      if (length(geneIDs) > 0) {
        geneNames <- rep(NA, length(geneIDs))
        for (iGene in 1:length(geneIDs)) {
          geneNames[iGene] <- name(species(model(sbml))[[geneIDs[iGene]]])
        }
        for (iMet in 1:length(metIDs)) {
          gsc[[metIDs[iMet]]] <- c(gsc[[metIDs[iMet]]], 
                                   geneNames)
        }
      }
    }
    if (length(gsc) == 0) {
      stop("no gene association found")
    }
    else {
      for (iMet in 1:length(gsc)) {
        tmp1 <- name(species(model(sbml))[[names(gsc)[iMet]]])
        tmp2 <- compartment(species(model(sbml))[[names(gsc)[iMet]]])
        names(gsc)[iMet] <- paste(tmp1, " (", tmp2, ")", 
                                  sep = "")
      }
    }
  }
  else if (type == "sif") {
    tmp <- try(gsc <- as.data.frame(read.delim(file, header = FALSE, 
                                               quote = "", as.is = TRUE), stringsAsFactors = FALSE), 
               silent = TRUE)
    if (class(tmp) == "try-error") {
      stop("argument file could not be read and converted into a data.frame")
    }
    if (ncol(gsc) != 3) {
      stop("sif file should contain three columns")
    }
    if (addUserInfo == "skip") 
      addInfo <- gsc[, c(1, 2)]
    gsc <- gsc[, c(3, 1)]
    tmp <- nrow(gsc)
    gsc <- unique(gsc)
    geneSets <- unique(gsc[, 2])
    gscList <- list()
    for (iGeneSet in 1:length(geneSets)) {
      gscList[[iGeneSet]] <- gsc[gsc[, 2] == geneSets[iGeneSet], 
                                 1]
    }
    names(gscList) <- geneSets
    gsc <- gscList
  }
  else if (type == "data.frame") {
    tmp <- try(gsc <- as.data.frame(file, stringsAsFactors = FALSE), 
               silent = TRUE)
    if (class(tmp) == "try-error") {
      stop("argument file could not be converted into a data.frame")
    }
    for (i in 1:ncol(gsc)) {
      gsc[, i] <- as.character(gsc[, i])
    }
    if (ncol(gsc) != 2) {
      stop("argument file has to contain exactly two columns")
    }
    tmp <- nrow(gsc)
    gsc <- unique(gsc)
    geneSets <- unique(gsc[, 2])
    gscList <- list()
    for (iGeneSet in 1:length(geneSets)) {
      gscList[[iGeneSet]] <- gsc[gsc[, 2] == geneSets[iGeneSet], 
                                 1]
    }
    names(gscList) <- geneSets
    gsc <- gscList
  }
  if (addUserInfo == "yes") {
    tmp <- try(addInfo <- as.data.frame(addInfo, stringsAsFactors = FALSE), 
               silent = TRUE)
    if (class(tmp) == "try-error") {
      stop("failed to convert additional info in argument 'addInfo' into a data.frame")
    }
  }
  
  addInfo <- as.data.frame(addInfo)
  
  if (class(addInfo) == "data.frame") {
    if (ncol(addInfo) != 2) 
      stop("additional info in argument 'file' or 'addInfo' has to contain 2 columns")
    tmp <- nrow(addInfo)
    addInfo <- unique(addInfo[addInfo[, 1] %in% names(gsc), 
    ])
  }
  else {
  }
  res <- list(gsc, addInfo)
  names(res) <- c("gsc", "addInfo")
  class(res) <- "GSC"
  return(res)
}

collect_gsa_hyper_results_clusters <- function (genes_list, clusters, gsc) 
{
  gene_universe <- unique(as.character(genes_list))
  gsa_results <- list()
  cluster_ids <- unique(clusters)
  for (cluster in cluster_ids) {
    cluster_genes <- unique(names(clusters[clusters == cluster]))
    gsaRes <- runGSAhyper(cluster_genes, gsc = gsc, universe = gene_universe, 
                          adjMethod = "BH")
    gsa_results[[cluster]] <- gsaRes
  }
  names(gsa_results) <- cluster_ids
  gsa_results
}

gsea_bar_plots <- function(GSAhyper_list, qval_cutoff, pattern, width, height, sample, gsc){
  
  for(cluster in names(GSAhyper_list)){
    
    print(cluster)
    
    GSAhyper_df <- as.data.frame(GSAhyper_list[[cluster]]$p.adj)
    GSAhyper_df$gene_set <- row.names(GSAhyper_df)
    colnames(GSAhyper_df) <- c("qval","gene_set")
    
    if(is.null(pattern) == FALSE){
      GSAhyper_df$gene_set <- stringr::str_replace(string = GSAhyper_df$gene_set, pattern = pattern, replace = "")
    }
    
    GSAhyper_df_cutoff <- GSAhyper_df %>% filter(qval < qval_cutoff) %>% arrange(desc(qval)) %>% 
      mutate(gene_set = factor(gene_set, levels = gene_set))
    
    plot_title <- paste0(sample,"_",as.character(cluster),"_",gsc,".png")
    print(plot_title)
    
    ggplot(GSAhyper_df_cutoff, aes(x = gene_set, y = -log10(qval))) + 
      geom_bar(stat = "identity") + 
      coord_flip() +
      theme_classic(base_size = 8)
    ggsave(plot_title, width = width, height = height)   
    
  }
  
}

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

# Set DelayedArray Parameters
DelayedArray:::set_verbose_block_processing(TRUE)
options(DelayedArray.block.size=1000e7)
options(stringsAsFactors = FALSE)

cds <- readRDS("cds_pilot_screen.RDS")
unique(colData(cds)$treatment)

ggplot(colData(cds) %>% 
         as.data.frame() %>% 
         arrange(desc(n.umi)) %>% 
         mutate(cell_rank = dplyr::row_number()), 
       aes(x= log10(cell_rank), y = log10(n.umi))) +
  geom_point(size = 0.1, stroke = 0) +
  geom_hline(yintercept = log10(100), color = "red", size = 0.1) +
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 6)) 
ggsave("QC_plots/Kneeplot.png", dpi = 600, height = 1, width = 1.5)

ggplot(colData(cds) %>% 
         as.data.frame(),
       aes(x = log10(n.umi))) +
  geom_density(size = 0.25) +
  geom_vline(xintercept = log10(100), color = "dimgrey", size = 0.1) +
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 6)) 
ggsave("QC_plots/UMI_count_distribution.png", dpi = 600, height = 1, width = 1.5)

ggplot(colData(cds) %>% 
         as.data.frame(),
       aes(x = log10(hash_umis))) +
  geom_density(size = 0.25) +
  geom_vline(xintercept = log10(10), color = "dimgrey", size = 0.1) +
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 6)) 
ggsave("QC_plots/Hash_UMI_count_distribution.png", dpi = 600, height = 1, width = 1.5)

ggplot(colData(cds) %>% 
         as.data.frame(),
       aes(x = log2(top_to_second_best_ratio))) +
  geom_density(size = 0.25) +
  geom_vline(xintercept = log2(5), color = "dimgrey", size = 0.1) +
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 6)) 
ggsave("QC_plots/Top_to_second_best_ratio_distribution.png", dpi = 600, height = 1, width = 1.5)

dim(cds)
dim(cds[,colData(cds)$hash_umis >= 10 &
          colData(cds)$top_to_second_best_ratio >= 5])

cds <- cds[,colData(cds)$hash_umis >= 10 &
             colData(cds)$top_to_second_best_ratio >= 5]

colData(cds)$vehicle <- sapply(colData(cds)$treatment, function(x){ifelse(x %in% c("DMSO"), TRUE, FALSE)})
cds <- estimate_size_factors(cds)
cds <- detect_genes(cds)

colData(cds)$dose_character <- factor(colData(cds)$dose, levels = c("0","5","10","100","500","1000","5000","10000"))

# Calculate the proliferation index for every cell 
cds = estimate_cell_cycle(cds, g1s_markers = cc.genes$s.genes, g2m_markers = cc.genes$g2m.genes)

expressed_genes_per_drug = colData(cds) %>%
  as.data.frame() %>%
  group_by(treatment, Cell.Line) %>%
  nest() %>%
  mutate(fraction_genes = purrr::map(data, .f = function(pdata_subset, cds) {
    cds_subset = cds[,as.character(pdata_subset$cell)]
    cds_subset = detect_genes(cds_subset)
    tibble(id = rowData(cds_subset)$id,
           gene_short_name = rowData(cds_subset)$gene_short_name,
           num_cells_expressed = rowData(cds_subset)$num_cells_expressed,
           fraction_cells_expressed = rowData(cds_subset)$num_cells_expressed / ncol(cds_subset))
  }, cds))

expressed_genes_per_drug = expressed_genes_per_drug %>%
  unnest(fraction_genes) %>%
  dplyr::select(everything(),-data)

expressed_genes_per_drug_id = expressed_genes_per_drug %>%
  filter(fraction_cells_expressed >= 0.01) %>%
  ungroup() %>%
  pull(id) %>%
  unique()

length(expressed_genes_per_drug_id)

expressed_genes <- row.names(fData(cds)[Matrix::rowSums(exprs(cds) > 0) > dim(cds)[2]*0.01 ,])

length(expressed_genes)

### Do a first pass to filter missassigned hash cells

cds <- preprocess_cds(cds,
                      method = "PCA",
                      num_dim = 25,
                      norm_method = "log",
                      use_genes = expressed_genes)

cds <- reduce_dimension(cds,
                        max_components = 2,
                        preprocess_method = "PCA",
                        reduction_method = "UMAP",
                        umap.metric = "cosine",
                        umap.n_neighbors = 15,
                        umap.min_dist = 0.1,
                        umap.fast_sgd=FALSE,
                        cores=1,
                        verbose = T)

colData(cds)$UMAP1 <- reducedDims(cds)[["UMAP"]][,1]
colData(cds)$UMAP2 <- reducedDims(cds)[["UMAP"]][,2]

cds <- cluster_cells(cds,
                     resolution = 5e-4,
                     reduction_method = "UMAP",
                     num_iter = 1,
                     random_seed = 2016L)

colData(cds)$Cluster <- clusters(cds, reduction_method = "UMAP")
length(unique(colData(cds)$Cluster))
colData(cds)$Partition <- partitions(cds, reduction_method = "UMAP")
length(unique(colData(cds)$Partition))

plot_cells(cds, color_cells_by = "Partition",label_cell_groups = FALSE, cell_size = 0.5) +
  theme_void() +
  facet_wrap(~Cell.Line, ncol = 2) +
  theme(legend.position = "right",
        text = element_text(size = 6),
        legend.key.width = unit(0.1,"line"),
        legend.key.height = unit(0.1,"line")) +
  scale_color_brewer("Partitions",palette = "Paired") +
  guides(guides(fill = guide_legend(override.aes = list(size=1)))) +
  labs(color="EGFRi")
ggsave("UMAPs/UMAP_by_cluster_faceted_by_Cell_line_for_filtering.png", dpi = 600, height = 3, width = 4)

### Below will vary slightly depending on clustering, keeping hash cells that are in clusters made or priamrily of one cell line 
A172_cells <- colData(cds) %>% as.data.frame() %>% filter(Cell.Line == "A172", Partition %in% c("3")) %>% pull(cell)
T98G_cells <- colData(cds) %>% as.data.frame() %>% filter(Cell.Line == "T98G", Partition %in% c("4","6")) %>% pull(cell)
U87MG_cells <- colData(cds) %>% as.data.frame() %>% filter(Cell.Line == "U87MG", Partition %in% c("5")) %>% pull(cell)
BT333_cells <- colData(cds) %>% as.data.frame() %>% filter(Cell.Line == "BT333", Partition %in% c("1","2")) %>% pull(cell)

GBM_cells_filtered <- c(A172_cells,T98G_cells,U87MG_cells,BT333_cells)
# saveRDS(GBM_cells_filtered, "GBM_cells_filtered.rds")

targets <- unique(colData(cds)$treatment)
targets <- targets[!(targets %in% c("DMSO"))] %>% sort()

diff_test.list <- list()

for(cell_line in unique(colData(cds)$Cell.Line)){
  
  diff_test.list[[cell_line]] <- list()
  
  for(egfri in targets){
    
    cds_subset <- cds[,colData(cds)$Cell.Line == cell_line]
    cds_subset <- cds_subset[,colData(cds_subset)$treatment == egfri | colData(cds_subset)$vehicle == TRUE]
    diff_test <- fit_models(cds_subset[expressed_genes,], model_formula_str = "~log(dose + 0.1)", cores = 1)
    diff_test <- coefficient_table(diff_test) %>% dplyr::select(-model, -model_summary)
    diff_test$treatment <- rep(egfri, dim(diff_test)[1])
    diff_test.list[[cell_line]][[egfri]] <- diff_test
    rm(diff_test, cds_subset)
    
    message("finished ", egfri)
    
  }
  
  message("Done with all tests in ", cell_line)
}

# saveRDS(diff_test.list, "diff_test.list.rds")
# diff_test.list <- readRDS("diff_test.list.rds")

#### Post cleaning dataset and DEG test
GBM_cells_filtered <- readRDS("GBM_cells_filtered.rds")
diff_test.list <- readRDS("diff_test.list.rds")

for(cell_line in unique(colData(cds)$Cell.Line)){
  
  diff_test.list[[cell_line]] <- do.call("rbind",diff_test.list[[cell_line]])
  diff_test.list[[cell_line]]$Cell.Line <- rep(cell_line, dim(diff_test.list[[cell_line]])[1])
}

diff_test_result <- do.call("rbind", diff_test.list)
diff_test_result$q_value <- p.adjust(diff_test_result$p_value, method = "BH")

EGFRi_de_sig_genes<- diff_test_result %>% 
  filter(term == "log(dose + 0.1)", q_value < 0.001, abs(normalized_effect)>0.05) %>%
  pull(id) %>%
  unique()

length(EGFRi_de_sig_genes)

cds <- preprocess_cds(cds,
                      method = "PCA",
                      num_dim = 25,
                      norm_method = "log",
                      use_genes = EGFRi_de_sig_genes)

cds <- align_cds(cds, 
                 residual_model_formula_str = "~Cell.Line + batch",
                 alignment_group = "Cell.Line")

cds <- reduce_dimension(cds,
                        max_components = 2,
                        preprocess_method = "Aligned",
                        reduction_method = "UMAP",
                        umap.metric = "cosine",
                        umap.n_neighbors = 15,
                        umap.min_dist = 0.1,
                        umap.fast_sgd=FALSE,
                        cores=1,
                        verbose = T)

colData(cds)$UMAP1 <- reducedDims(cds)[["UMAP"]][,1]
colData(cds)$UMAP2 <- reducedDims(cds)[["UMAP"]][,2]

cds <- cluster_cells(cds,
                     resolution = 7.5e-5,
                     reduction_method = "Aligned",
                     num_iter = 1,
                     random_seed = 2016L)

colData(cds)$Cluster <- clusters(cds, reduction_method = "Aligned")
length(unique(colData(cds)$Cluster))
colData(cds)$Partition <- partitions(cds, reduction_method = "Aligned")
length(unique(colData(cds)$Partition))

colData(cds)$treatment <- factor(colData(cds)$treatment, 
                                 levels = c("DMSO","Afatinib","Brigatinib",
                                            "CUDC-101","EAI045","Neratinib","Osimertinib"))

plot_cells(cds, color_cells_by = "treatment",label_cell_groups = FALSE, cell_size = 0.1) +
  theme_void() +
  theme(legend.position = "right",
        text = element_text(size = 9),
        legend.text = element_text(size = 6),
        legend.key.width = unit(0.1,"line"),
        legend.key.height = unit(0.1,"line")) +
  scale_color_manual("Exposure",
                     values = c("Afatinib" = "firebrick1","Brigatinib" = "navy",
                                "CUDC-101" = "forestgreen","EAI045" = "deepskyblue3",
                                "Neratinib" = "brown4","Osimertinib" = "darkorange1",
                                "DMSO" = "grey80")) +
  guides(color = guide_legend(override.aes = list(size=2)))
ggsave("UMAPs/Fig1_UMAP_by_treatment.png", dpi = 600, height = 1, width = 2)

plot_cells(cds, color_cells_by = "Cell.Line",label_cell_groups = FALSE, cell_size = 0.1) +
  theme_void() +
  theme(legend.position = "right",
        text = element_text(size = 6),
        legend.key.width = unit(0.1,"line"),
        legend.key.height = unit(0.1,"line")) +
  scale_color_manual("GBM line",
                     values = c("A172" = "firebrick1",
                                "T98G" = "forestgreen",
                                "U87MG" = "brown4",
                                "BT333" = "navy")) +
  guides(color = guide_legend(override.aes = list(size=2)))
ggsave("UMAPs/Fig1_UMAP_by_GBM_line.png", dpi = 600, height = 1, width = 2)

plot_cells(cds, color_cells_by = "log10_dose",label_cell_groups = FALSE, cell_size = 0.1) +
  theme_void() +
  theme(legend.position = "right",
        text = element_text(size = 9),
        legend.text = element_text(size = 6),
        legend.key.width = unit(0.5,"line"),
        legend.key.height = unit(0.5,"line")) +
  viridis::scale_color_viridis("log10(dose)\n(nM)", option = "magma") 
ggsave("UMAPs/Fig1_UMAP_by_dose.png", dpi = 600, height = 1, width = 2.1)

plot_cells(cds, color_cells_by = "Cluster",label_cell_groups = FALSE, cell_size = 0.1) +
  theme_void() +
  theme(legend.position = "right",
        text = element_text(size = 9),
        legend.text = element_text(size = 6),
        legend.key.width = unit(0.1,"line"),
        legend.key.height = unit(0.1,"line")) +
  scale_color_brewer("Cluster", palette = "Paired") +
  guides(color = guide_legend(override.aes = list(size=2)))
ggsave("UMAPs/Fig1_UMAP_by_Cluster.png", dpi = 600, height = 1, width = 1.8)

plot_cells(cds, color_cells_by = "treatment",label_cell_groups = FALSE, cell_size = 0.01) +
  theme_void() +
  facet_wrap(~treatment, ncol = 3) +
  theme(legend.position = "right",
        text = element_text(size = 6),
        legend.key.width = unit(0.1,"line"),
        legend.key.height = unit(0.1,"line")) +
  scale_color_manual("EGFRi",
                     values = c("Afatinib" = "firebrick2","Brigatinib" = "navy",
                                "CUDC-101" = "forestgreen","EAI045" = "deepskyblue3",
                                "Neratinib" = "brown4","Osimertinib" = "darkorange3",
                                "DMSO" = "grey80")) +
  guides(guides(color = guide_legend(override.aes = list(size=2)))) 
ggsave("UMAPs/Fig1_UMAP_by_treatment_faceted.png", dpi = 600, height = 1, width = 2)

colData(cds) %>%
  as.data.frame() %>%
  filter(dose %in% c(0,10000)) %>%
  group_by(Cell.Line, treatment, Cluster) %>%
  summarise(cells_per_cluster = n()) %>%
  group_by(Cell.Line, treatment) %>% 
  mutate(freq_per_cluster = cells_per_cluster/sum(cells_per_cluster)) %>%
  ungroup() %>%
  ggplot() +
  geom_bar(aes(x = treatment, y = freq_per_cluster, fill = Cluster), stat = "identity", color = "black", size = 0.1) +
  facet_wrap(~Cell.Line, ncol = 4) +
  monocle3:::monocle_theme_opts() +
  scale_fill_brewer("PCA\ncluster", palette = "Paired") +
  scale_x_discrete(limits=rev) +
  theme(text = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        legend.position = "none",
        legend.key.width = unit(0.1,"line"),
        legend.key.height = unit(0.1,"line")) +
  xlab("Exposure") +
  ylab("Frequency per Cluster") +
  coord_flip() +
  guides(guides(fill = guide_legend(override.aes = list(size=2)))) 
ggsave("UMAPs/Fig1_UMAP_by_treatment_cluster_distribution.png", dpi = 600, height = 1.3, width = 3.25)

ggplot(colData(cds) %>% as.data.frame(), aes(x = Cluster, y = proliferation_index)) +
  geom_violin() +
  facet_wrap(~Cell.Line) +
  stat_summary(fun=mean, geom="point", size = 1, stroke = 0) +
  monocle3:::monocle_theme_opts()

##
cluster_diff_test.list <- list()
for(cluster in unique(colData(cds)$Cluster) %>% sort()){
  
  message("Identifying differentially expressed genes for cluster ", cluster)
  cds_subset <- cds
  colData(cds_subset)$Cluster_tested_ <- sapply(colData(cds_subset)$Cluster, function(x){ifelse(x == cluster, cluster, "Out")})
  colData(cds_subset)$Cluster_tested_ <- factor(colData(cds_subset)$Cluster_tested_, levels = c("Out", cluster))
  
  cluster_diff_test.list[[cluster]] <- fit_models(cds_subset[expressed_genes_per_drug_id,], model_formula_str = "~Cluster_tested_ + Cell.Line")
  cluster_diff_test.list[[cluster]] <- coefficient_table(cluster_diff_test.list[[cluster]]) %>% dplyr::select(-model,-model_summary)
  cluster_diff_test.list[[cluster]]$Cluster <- rep(cluster, nrow(cluster_diff_test.list[[cluster]]))
  
  rm(cds_subset)
}

cluster_diff_test_result <- do.call("rbind",cluster_diff_test.list)

# saveRDS(cluster_diff_test_result, "cluster_diff_test_result_cluster_vs_all.rds")
cluster_diff_test_result <- readRDS("cluster_diff_test_result_cluster_vs_all.rds")
# saveRDS(cds, "Fig1_postprocessed_cds.rds")

cluster_diff_test_result %>%
  filter(grepl("Cluster", term)) %>%
  ggplot() +
  geom_point(aes(x = normalized_effect, y = -log10(q_value), color = q_value < 0.001 & abs(normalized_effect) > 0.5), size = 0.5, stroke = 0.01) +
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 9),
        legend.position = "none",
        legend.key.width = unit(0.1,"line"),
        legend.key.height = unit(0.1,"line")) +
  scale_color_manual("FDR < 1%", values = c("TRUE" = "firebrick3", "FALSE" = "black")) +
  guides(guides(fill = guide_legend(override.aes = list(size=2)))) 

cluster_up_sig_genes <- cluster_diff_test_result %>%
  filter(grepl("Cluster", term), q_value < 0.001, normalized_effect > 0.5) %>%
  dplyr::select(Cluster,id) %>%
  distinct()

cluster_dn_sig_genes <- cluster_diff_test_result %>%
  filter(grepl("Cluster", term), q_value < 0.001, normalized_effect < -0.5) %>%
  dplyr::select(Cluster,id) %>%
  distinct()

hallmarksGSC <- loadGSCSafe(file="h.all.v6.0.symbols.gmt")
OncogenicSignaturesGSC <- loadGSCSafe(file="c6.all.v6.0.OncogenicSignatures.symbols.gmt")

Ensembl_GSAlist <- rowData(cds) %>% 
  as.data.frame() %>%
  filter(id %in% expressed_genes_per_drug_id) %>%
  pull(gene_short_name)
Ensembl_GSAlist <- as.matrix(Ensembl_GSAlist)

rownames(Ensembl_GSAlist)<- rowData(cds) %>% 
  as.data.frame() %>%
  filter(id %in% expressed_genes_per_drug_id) %>%
  row.names()

colnames(Ensembl_GSAlist) <- c("gene_short_name")
Ensembl_GSAlist<-Ensembl_GSAlist[,1]
Ensembl_GSAlist<-toupper(Ensembl_GSAlist)
length(Ensembl_GSAlist)

cluster_up_sig_genes_id <- cluster_up_sig_genes$Cluster %>% as.matrix()
row.names(cluster_up_sig_genes_id) <- cluster_up_sig_genes$id

cluster_dn_sig_genes_id <- cluster_dn_sig_genes$Cluster %>% as.matrix()
row.names(cluster_dn_sig_genes_id) <- cluster_dn_sig_genes$id

dir.create("GSA/Cluster_GSEA/Up_genes/OncogenicSignatures")
gsea_bar_plots(OncogenicSignatures_GSAhyper_up, qval_cutoff = 0.05, pattern = "OncogenicSignatures_", 
               width = 8, height = 10, sample = "GSA/Cluster_GSEA/Up_genes/OncogenicSignatures/Cluster_module", gsc = "Oncogenic_signatures")

OncogenicSignatures_GSAhyper_dn <- collect_gsa_hyper_results_clusters(Ensembl_GSAlist,
                                                                      replace_gene_names_vec(cluster_dn_sig_genes_id,
                                                                                             Ensembl_GSAlist),
                                                                      OncogenicSignaturesGSC)

dir.create("GSA/Cluster_GSEA/Dn_genes/OncogenicSignatures")
gsea_bar_plots(OncogenicSignatures_GSAhyper_dn, qval_cutoff = 0.05, pattern = "OncogenicSignatures_", 
               width = 8, height = 10, sample = "GSA/Cluster_GSEA/Dn_genes/OncogenicSignatures/Cluster_module", gsc = "Oncogenic_signatures")

#### Drug module (signatures) -- see DEG_signatures.R to generate signature list

pilot_screen_EGFRi_intersection_sigantures_list <- readRDS("pilot_screen_EGFRi_intersection_sigantures_list.RDS")

colData(cds)$Shared_EGFR <- calculate_aggregate_expression_score(cds, signature_genes = pilot_screen_EGFRi_intersection_sigantures_list[["Afat/Brigat/CUDC/Nerat/Osimert"]])
colData(cds)$CUDC <- calculate_aggregate_expression_score(cds, signature_genes = pilot_screen_EGFRi_intersection_sigantures_list[["CUDC"]])
colData(cds)$CUDC_Osimert <- calculate_aggregate_expression_score(cds, signature_genes = pilot_screen_EGFRi_intersection_sigantures_list[["CUDC/Osimert"]])
colData(cds)$Brigat_CUDC <- calculate_aggregate_expression_score(cds, signature_genes = pilot_screen_EGFRi_intersection_sigantures_list[["Brigat/CUDC"]])
colData(cds)$Afat_CUDC_Osimert <- calculate_aggregate_expression_score(cds, signature_genes = pilot_screen_EGFRi_intersection_sigantures_list[["Afat/CUDC/Osimert"]])

dmso_median_shared_egfr <- colData(cds) %>% as.data.frame() %>% filter(treatment == "DMSO") %>% pull(Shared_EGFR) %>% median()

ggplot(colData(cds) %>% as.data.frame() %>% filter(dose %in% c(0,10000), treatment != "EAI045", Cell.Line == "BT333"), 
       aes(x = factor(treatment, levels = c("DMSO","Afatinib","Brigatinib","CUDC-101","Neratinib","Osimertinib")), y = Shared_EGFR, fill = treatment)) +
  geom_violin(size = 0.1, scale = "width") +
  geom_hline(yintercept = dmso_median_shared_egfr, size = 0.1, linetype = "dashed") +
  stat_summary(fun.y = "mean", geom = "point", stroke = 0, size = 1) +
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 9),
        axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.position = "none",
        legend.key.width = unit(0.5,"line"),
        legend.key.height = unit(0.5,"line"),
        plot.title = element_text(hjust = 0.5, size = 9)) +
  ggtitle("Shared EGFR\nmodule") +
  xlab("Exposure") +
  ylab("Aggregate Expression") +
  scale_fill_manual("Exposure",
                    values = c("Afatinib" = "firebrick2","Brigatinib" = "navy",
                               "CUDC-101" = "forestgreen",
                               "Neratinib" = "brown4","Osimertinib" = "darkorange3",
                               "DMSO" = "grey80")) +
  guides(fill = guide_legend(override.aes = list(size=1))) 
ggsave("Signature_violins/Shared_EGFR_signature_by_top_dose.png", width = 1.25, height = 1.75, dpi = 900)

ggplot(colData(cds) %>% as.data.frame() %>% filter(treatment %in% c("DMSO","Osimertinib"), Cell.Line == "BT333"), 
       aes(x = factor(as.character(dose), levels = c("0","5","10","100","500","1000","5000","10000")), y = Shared_EGFR, 
           fill = factor(as.character(dose), levels = c("0","5","10","100","500","1000","5000","10000")))) +
  geom_violin(size = 0.1, scale = "width") +
  geom_hline(yintercept = dmso_median_shared_egfr, size = 0.1, linetype = "dashed") +
  stat_summary(fun.y = "mean", geom = "point", stroke = 0, size = 1) +
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 9),
        axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.position = "none",
        legend.key.width = unit(0.5,"line"),
        legend.key.height = unit(0.5,"line"),
        plot.title = element_text(hjust = 0.5, size = 9)) +
  ggtitle("Shared EGFR\nmodule") +
  xlab("Osimertinib\n(nM)") +
  ylab("Aggregate Expression") +
  viridis::scale_fill_viridis("Dose (nM)", discrete = "TRUE", option = "magma") +
  guides(fill = guide_legend(override.aes = list(size=1))) 
ggsave("Signature_violins/Shared_EGFR_signature_by_Osimertinib.png", width = 1.25, height = 1.75, dpi = 900)

dmso_median_cudc <- colData(cds) %>% as.data.frame() %>% filter(treatment == "DMSO") %>% pull(CUDC) %>% median()

ggplot(colData(cds) %>% as.data.frame() %>% filter(dose %in% c(0,10000), treatment != "EAI045", Cell.Line == "BT333"), 
       aes(x = factor(treatment, levels = c("DMSO","Afatinib","Brigatinib","CUDC-101","Neratinib","Osimertinib")), y = CUDC, fill = treatment)) +
  geom_violin(size = 0.1, scale = "width") +
  geom_hline(yintercept = dmso_median_cudc, size = 0.1, linetype = "dashed") +
  stat_summary(fun.y = "mean", geom = "point", stroke = 0, size = 1) +
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 9),
        axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.position = "none",
        legend.key.width = unit(0.5,"line"),
        legend.key.height = unit(0.5,"line"),
        plot.title = element_text(hjust = 0.5, size = 9)) +
  ggtitle("CUDC-101\nmodule") +
  xlab("Exposure") +
  ylab("Aggregate Expression") +
  scale_fill_manual("Exposure",
                    values = c("Afatinib" = "firebrick2","Brigatinib" = "navy",
                               "CUDC-101" = "forestgreen",
                               "Neratinib" = "brown4","Osimertinib" = "darkorange3",
                               "DMSO" = "grey80")) +
  guides(fill = guide_legend(override.aes = list(size=1))) 
ggsave("Signature_violins/CUDC_signature_by_top_dose.png", width = 1.25, height = 1.75, dpi = 900)

dmso_median_cudc_osimert <- colData(cds) %>% as.data.frame() %>% filter(treatment == "DMSO") %>% pull(CUDC_Osimert) %>% median()

ggplot(colData(cds) %>% as.data.frame() %>% filter(dose %in% c(0,10000), treatment != "EAI045", Cell.Line == "BT333"), 
       aes(x = factor(treatment, levels = c("DMSO","Afatinib","Brigatinib","CUDC-101","Neratinib","Osimertinib")), y = CUDC_Osimert, fill = treatment)) +
  geom_violin(size = 0.1, scale = "width") +
  geom_hline(yintercept = dmso_median_cudc_osimert, size = 0.1, linetype = "dashed") +
  stat_summary(fun.y = "mean", geom = "point", stroke = 0, size = 1) +
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 9),
        axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.position = "none",
        legend.key.width = unit(0.5,"line"),
        legend.key.height = unit(0.5,"line"),
        plot.title = element_text(hjust = 0.5, size = 9)) +
  ggtitle("CUDC-101/Osim\nmodule") +
  xlab("Exposure") +
  ylab("Aggregate Expression") +
  scale_fill_manual("Exposure",
                    values = c("Afatinib" = "firebrick2","Brigatinib" = "navy",
                               "CUDC-101" = "forestgreen",
                               "Neratinib" = "brown4","Osimertinib" = "darkorange3",
                               "DMSO" = "grey80")) +
  guides(fill = guide_legend(override.aes = list(size=1))) 
ggsave("Signature_violins/CUDC_Osimertinib_signature_by_top_dose.png", width = 1.25, height = 1.75, dpi = 900)

#### Drug UMAPs based on correlation coefficient of DEG beta coefficients

deg_list <- diff_test_result %>% 
  filter(grepl(pattern = "log", term), 
         q_value < 0.001) %>%
  pull(id) %>%
  unique()

deg_precorr <- list()
deg_corr <- list()
for (cell_type in c("A172", "T98G", "U87MG", "BT333")){
  
  deg_precorr[[cell_type]] <- diff_test_result %>%
    filter(Cell.Line == cell_type) %>%
    filter(grepl("log", term),id %in% deg_list) %>%
    # mutate(condition = paste0(Cell.Line,"_",treatment)) %>%
    dplyr::select(condition = treatment, id, normalized_effect) %>%
    tidyr::spread(key = condition, value = normalized_effect) %>%
    as.data.frame()
  
  row.names(deg_precorr[[cell_type]]) <- deg_precorr[[cell_type]]$id
  deg_precorr[[cell_type]]$id <- NULL
  
  deg_corr[[cell_type]] <- cor(deg_precorr[[cell_type]], method = "pearson")
}

deg_corr_df <- do.call("cbind", deg_corr)
row.names(deg_corr_df) <- sapply(row.names(deg_corr_df),function(x){stringr::str_split(x, pattern = "_")[[1]][2]})
# saveRDS(deg_corr_df, "correlation_coefficient_matrix.RDS")

deg_corr_df <- readRDS("correlation_coefficient_matrix.RDS")
gene_meta <- data.frame(id = c("Afatinib", "Brigatinib", "CUDC-101", "EAI045", "Neratinib", "Osimertinib"), gene_short_name = NA)
row.names(gene_meta) <- gene_meta$id
cell_meta <- data.frame(id = colnames(deg_corr_df), id2 = colnames(deg_corr_df)) %>%
  separate(id2, into = c("Cell.Line", "treatment"), sep = "_")
row.names(cell_meta) <- cell_meta$id

cds_corr <- new_cell_data_set(deg_corr_df,
                              gene_metadata = gene_meta,
                              cell_metadata = cell_meta)

reducedDims(cds_corr)[["PCA"]] <- t(exprs(cds_corr))
cds_corr <- align_cds(cds_corr, 
                      alignment_group = "Cell.Line", 
                      residual_model_formula_str = "~Cell.Line")

cds_corr <- reduce_dimension(cds_corr, umap.min_dist = 0.1, umap.n_neighbors = 10)
colData(cds_corr)$UMAP1 <- reducedDims(cds_corr)[["UMAP"]][,1]
colData(cds_corr)$UMAP2 <- reducedDims(cds_corr)[["UMAP"]][,2]
# saveRDS("correlation_coefficient_cds.rds")

cor_cds <- readRDS("correlation_coefficient_cds.rds")

ggplot(colData(cor_cds) %>% as.data.frame(), aes(x = UMAP1, y = UMAP2, fill = treatment, shape = Cell.Line), label_cell_groups = FALSE) +
  geom_point(color = "black", size = 3, stroke = 0.3) +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "right",
        text = element_text(size = 6),
        legend.title = element_text(size = 9),
        legend.key.width = unit(0.1,"line"),
        legend.key.height = unit(0.1,"line")) +
  scale_fill_manual("EGFRi",
                     values = c("Afatinib" = "firebrick2","Brigatinib" = "navy",
                                "CUDC-101" = "forestgreen","EAI045" = "deepskyblue3",
                                "Neratinib" = "brown4","Osimertinib" = "darkorange3",
                                "DMSO" = "grey80")) +
  scale_shape_manual("GBM line",
                     values = c("A172" = 21,
                                "T98G" = 22,
                                "U87MG" = 23,
                                "BT333" = 24)) +
  guides(fill = guide_legend(override.aes = list(size=2))) +
  guides(shape = guide_legend(override.aes = list(size=2)))
ggsave("Drug_UMAPs/Drug_UMAP_by_treatment.png", width = 1.75, height = 1.25, dpi = 900)

ggplot(colData(cor_cds) %>% as.data.frame(), aes(x = UMAP1, y = UMAP2, fill = treatment), label_cell_groups = FALSE) +
  geom_point(color = "black", size = 3, stroke = 0, shape = 21) +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "right",
        text = element_text(size = 6),
        legend.title = element_text(size = 9),
        legend.key.width = unit(0.1,"line"),
        legend.key.height = unit(0.1,"line")) +
  scale_fill_manual("EGFRi",
                    values = c("Afatinib" = "firebrick2","Brigatinib" = "navy",
                               "CUDC-101" = "forestgreen","EAI045" = "deepskyblue3",
                               "Neratinib" = "brown4","Osimertinib" = "darkorange3",
                               "DMSO" = "grey80")) +
  guides(fill = guide_legend(override.aes = list(size=2))) 
ggsave("Drug_UMAPs/Drug_UMAP_by_treatment_for_legend.png", width = 2, height = 2, dpi = 900)

proliferation_index_summary <- colData(cds) %>% 
  as.data.frame() %>%
  group_by(Cell.Line, treatment, dose, batch) %>%
  summarize(mean_proliferation_index = mean(proliferation_index)) %>%
  group_by(Cell.Line, treatment, dose) %>%
  summarize(mean_proliferation_index = mean(mean_proliferation_index)) %>%
  group_by(Cell.Line) %>%
  mutate(scaled_mean_proliferation_index = mean_proliferation_index/mean_proliferation_index[treatment == "DMSO"],
         sample_id = paste0(Cell.Line,"_",treatment)) %>%
  ungroup() %>%
  filter(dose == 10000)

colData(cor_cds)$proliferation_index <- sapply(colData(cor_cds)$id, 
                                              function(x){proliferation_index_summary %>% filter(sample_id == x) %>% pull(scaled_mean_proliferation_index)})

ggplot(colData(cor_cds) %>% as.data.frame(), aes(x = UMAP1, y = UMAP2, fill = proliferation_index, shape = Cell.Line), label_cell_groups = FALSE) +
  geom_point(color = "black", size = 3, stroke = 0.3) +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "right",
        text = element_text(size = 6),
        legend.title = element_text(size = 9),
        legend.key.width = unit(0.5,"line"),
        legend.key.height = unit(0.3,"line")) +
  viridis::scale_fill_viridis("PI",option = "viridis") +
  scale_shape_manual("GBM line",
                     values = c("A172" = 21,
                                "T98G" = 22,
                                "U87MG" = 23,
                                "BT333" = 24)) +
  guides(shape = guide_legend(override.aes = list(size=2)))
ggsave("Drug_UMAPs/Drug_UMAP_by_proliferaiton_index.png", width = 1.75, height = 1.25, dpi = 900)

colData(cds)$KRAS <- calculate_aggregate_expression_score(cds, signature_genes = hallmarksGSC$gsc$HALLMARK_KRAS_SIGNALING_UP)

KRAS_summary <- colData(cds) %>% 
  as.data.frame() %>%
  group_by(Cell.Line) %>%
  mutate(KRAS = scale(KRAS)) %>%
  group_by(Cell.Line, treatment, dose, batch) %>%
  summarize(mean_KRAS = mean(KRAS)) %>%
  group_by(Cell.Line, treatment, dose) %>%
  summarize(mean_KRAS = mean(mean_KRAS)) %>%
  group_by(Cell.Line) %>%
  mutate(scaled_mean_proliferation_index = mean_KRAS/mean_KRAS[treatment == "DMSO"],
         sample_id = paste0(Cell.Line,"_",treatment)) %>%
  ungroup() %>%
  filter(dose == 10000)

colData(cor_cds)$KRAS <- sapply(colData(cor_cds)$id, 
                                function(x){KRAS_summary %>% filter(sample_id == x) %>% pull(mean_KRAS)})

ggplot(colData(cor_cds) %>% as.data.frame(), aes(x = UMAP1, y = UMAP2, fill = KRAS, shape = Cell.Line), label_cell_groups = FALSE) +
  geom_point(color = "black", size = 3, stroke = 0.3) +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "right",
        text = element_text(size = 6),
        legend.title = element_text(size = 9),
        legend.key.width = unit(0.5,"line"),
        legend.key.height = unit(0.3,"line")) +
  viridis::scale_fill_viridis("KRAS\nscore",option = "plasma", direction = -1) +
  scale_shape_manual("GBM line",
                     values = c("A172" = 21,
                                "T98G" = 22,
                                "U87MG" = 23,
                                "BT333" = 24)) +
  guides(shape = guide_legend(override.aes = list(size=2)))
ggsave("Drug_UMAPs/Drug_UMAP_by_KRAS.png", width = 1.75, height = 1.25, dpi = 900)
####################################################################################################
