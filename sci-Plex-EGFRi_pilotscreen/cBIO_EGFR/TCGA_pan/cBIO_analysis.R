library(tidyverse)
library(survival)
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
library(condsurv)

this.dir <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(this.dir))
rm(this.dir)

dir.create("figures")

clin_data <- read_tsv(file = "combined_study_clinical_data_withgrade.tsv",
                      col_names = TRUE)
survival <- read_tsv("KM_Plot__Overall_Survival__(months).txt", col_names = TRUE)
EGFR_cn <- read_tsv("Log2 copy-number values_EGFR.txt", col_names = TRUE) %>%
  drop_na() %>%
  arrange(desc(EGFR)) %>%
  mutate(rank = row_number()) %>%
  mutate(EGFR_cn = case_when(
    EGFR >= 2 ~ "Amp",
    TRUE ~ "WT"
  ))

EGFR_mut <- read_tsv("mutations_EGFR.txt", col_names = TRUE)

sig_genes_mRNA <- read_tsv(file = "mRNA Expression, RSEM (Batch normalized from Illumina HiSeq_RNASeqV2)_UPSIGNATURE.txt")
sig_genes_mRNA_agg <- sig_genes_mRNA %>%  
  pivot_longer(cols = 3:ncol(sig_genes_mRNA), names_to = "gene", values_to = "count") %>%
  group_by(STUDY_ID, SAMPLE_ID) %>%
  summarise(agg_score = sum(count))

all_data <- left_join(sig_genes_mRNA_agg, 
                           clin_data, 
                           by = c("STUDY_ID" = "Cancer Study",
                                  "SAMPLE_ID" = "Sample ID")) %>%
  left_join(survival, by = c("STUDY_ID" = "Study ID",
                             "Patient ID" = "Patient ID")) %>%
  left_join(EGFR_cn, by = c("STUDY_ID" = "STUDY_ID",
                            "SAMPLE_ID" = "SAMPLE_ID")) %>%
  left_join(EGFR_mut, by = c("STUDY_ID" = "STUDY_ID",
                             "SAMPLE_ID" = "SAMPLE_ID")) %>%
  mutate(OS_STATUS = case_when(
    OS_STATUS == "1:DECEASED" ~ 1,
    OS_STATUS == "0:LIVING" ~ 0
  ))

factor_order <- all_data %>% 
  filter(!is.na(agg_score)) %>% 
  group_by(`Cancer Type`) %>% 
  summarise(count= n(), med = median(agg_score)) %>%
  arrange(desc(med))

ggplot(all_data %>% filter(`Cancer Type` %in% c("Glioblastoma", "Glioma", 
                                                "Breast Cancer", "Leukemia",
                                                "Non-Small Cell Lung Cancer",
                                                     "Colorectal Cancer",
                                                     "Melanoma")) %>%
         filter(!is.na(agg_score)) %>%
         mutate(test = factor(`Cancer Type`, levels = factor_order$`Cancer Type`)), 
       aes(x = test, y = log10(agg_score))) +
  geom_boxplot(aes(color = `Cancer Type`), show.legend = FALSE, size = 0.25, outlier.size = 0.5) +
  viridis::scale_color_viridis(discrete = TRUE, option = "plasma") +
  ylab("Log10 Aggregate Score") +
  monocle3:::monocle_theme_opts() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        text = element_text(size = 6),
        axis.ticks = element_line(linewidth = unit(0.1, "mm")))
ggsave("figures/pan_cancer_EGFRi_log_agg_score_comparison_ranked.png",
       dpi = 600, height = 2, width = 2.5)

gbm_TCGA <- all_data %>%
  filter(`Cancer Type` %in% c("Glioblastoma", "Glioma")) %>%
  mutate(Mutation_Group = case_when(
    EGFR_cn == "WT" & EGFR.y == "WT" ~ "WT",
    TRUE ~ "Mut"
  )) %>%
  mutate(`Neoplasm Histologic Grade` = case_when(
    is.na(`Neoplasm Histologic Grade`) == TRUE ~ "G4",
    TRUE ~ `Neoplasm Histologic Grade`
  ))

factor_order <- gbm_TCGA %>% 
  filter(!is.na(agg_score)) %>% 
  group_by(`Neoplasm Histologic Grade`) %>% 
  summarise(count= n(), med = median(agg_score)) %>%
  arrange(desc(med))

ggplot(gbm_TCGA %>% 
         filter(!is.na(agg_score)) %>%
         mutate(test = factor(`Neoplasm Histologic Grade`, levels = factor_order$`Neoplasm Histologic Grade`)), 
       aes(x = test, y = log10(agg_score))) +
  geom_boxplot(aes(color = `Neoplasm Histologic Grade`), show.legend = FALSE, size = 0.25, outlier.size = 0.5) +
  viridis::scale_color_viridis(discrete = TRUE, option = "D", direction = -1) +
  ylab("Log10 Aggregate Score") +
  xlab("Glioma Histologic Grade") +
  monocle3:::monocle_theme_opts() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        # axis.title.x = element_blank(),
        text = element_text(size = 7),
        axis.ticks = element_line(linewidth = unit(0.1, "mm")))
ggsave("figures/glioma_grade_EGFRi_log_agg_score_comparison_ranked.png",
       dpi = 600, height = 1.5, width = 1.5)


counts_half_facet <- gbm_TCGA %>% 
  group_by(`Neoplasm Histologic Grade`) %>%
  filter(!is.na(agg_score)) %>%
  arrange(desc(agg_score)) %>%
  mutate(rank = row_number()) %>%
  mutate(signature_group = case_when(
    rank <= 0.5*max(rank) ~ "high expression",
    TRUE ~ "low expression"
  )) %>%
  group_by(`Neoplasm Histologic Grade`, signature_group) %>%
  summarise(count = n())

gbm_TCGA_facet <- gbm_TCGA %>% 
  group_by(`Neoplasm Histologic Grade`) %>%
  filter(!is.na(agg_score)) %>%
  arrange(desc(agg_score)) %>%
  mutate(rank = row_number()) %>%
  mutate(signature_group = case_when(
    rank <= 0.5*max(rank) ~ "high expression",
    TRUE ~ "low expression"
  ))

for (g in c("G2", "G3", "G4")) {
  survfit2(Surv(OS_MONTHS, OS_STATUS) ~ signature_group, data = gbm_TCGA_facet %>% 
             filter(`Neoplasm Histologic Grade` == g)) %>% 
    ggsurvfit(size = 0.4) +
    ggtitle(paste("Grade", g)) +
    labs(
      x = "Months",
      y = "Overall survival\nprobability"
    ) +
    # add_risktable() +
    add_quantile(y_value = 0.6, color = "gray50", linewidth = 0.25) +
    add_pvalue(location = "annotation", size = 2, rho = 0) +
    scale_color_manual(labels = c(paste0("High Signature Expression (n = ", counts_half_facet %>%
                                           filter(`Neoplasm Histologic Grade` == g & signature_group == "high expression") %>% 
                                           pull(count),
                                         ")"), 
                                  paste0("Low Signature Expression (n = ", counts_half_facet %>%
                                           filter(`Neoplasm Histologic Grade` == g & signature_group == "low expression") %>% 
                                           pull(count),
                                         ")")),
                       breaks = c("high expression", "low expression"),
                       values = viridis::plasma(n = 5, direction = -1)[c(2,4)]) +
    scale_y_continuous(expand = c(0,0)) +
    theme(text = element_text(size = 6),
          plot.title = element_text(hjust = 0.5),
          axis.text = element_text(size = 6),
          axis.title = element_text(size = 6),
          axis.ticks = element_line(linewidth = unit(0.1, "mm")),
          legend.key.height = unit(0.1, "cm"),
          legend.key.width = unit(0.3, "cm"),
          legend.box.margin = unit(c(t = -3,0,0,-10), "mm")) +
    guides(color = guide_legend(byrow = TRUE, nrow = 2, label.theme = element_text(size = 6))) +
    monocle3:::monocle_theme_opts()
  ggsave(paste0("figures/km_curve_top50_vs_bottom50_",g ,"_logrank_long.png"),
         dpi = 600, height = 1.8, width = 2)
}


survfit2(Surv(OS_MONTHS, OS_STATUS) ~ signature_group, data = gbm_TCGA %>% 
           filter(!is.na(agg_score)) %>%
           filter(`Cancer Type` == "Glioblastoma") %>%
           arrange(desc(agg_score)) %>%
           mutate(rank = row_number()) %>%
           mutate(signature_group = case_when(
             rank <= 0.5*max(rank) ~ "high expression",
             TRUE ~ "low expression"
           ))) %>% 
  ggsurvfit() +
  labs(
    x = "Months",
    y = "Overall survival probability"
  ) +
  # add_risktable() +
  add_quantile(y_value = 0.6, color = "gray50", linewidth = 0.75) +
  add_pvalue(location = "annotation", rho = 1) +
  # scale_color_manual(labels = c("High Signature\nn=73", "Low Signature\nn=73"), 
  #                    breaks = c("High Signature", "Low Signature"), 
  #                    values = c("red", "grey60")) +
  scale_y_continuous(expand = c(0,0)) +
  theme(text = element_text(size = 6)) +
  monocle3:::monocle_theme_opts()
ggsave("figures/km_curve_top50_vs_bottom50_gbm_alone_gwpvalue.png",
       dpi = 600, height = 3, width = 4)


counts_half <- gbm_TCGA %>% 
  filter(!is.na(agg_score)) %>%
  arrange(desc(agg_score)) %>%
  mutate(rank = row_number()) %>%
  mutate(signature_group = case_when(
    rank <= 0.5*max(rank) ~ "high expression",
    TRUE ~ "low expression"
  )) %>%
  group_by(signature_group) %>%
  summarise(count = n())

survfit2(Surv(OS_MONTHS, OS_STATUS) ~ signature_group, data = gbm_TCGA %>% 
           filter(!is.na(agg_score)) %>%
           arrange(desc(agg_score)) %>%
           mutate(rank = row_number()) %>%
           mutate(signature_group = case_when(
             rank <= 0.5*max(rank) ~ "high expression",
             TRUE ~ "low expression"
           ))) %>% 
  ggsurvfit(size = 0.4) +
  labs(
    x = "Months",
    y = "Overall survival\nprobability"
  ) +
  # add_risktable() +
  add_quantile(y_value = 0.6, color = "gray50", linewidth = 0.5) +
  add_pvalue(location = "annotation", size = 2) +
  scale_color_manual(labels = c("High Signature Expression (n = 337)", "Low Signature Expression (n = 337)"),
                     breaks = c("high expression", "low expression"),
                     values = viridis::plasma(n = 5, direction = -1)[c(2,4)]) +
  scale_y_continuous(expand = c(0,0)) +
  theme(text = element_text(size = 6),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.ticks = element_line(linewidth = unit(0.1, "mm")),
        legend.key.height = unit(0.1, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.box.margin = unit(c(t = -3,0,0,-10), "mm")) +
  guides(color = guide_legend(byrow = TRUE, nrow = 2, label.theme = element_text(size = 6))) +
  monocle3:::monocle_theme_opts()


# =============================================================================
# random set of 48 genes for comparison
# =============================================================================

clin_data <- read_tsv(file = "combined_study_clinical_data_withgrade.tsv",
                      col_names = TRUE)
survival <- read_tsv("KM_Plot__Overall_Survival__(months).txt", col_names = TRUE)
EGFR_cn <- read_tsv("Log2 copy-number values_EGFR.txt", col_names = TRUE) %>%
  drop_na() %>%
  arrange(desc(EGFR)) %>%
  mutate(rank = row_number()) %>%
  mutate(EGFR_cn = case_when(
    EGFR >= 2 ~ "Amp",
    TRUE ~ "WT"
  ))
# ggplot(EGFR_cn, aes(x = rank, y = EGFR)) +
#   geom_point() +
#   theme(axis.ticks = element_blank(),
#         axis.text.x = element_blank()) +
#   monocle3:::monocle_theme_opts()

EGFR_mut <- read_tsv("mutations_EGFR.txt", col_names = TRUE)

sig_genes_mRNA <- read_tsv(file = "mRNA Expression, RSEM (Batch normalized from Illumina HiSeq_RNASeqV2)_random_gene_set.txt")
sig_genes_mRNA_agg <- sig_genes_mRNA %>%  
  pivot_longer(cols = 3:ncol(sig_genes_mRNA), names_to = "gene", values_to = "count") %>%
  group_by(STUDY_ID, SAMPLE_ID) %>%
  summarise(agg_score = sum(count))

all_data <- left_join(sig_genes_mRNA_agg, 
                      clin_data, 
                      by = c("STUDY_ID" = "Cancer Study",
                             "SAMPLE_ID" = "Sample ID")) %>%
  left_join(survival, by = c("STUDY_ID" = "Study ID",
                             "Patient ID" = "Patient ID")) %>%
  left_join(EGFR_cn, by = c("STUDY_ID" = "STUDY_ID",
                            "SAMPLE_ID" = "SAMPLE_ID")) %>%
  left_join(EGFR_mut, by = c("STUDY_ID" = "STUDY_ID",
                             "SAMPLE_ID" = "SAMPLE_ID")) %>%
  mutate(OS_STATUS = case_when(
    OS_STATUS == "1:DECEASED" ~ 1,
    OS_STATUS == "0:LIVING" ~ 0
  ))

factor_order <- all_data %>% 
  filter(!is.na(agg_score)) %>% 
  group_by(`Cancer Type`) %>% 
  summarise(count= n(), med = median(agg_score)) %>%
  arrange(desc(med))

ggplot(all_data %>% filter(`Cancer Type` %in% c("Glioblastoma", "Glioma", 
                                                "Breast Cancer", "Leukemia",
                                                "Non-Small Cell Lung Cancer",
                                                "Colorectal Cancer",
                                                "Melanoma")) %>%
         filter(!is.na(agg_score)) %>%
         mutate(test = factor(`Cancer Type`, levels = factor_order$`Cancer Type`)), 
       aes(x = test, y = log10(agg_score))) +
  geom_boxplot(aes(color = `Cancer Type`), show.legend = FALSE, size = 0.25, outlier.size = 0.5) +
  viridis::scale_color_viridis(discrete = TRUE, option = "plasma") +
  ylab("Aggregate Score") +
  monocle3:::monocle_theme_opts() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        text = element_text(size = 6))
ggsave("figures/random_genes/pan_cancer_EGFRi_log_agg_score_comparison_ranked.png",
       dpi = 600, height = 2, width = 2.5)

gbm_TCGA <- all_data %>%
  filter(`Cancer Type` %in% c("Glioblastoma", "Glioma")) %>%
  # filter(!is.na(agg_score)) %>%
  # filter(!is.na(EGFR_cn)) %>%
  # filter(!is.na(EGFR.y)) %>%
  # arrange(desc(agg_score)) %>%
  # mutate(rank = row_number()) %>%
  # mutate(Signature_Group = case_when(
  #   rank < 0.25*max(rank) ~ "High Signature",
  #   rank > 0.75*max(rank) ~ "Low Signature"
  # )) %>%
  mutate(Mutation_Group = case_when(
    EGFR_cn == "WT" & EGFR.y == "WT" ~ "WT",
    TRUE ~ "Mut"
  )) %>%
  mutate(`Neoplasm Histologic Grade` = case_when(
    is.na(`Neoplasm Histologic Grade`) == TRUE ~ "G4",
    TRUE ~ `Neoplasm Histologic Grade`
  ))

factor_order <- gbm_TCGA %>% 
  filter(!is.na(agg_score)) %>% 
  group_by(`Neoplasm Histologic Grade`) %>% 
  summarise(count= n(), med = median(agg_score)) %>%
  arrange(desc(med))

ggplot(gbm_TCGA %>% 
         filter(!is.na(agg_score)) %>%
         mutate(test = factor(`Neoplasm Histologic Grade`, levels = factor_order$`Neoplasm Histologic Grade`)), 
       aes(x = test, y = log10(agg_score))) +
  geom_boxplot(aes(color = `Neoplasm Histologic Grade`), show.legend = FALSE, size = 0.25, outlier.size = 0.5) +
  viridis::scale_color_viridis(discrete = TRUE, option = "D", direction = -1) +
  ylab("Aggregate Score") +
  xlab("Grade") +
  monocle3:::monocle_theme_opts() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        # axis.title.x = element_blank(),
        text = element_text(size = 6))
ggsave("figures/random_genes/glioma_grade_EGFRi_log_agg_score_comparison_ranked.png",
       dpi = 600, height = 1.5, width = 2)


ggplot(gbm_TCGA) +
  geom_histogram(aes(x = agg_score, fill = Signature_Group))

survfit2(Surv(OS_MONTHS, OS_STATUS) ~ signature_group, data = gbm_TCGA %>% 
           filter(!is.na(agg_score)) %>%
           arrange(desc(agg_score)) %>%
           mutate(rank = row_number()) %>%
           mutate(signature_group = case_when(
             rank <= 0.5*max(rank) ~ "high expression",
             TRUE ~ "low expression"
           ))) %>% 
  ggsurvfit() +
  labs(
    x = "Months",
    y = "Overall survival probability"
  ) +
  # add_risktable() +
  add_quantile(y_value = 0.6, color = "gray50", linewidth = 0.75) +
  add_pvalue(location = "annotation") +
  # scale_color_manual(labels = c("High Signature\nn=73", "Low Signature\nn=73"), 
  #                    breaks = c("High Signature", "Low Signature"), 
  #                    values = c("red", "grey60")) +
  scale_y_continuous(expand = c(0,0)) +
  theme(text = element_text(size = 6)) +
  monocle3:::monocle_theme_opts()
ggsave("figures/random_genes/km_curve_top50_vs_bottom50_allglioma.png",
       dpi = 600, height = 3, width = 4)



