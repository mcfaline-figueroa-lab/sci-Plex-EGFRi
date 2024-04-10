library(tidyverse)
library(survival)
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
library(condsurv)
library(ComplexHeatmap)
library(ggsankey)

# ==========================================================================================
# Sankey plots of patient cohort
# ==========================================================================================

clin_data <- read_tsv(file = "gbm_tcga_pub2013_clinical_data.tsv",
                      col_names = TRUE)
survival <- read_tsv("KM_Plot__Overall_(months).txt", 
                     col_names = TRUE)
EGFR_cn <- read_tsv("cna.txt", 
                    col_names = TRUE) %>%
  arrange(desc(EGFR)) %>%
  mutate(rank = row_number()) %>%
  mutate(EGFR_cn = case_when(
    EGFR == "NP" ~ NA,
    EGFR >= 2 ~ "Amp",
    TRUE ~ "WT"
  ))

EGFR_mut <- read_tsv("mutations.txt", 
                     col_names = TRUE)

sig_genes_mRNA <- read_tsv(file = "mRNA expression (RNA Seq V2 RSEM).txt")

sig_genes_mRNA_agg <- sig_genes_mRNA %>%  
  pivot_longer(cols = 3:ncol(sig_genes_mRNA), names_to = "gene", values_to = "count") %>%
  mutate(count = log10(count + 1)) %>%
  group_by(STUDY_ID, SAMPLE_ID) %>%
  summarise(agg_score = sum(count))
  # mutate(log_agg_score = log10(agg_score))

all_data <- clin_data %>%
  left_join(EGFR_cn %>% select(-EGFR), by = c("Sample ID" = "SAMPLE_ID")) %>%
  left_join(EGFR_mut %>% dplyr::rename(EGFR_mut = EGFR), by = c("Sample ID" = "SAMPLE_ID")) %>%
  left_join(sig_genes_mRNA_agg, by = c("Sample ID" = "SAMPLE_ID")) %>%
  left_join(survival, by = c("Patient ID" = "Patient ID"))

df <- all_data %>%
  mutate(EGFR_cn = case_when(
    is.na(EGFR_cn) == TRUE ~ "NONE",
    TRUE ~ "AVAIL"),
    EGFR_mut = case_when(
      EGFR_mut == "NS" ~ "NONE",
      TRUE ~ "AVAIL"
    ),
    Survival = case_when(
      is.na(OS_STATUS) == TRUE ~ "NONE",
      TRUE ~ "AVAIL"
    ),
    RNA_seq = case_when(
      is.na(agg_score) == TRUE ~ "NONE",
      TRUE ~ "AVAIL"
    )) %>%
  select(`Copy Number` = EGFR_cn, `Mutation` = EGFR_mut, Survival, `RNAseq` = RNA_seq) %>%
  mutate(`Patients in Study` = "AVAIL") %>%
  make_long(`Patients in Study`, Survival, `Copy Number`, `Mutation`, `RNAseq`)

dagg <- df%>%
  dplyr::group_by(x, node)%>%
  tally()

df2 <- merge(df, dagg, by.x = c('x','node'), by.y = c('x', 'node'), all.x = TRUE) %>% 
  mutate(color_manual = case_when(
    # x == "EGFR_mut" & node == "AVAIL" & next_x == "RNA_seq" & next_node == "AVAIL" ~ "red",
    x == "RNAseq" & node == "AVAIL" ~ "red",
    TRUE ~ "grey90"
  ))

ggplot(df2, aes(x = x, 
               next_x = next_x, 
               node = node, 
               next_node = next_node,
               fill = color_manual,
               label = paste0(node,"\nn = ", n))) +
  # geom_sankey(flow.alpha = 0.5, node.color = "black", show.legend = FALSE) +
  geom_sankey(flow.alpha = 0.5,
              show.legend = FALSE, 
              alpha = 0.8) +
  # scale_fill_manual(values = c(viridis::magma(n = 7)[5:6])) +
  scale_fill_manual(values = c("grey60", "red"))+
  geom_sankey_text(show.legend = FALSE, 
                   angle = 90, 
                   size = 2, 
                   position = position_nudge(x = 0)) +
  theme_sankey(base_size = 16) +
  # scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0, 0, 0.05, 1)) +
  xlab("cBio Data Mode") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 6))
ggsave("patient_cohort_sankey.png",
       dpi= 600, height = 2.5, width = 3)

ggplot(df2, aes(x = x, 
                next_x = next_x, 
                node = node, 
                next_node = next_node,
                label = paste0(node," n=", n),
                fill = color_manual)) +
  geom_sankey(flow.alpha = 0.5, node.color = "black", show.legend = FALSE) +
  geom_sankey_label(show.legend = FALSE) +
  theme_sankey(base_size = 16) +
  xlab("cBio Data")

sig_data <- all_data %>%
  filter(!is.na(agg_score)) %>%
  mutate(OS_STATUS = case_when(
    OS_STATUS == "1:DECEASED" ~ 1,
    OS_STATUS == "0:LIVING" ~ 0
  )) %>%
  mutate(Mutation_Group = case_when(
    EGFR_cn == "WT" & EGFR_mut == "WT" ~ "WT",
    TRUE ~ "Mut"
  )) %>%
  arrange(desc(agg_score)) %>%
  mutate(rank = row_number()) %>%
  mutate(quartiles = case_when(
    rank <= 0.25*max(rank) ~ 1,
    rank > 0.25*max(rank) & rank <= 0.5*max(rank) ~ 2,
    rank > 0.5*max(rank) & rank <= 0.75*max(rank) ~ 3,
    rank > 0.75*max(rank) ~ 4
  )) %>%
  mutate(tertiles = case_when(
    rank <= 0.33*max(rank) ~ 1,
    rank > 0.33*max(rank) & rank <= 0.66*max(rank) ~ 2,
    rank > 0.66*max(rank) ~ 3
  )) %>%
  mutate(halves = case_when(
    rank <= 0.5*max(rank) ~ 1,
    rank > 0.5*max(rank) ~ 2
  ))

quartile_y <- sig_data %>% group_by(quartiles) %>% summarise(min = min(agg_score))
halves_y <- sig_data %>% group_by(halves) %>% summarise(min = min(agg_score))

density <- density(sig_data$agg_score)
df <- data.frame("agg_score" = density$x,
                 "density" = density$y) %>%
  mutate(quartiles = case_when(
    agg_score >= quartile_y$min[1] ~ "1",
    agg_score < quartile_y$min[1] & agg_score >= quartile_y$min[2] ~ "2",
    agg_score < quartile_y$min[2] & agg_score >= quartile_y$min[3] ~ "3",
    agg_score < quartile_y$min[3] ~ "4"
  )) %>%
  mutate(halves = case_when(
    agg_score >= halves_y$min[1] ~ "1",
    agg_score < halves_y$min[1] ~ "2"
  ))


ggplot(data = df %>% mutate(halves = factor(halves, levels = c("2","1"))), aes(x = agg_score, y = density)) +
  geom_area(data = filter(df, halves == 1), aes(fill = halves), alpha = 0.3) +
  geom_area(data = filter(df, halves == 2), aes(fill = halves), alpha = 0.3) +
  geom_line(linewidth = 0.2) +
  ylab("Density") +
  xlab("EGFRi Aggregate Log10 Counts") +
  scale_y_continuous(expand = c(0, 0, 0.01, 0)) +
  # scale_x_continuous(breaks = c(-3, -2, -1, 0, 1, 2, 3, 4)) +
  scale_fill_manual(values = viridis::plasma(n = 5, direction = -1)[c(2,4)],
                    labels = c("Bottom 50%\nn = 76",
                               "Top 50%\nn = 76"),
                    breaks = c(2,1)) +
  theme(
    text = element_text(size = 6),
    legend.key.size = unit(0.5, "line"),
    legend.title = element_blank(),
    legend.box.spacing = unit(0, "line"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank(),
    # text = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank()) +
  monocle3:::monocle_theme_opts()
ggsave("patient_density_agg_expression_halves.png",
       dpi = 600, height = 1.5, width = 1.5)

counts <- sig_data %>% 
  group_by(quartiles) %>% summarise(count = n())
counts_half <- sig_data %>%
  group_by(halves) %>%
  summarise(count = n())

survfit2(Surv(OS_MONTHS, OS_STATUS) ~ group,
         data = sig_data %>% 
           filter(halves %in% c(1,2)) %>%
           mutate(group = case_when(
             halves %in% c(1) ~ "Top 50%",
             halves %in% c(2) ~ "Bottom 50%"
           ))
         ) %>%
  ggsurvfit(show.legend = T, size = 0.3) +
  labs(x = "Months",
       y = "Overall survival probability") +
  add_quantile(y_value = 0.6,
               color = "black",
               linewidth = 0.3) +
  add_pvalue(location = "annotation", rho = 0, size = 1.5, x = 40) +
  scale_color_manual(
    labels = c(paste("Bottom 50%\nn =", counts_half %>% filter(halves == 2) %>% select(count) %>% pull()),
               paste("Top 50%\nn =", counts_half %>% filter(halves == 1) %>% select(count) %>% pull())),
    breaks = c("Bottom 50%", "Top 50%"),
    values = viridis::plasma(n = 5, direction = -1)[c(2,4)]
  ) +
  scale_y_continuous(expand = c(0, 0, 0.01, 0)) +
  scale_x_continuous(expand = c(0, 0, 0, 0)) +
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 6, color = "black"),
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 6, color = "black"),
    legend.text = element_text(size = 6),
    legend.position = "bottom",
    legend.key.width = unit(1, "line"),
    legend.key.height = unit(1, "line"),
    legend.spacing.y = unit(0.2, "mm"),
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-10,40,0,0),
    axis.line = element_line(linewidth = 0.1),
    axis.ticks = element_blank()
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))
ggsave("km_halves_log_count_agg_logrankp.png",
       dpi = 600, height = 1.5, width = 1.5)

