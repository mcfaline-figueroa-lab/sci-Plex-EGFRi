library(tidyverse)
library(monocle3)
library(ComplexHeatmap)

this.dir <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(this.dir))
rm(this.dir)

source("calculate_aggreg_expression.R")

cds <- readRDS("sci-Plex-EGFRi-combination-PI3Ki-cds.RDS")

# ================================================================================
# Calculate adaptive resistance signature
# ================================================================================
adaptive_up_genes <- read_csv("adaptive_program.csv",
                              col_names = TRUE) %>%
  filter(grepl(pattern = "pregulated", gene_module))

colData(cds)$adaptive_up <- calculate_aggregate_expression_score(cds, 
                                                                 signature_genes = adaptive_up_genes$gene_short_name,
                                                                 from_id = FALSE)

# ====================================================================================
# Plotting adaptive up signature
# ====================================================================================

# ridgeplot showing single-cell distribution of signature
adaptive_plot <- colData(cds) %>% as.data.frame() %>%
  group_by(cell_line) %>%
  mutate(dose1 = factor(dose1, levels = c("0", "10", "100", "500", "1000", "5000", "10000"))) %>%
  mutate(adaptive_up=scale(adaptive_up) %>% as.vector()) %>%
  mutate(mean_DMSO = mean(adaptive_up[drug1 == "DMSO" & drug2 == "DMSO"]))

adaptive_plot_DMSO <- adaptive_plot %>%
  filter(drug1 == "DMSO" & dose2 == "0")

adaptive_plot_MTX <- adaptive_plot %>%
  filter(drug1 == "MTX-211")

ggplot(adaptive_plot %>% 
         filter(drug1 %in% c("Osimertinib") & dose1 %in% c(0, 10, 100, 1000, 10000)) %>%
         bind_rows(adaptive_plot_MTX %>% mutate(dose2 = "MTX-211") 
                   %>% filter(dose1 %in% c(0, 10, 100, 1000, 10000))) %>%
         bind_rows(adaptive_plot_DMSO %>% mutate(dose2 = "MTX-211")) %>%
         bind_rows(adaptive_plot_DMSO) %>%
         bind_rows(adaptive_plot_DMSO %>% mutate(dose2 = "10")) %>%
         bind_rows(adaptive_plot_DMSO %>% mutate(dose2 = "100")) %>%
         bind_rows(adaptive_plot_DMSO %>% mutate(dose2 = "1000")) %>%
         bind_rows(adaptive_plot_DMSO %>% mutate(dose2 = "10000")) %>%
         filter(cell_line == "BT228" & dose2 %in% c("MTX-211","0", "100", "10000")) %>%
         filter(dose1 == 0 | dose1 == 10000) %>%
         mutate(dose2 = factor(dose2, levels = c("MTX-211","0", "100", "10000"))) %>%
         mutate(dose1_fill = case_when(
           drug1 == "DMSO" & drug2 == "DMSO" ~ "black",
           drug1 == "MTX-211" ~ "red",
           drug1 == "Osimertinib" & dose2 == "0" ~ "blue",
           drug1 == "Osimertinib" & dose2 == "100" ~ "yellow1",
           drug1 == "Osimertinib" & dose2 == "10000" ~ "yellow2"
         )) %>%
         mutate(adaptive_up = case_when(
           adaptive_up < -3 ~ -3,
           TRUE ~ adaptive_up
         )), 
       aes(x = adaptive_up,
           y = 1,
           fill = dose1_fill,
           alpha = dose1_fill,
           group = dose1_fill)) +
  facet_grid(~dose2) +
  ggridges::geom_density_ridges(panel_scaling = T, 
                                size = 0.2, 
                                show.legend = F,
                                aes(height = after_stat(ndensity))) +
  scale_alpha_manual(values = c(0.4, 0.5, 0.5, 0.4, 0.75),
                     breaks = c("black", "red", "blue", "yellow1", "yellow2"),
                     labels = c("DMSO", "MTX-211", "0", "100", "10000")) +
  scale_fill_manual(values = c("black", "red", "blue", "goldenrod1", "goldenrod1"),
                    breaks = c("black", "red", "blue", "yellow1", "yellow2"),
                    labels = c("DMSO", "MTX-211", "0", "100", "10000")) +
  scale_x_continuous(limits = c(-3, 3.5)) +
  xlab("Scaled Adaptive Resistance Score") +
  ylab("Density") +
  theme(axis.ticks = element_line(linewidth = 0.3),
        axis.ticks.length = unit(0.5, "mm"),
        text = element_text(size = 5),
        axis.line = element_line(linewidth = 0.01),
        legend.position = "bottom",
        legend.margin = margin(-4, 0, 0, 0),
        legend.key.size = unit(2, "mm"),
        legend.key.width = unit(6, "mm"),
        legend.key.height = unit(2, "mm"),
        legend.spacing.x = unit(3, "mm"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  monocle3:::monocle_theme_opts()
ggsave("ridge_plot_BT228_adaptive_lowandhigh_filterpaxdose_includeMTX_fixcolors_scaled.png", 
       dpi = 600, width = 2.25, height = 1.15)

# line plot over dose
plot_data <- adaptive_plot %>% 
  filter(drug1 %in% c("Osimertinib", "DMSO", "MTX-211", "Paxilisib")) %>%
  filter(dose2 == "0" | dose1 == dose2) %>%
  mutate(plot_fill = case_when(
    dose2 != "0" ~ "Osi + Pax",
    TRUE ~ drug1
  )) %>%
  mutate(plot_fill = factor(plot_fill, levels = c("DMSO", "MTX-211", "Osimertinib", "Paxilisib", "Osi + Pax"))) %>%
  group_by(cell_line,plot_fill, dose1) %>%
  summarise(adaptive_up = mean(adaptive_up)) %>%
  mutate(dose1 = as.double(as.character(dose1))) %>%
  mutate(dose1 = log10(dose1 + 1))

ggplot(adaptive_plot %>% 
         filter(cell_line == "BT228") %>%
         filter(drug1 %in% c("Osimertinib", "DMSO", "MTX-211", "Paxilisib")) %>%
         filter(dose2 == "0" | dose1 == dose2) %>%
         mutate(plot_fill = case_when(
           dose2 != "0" ~ "Osi + Pax",
           TRUE ~ drug1
         )) %>%
         mutate(dose1 = as.double(as.character(dose1))) %>%
         mutate(dose1 = log10(dose1 + 1)) %>%
         mutate(plot_fill = factor(plot_fill, levels = c("DMSO", "MTX-211", "Osimertinib", "Paxilisib", "Osi + Pax"))),
       aes(x = dose1, y = adaptive_up, color = plot_fill, group = plot_fill)) +
  # facet_wrap(~cell_line, scales = "free_y") +
  # geom_smooth(se = F, method = "loess", formula = "y ~ x") +
  geom_point(data = plot_data %>% filter(cell_line == "BT228" & plot_fill != "DMSO"), aes(x = dose1, y = adaptive_up)) +
  # geom_line() +
  # scale_fill_manual(values = RColorBrewer::brewer.pal(n = 5, "Set2")) +
  # scale_color_manual(values = RColorBrewer::brewer.pal(n = 5, "Set2")) +
  # guides(color = "none",
  #        fill = guide_legend(title = "Drug", override.aes = list(shape = 21, size = 3, alpha = 1))) +
  ylab("Adaptive Aggregate Score") +
  # scale_y_continuous(limits = c(-3, 3)) +
  theme(text = element_text(size = 7),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  monocle3:::monocle_theme_opts()

test <- adaptive_plot %>%
  group_by(cell_line, drug1, dose1, drug2, dose2) %>%
  summarise(count = n())


colData(cds)$MKI67 <- calculate_aggregate_expression_score(cds, 
                                                           signature_genes = "ERBB3",
                                                           from_id = FALSE)

mean_proliferation <- colData(cds) %>% as.data.frame() %>%
  group_by(cell_line, drug1, dose1, drug2, dose2, hash_plate) %>%
  dplyr::summarise(mean_prolif = mean(MKI67)) %>%
  group_by(cell_line, drug1, dose1, drug2, dose2) %>%
  dplyr::summarise(mean_prolif = mean(mean_prolif), .groups = "drop") %>%
  # group_by(cell_line) %>%
  # mutate(mean_prolif_norm = mean_prolif/mean_prolif[drug1 == "DMSO" & drug2 == "DMSO"]) %>%
  ungroup() %>%
  filter(drug1 %in% c("Osimertinib", "DMSO", "MTX-211", "Paxilisib")) %>%
  filter(dose2 == "0" | dose1 == dose2) %>%
  mutate(plot_fill = case_when(
    dose2 != "0" ~ "Osi + Pax",
    TRUE ~ drug1
  )) %>%
  mutate(dose1 = factor(dose1, levels = c("0", "10", "100", "500", "1000", "5000", "10000"))) %>%
  mutate(plot_fill = case_when(
    dose2 != "0" ~ "Osi + Pax",
    TRUE ~ drug1
  )) %>%
  mutate(plot_fill = factor(plot_fill, levels = c("DMSO", "MTX-211", "Osimertinib", "Paxilisib", "Osi + Pax"))) %>%
  mutate(dose1 = as.double(as.character(dose1))) %>%
  mutate(dose1 = log10(dose1 + 1))

prolif_controls <- mean_proliferation %>%
  filter(drug1 == "DMSO" & drug2 == "DMSO") %>%
  dplyr::select(cell_line, mean_prolif)

mean_proliferation_plot <- colData(cds) %>% as.data.frame() %>%
  filter(drug1 %in% c("Osimertinib", "MTX-211", "Paxilisib")) %>%
  filter(dose2 == "0" | dose1 == dose2) %>%
  mutate(plot_fill = case_when(
    dose2 != "0" ~ "Osi + Pax",
    TRUE ~ drug1
  )) %>%
  mutate(plot_fill = factor(plot_fill, levels = c("DMSO", "MTX-211", "Osimertinib", "Paxilisib", "Osi + Pax"))) %>%
  # left_join(prolif_controls, by = c("cell_line")) %>%
  # mutate(plot_MKI67 = proliferation_index/mean_prolif) %>%
  mutate(dose1 = as.double(as.character(dose1))) %>%
  mutate(dose1 = log10(dose1 + 1))
  
ggplot(mean_proliferation_plot,
       aes(x = dose1, y = mean_prolif, color = plot_fill, group = plot_fill)) +
  facet_wrap(~cell_line, scales = "free_y") +
  geom_smooth(se = F, method = "loess",
              formula = "y ~ log(x)") +
  # geom_smooth(se = F, method = "nls", 
  #             formula = "y ~ min + ((max - min) / (1 + exp(hill_coefficient * (ec50 - x))))",
  #             method.args = list(start = list(min = 1.67, max = 397, ec50 = -7, hill_coefficient = 1)) ) +
  # geom_point(data = mean_proliferation %>% filter(drug1 != "DMSO"), aes(x = dose1, y = mean_prolif_norm)) +
  # geom_line() +
  # scale_fill_manual(values = RColorBrewer::brewer.pal(n = 5, "Set2")) +
  # scale_color_manual(values = RColorBrewer::brewer.pal(n = 5, "Set2")) +
  # guides(color = "none",
  #        fill = guide_legend(title = "Drug", override.aes = list(shape = 21, size = 3, alpha = 1))) +
  ylab("Adaptive Aggregate Score") +
  # scale_y_continuous(limits = c(-3, 3)) +
  theme(text = element_text(size = 7),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  monocle3:::monocle_theme_opts()


# linear regression over dose
adaptive_plot <- colData(cds) %>% as.data.frame() %>%
  dplyr::group_by(cell_line) %>%
  dplyr::mutate(adaptive_up = scale(adaptive_up)) %>%
  dplyr::mutate(dose1 = factor(dose1, levels = c("0", "10", "100", "500", "1000", "5000", "10000"))) %>%
  filter(drug1 %in% c("Osimertinib", "DMSO", "MTX-211", "Paxilisib")) %>%
  filter(dose2 == "0" | dose2 == "100" | dose2 == "1000" | dose2 == "10000") %>%
  mutate(plot_fill = case_when(
    dose2 == "100" ~ "Osi + Pax",
    dose2 == "1000" ~ "Osi + Pax",
    dose2 == "10000" ~ "Osi + Pax",
    TRUE ~ drug1
  )) %>%
  mutate(plot_fill = factor(plot_fill, levels = c("DMSO", "MTX-211", "Osimertinib", "Paxilisib", "Osi + Pax"))) %>%
  mutate(test_facet = case_when(
    dose2 == "100" ~ "0.1uM",
    dose2 == "1000" ~ "1uM",
    dose2 == "10000" ~ "10uM",
    TRUE ~ "0"
  )) %>%
  mutate(dose1 = as.double(as.character(dose1)))

adaptive_plot_new <- adaptive_plot %>%
  filter(plot_fill %in% c("Osi + Pax")) %>%
  bind_rows(adaptive_plot %>%
              filter(plot_fill %in% c("Osimertinib", "MTX-211")) %>%
              mutate(test_facet = "0.1uM")) %>%
  bind_rows(adaptive_plot %>%
              filter(plot_fill %in% c("Osimertinib", "MTX-211")) %>%
              mutate(test_facet = "1uM")) %>%
  bind_rows(adaptive_plot %>%
              filter(plot_fill %in% c("Osimertinib", "MTX-211")) %>%
              mutate(test_facet = "10uM"))

ggplot(adaptive_plot_new %>%
         filter(cell_line == "BT228") %>%
         group_by(cell_line) %>%
         dplyr::mutate(test_facet = factor(test_facet, levels = c("0.1uM", "1uM", "10uM"))),
       aes(x = dose1, y = adaptive_up, color = plot_fill, group = plot_fill)) +
  geom_smooth(aes(linetype = plot_fill),
              se = F, method = "glm",
              formula = "y ~ log10(x)") +
  scale_color_manual("Exposure(s)",
                     values = c(RColorBrewer::brewer.pal(n = 3, "Set1")[c(1,2)], "goldenrod1")) +
  scale_linetype_manual("Exposure(s)",
                        values = c(4,4,1)) +
  scale_x_log10() +
  ylab("Scaled Adaptive Score") +
  xlab("Log10 Dose") +
  guides(color = guide_legend(byrow = T, direction = "horizontal", override.aes = list(size = 2)),
         shape = guide_legend(byrow = T, direction = "horizontal", override.aes = list(size = 2))) +
  theme(text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 10),
        legend.key.height = unit(4, "mm"),
        legend.margin = margin(t = -5),
        legend.position = "bottom") +
  monocle3:::monocle_theme_opts()

ggsave("adaptive_score_glm_line_dotted_BT228.png",
       dpi = 600, height = 2, width = 5.25)


ggplot(adaptive_plot_new %>%
         filter(cell_line == "BT333") %>%
         group_by(cell_line) %>%
         dplyr::mutate(test_facet = factor(test_facet, levels = c("0.1uM", "1uM", "10uM"))),
       aes(x = dose1, y = adaptive_up, color = plot_fill, group = plot_fill)) +
  geom_smooth(aes(linetype = plot_fill),
              se = F, method = "glm",
              formula = "y ~ log10(x)") +
  scale_color_manual("Exposure(s)",
                     values = c(RColorBrewer::brewer.pal(n = 3, "Set1")[c(1,2)], "goldenrod1")) +
  scale_linetype_manual("Exposure(s)",
                        values = c(4,4,1)) +
  scale_x_log10() +
  ylab("Scaled Adaptive Score") +
  xlab("Log10 Dose") +
  guides(color = guide_legend(byrow = T, direction = "horizontal", override.aes = list(size = 2)),
         shape = guide_legend(byrow = T, direction = "horizontal", override.aes = list(size = 2))) +
  theme(text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 10),
        legend.key.height = unit(4, "mm"),
        legend.margin = margin(t = -5),
        legend.position = "bottom") +
  monocle3:::monocle_theme_opts()

ggsave("adaptive_score_glm_line_dotted_BT333.png",
       dpi = 600, height = 2, width = 5.25)

