library(tidyverse)
library(monocle3)

this.dir <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(this.dir))
rm(this.dir)

# ================================================================================================
# Calculating percent cells of total for specific inhibitors
# ================================================================================================

cds_24hr <- readRDS("sci-Plex-EGFRi-largescreen-cds.RDS")
EGFRi_24hr_col <- colData(cds_24hr) %>% as.data.frame()

cds_7d <- readRDS("sci-Plex-EGFRi-7d-largescreen-cds")
EGFRi_7d_col <- colData(cds_7d) %>% as.data.frame()

joint_colData <- rbind(EGFRi_24hr_col %>% 
                         dplyr::select(cell_ID, cell_line, replicate, drug, dose, proliferation_index) %>%
                         mutate(time_point = "24hr"),
                       EGFRi_7d_col %>% 
                         dplyr::select(cell_ID, cell_line, replicate, drug, dose, proliferation_index) %>%
                         mutate(time_point = "7d"))

proportion_summary_10uM <- joint_colData %>%
  mutate(dose = case_when(
    drug == "Panitumumab" ~ dose * 10,
    TRUE ~ dose
  )) %>%
  filter(dose == 10000) %>%
  group_by(cell_line, drug, time_point, replicate) %>%
  dplyr::summarize(n = n()) %>%
  group_by(cell_line, time_point, replicate) %>%
  dplyr::mutate(total_cells = sum(n), frequency = n/sum(n) * 100) %>%
  group_by(cell_line, drug, time_point) %>%
  dplyr::summarize(mean_frequency = mean(frequency), sd_frequency = sd(frequency))

ggplot(proportion_summary_10uM %>% filter(drug %in% c("MTX-211")), 
       aes(x = time_point,y = mean_frequency)) +
  geom_bar(stat = "identity", 
           fill = "#54F2F2",
           width = 0.75,
           show.legend = F,
           aes(alpha = time_point),
           color = "black",
           linewidth = 0.2) +
  geom_errorbar(aes(ymin = mean_frequency - sd_frequency, ymax = mean_frequency + sd_frequency),
                width = 0.2, linewidth = 0.2) +
  scale_alpha_manual(values = c(0.4, 0.8)) +
  scale_y_continuous(expand = c(0.01,0)) +
  facet_grid(~cell_line, scales = "free_y") +
  ylab("% of Timepoint-Dose Cells") +
  xlab("Drug Exposure Time") +
  theme(text = element_text(size = 6),
        axis.ticks = element_line(linewidth = 0.3),
        axis.ticks.length = unit(0.7, "mm")) +
  monocle3:::monocle_theme_opts()
ggsave("proportion_MTX-211_timepoint.png",
       dpi = 600, height = 1.2, width = 1.4)

ggplot(proportion_summary_10uM %>% filter(drug %in% c("Osimertinib")), 
       aes(x = time_point,y = mean_frequency)) +
  geom_bar(stat = "identity", 
           fill = "#F76F8E",
           width = 0.75,
           show.legend = F,
           aes(alpha = time_point),
           color = "black",
           linewidth = 0.2) +
  geom_errorbar(aes(ymin = mean_frequency - sd_frequency, ymax = mean_frequency + sd_frequency),
                width = 0.2, linewidth = 0.2) +
  scale_alpha_manual(values = c(0.4, 0.8)) +
  scale_y_continuous(expand = c(0.01,0)) +
  facet_grid(~cell_line, scales = "free_y") +
  ylab("% of Timepoint-Dose Cells") +
  xlab("Drug Exposure Time") +
  theme(text = element_text(size = 6),
        axis.ticks = element_line(linewidth = 0.3),
        axis.ticks.length = unit(0.7, "mm")) +
  monocle3:::monocle_theme_opts()
ggsave("proportion_osimertinib_timepoint.png",
       dpi = 600, height = 1.2, width = 1.4)


proportion_delta_10uM <- proportion_summary_10uM %>%
  mutate(mean_frequency = mean_frequency * 100) %>%
  select(-sd_frequency) %>%
  pivot_wider(id_cols = c("cell_line", "drug"), 
              names_from = "time_point", 
              values_from = "mean_frequency") %>%
  mutate(`7d` = case_when(
    is.na(`7d`) ~ 0,
    TRUE ~ `7d`
  )) %>%
  mutate(delta = `7d` - `24hr`,
         percent_change = (`7d` - `24hr`)/`24hr`) 

saveRDS(proportion_delta_10uM, "sci-Plex-EGFRi-7d-largescreen-delta-proportion-10uM.RDS")

