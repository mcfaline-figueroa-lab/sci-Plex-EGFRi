library(tidyverse)
library(monocle3)
this.dir <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(this.dir))
rm(this.dir)

temp_df <- read_csv("EGFRi_PI3Kicell_dists_attention_singlecell.csv") %>%
  filter(target_drug == "DMSO_0_MTX-211_10000")

temp_df_plot <- temp_df %>% 
  filter(sample_drug %in% c("DMSO_0_MTX-211_100",
                            "DMSO_0_Osimertinib_10",
                            "DMSO_0_Osimertinib_100",
                            "DMSO_0_Osimertinib_1000",
                            "Paxilisib_10_Osimertinib_10",
                            "Paxilisib_1000_Osimertinib_1000",
                            "Paxilisib_100_Osimertinib_100")) %>%
  separate(sample_drug,
           into = c("drug1", "dose1", "drug2", "dose2"),
           sep = "_",
           remove = F) %>%
  mutate(plot_fill = case_when(
    drug1 == "DMSO" & drug2 == "MTX-211" ~ "red",
    drug1 == "DMSO" & drug2 == "Osimertinib" ~ "blue",
    drug1 == "Paxilisib" & drug2 == "Osimertinib" ~ "yellow",
  )) %>%
  mutate(plot_alpha = case_when(
    dose2 == "10" ~ "light",
    dose2 == "100" ~ "medium",
    TRUE ~ "full"
  ))

temp_df_plot_mean <- temp_df_plot %>%
  group_by(sample_drug) %>%
  summarise(mean_distance = mean(`Counterfactual Distance`),
            median_distance = median(`Counterfactual Distance`))

ggplot(data = temp_df_plot, aes(x = sample_drug, y = `Counterfactual Distance`)) +
  geom_violin(aes(fill = plot_fill,
                  alpha = plot_alpha),
              linewidth = 0.25) +
  geom_point(data = temp_df_plot_mean,
             aes(x = sample_drug, y = mean_distance),
             size = 0.25) +
  geom_hline(yintercept = temp_df_plot_mean %>% filter(sample_drug == "DMSO_0_MTX-211_100") %>% pull(mean_distance),
             linetype = 2,
             linewidth = 0.25) +
  scale_alpha_manual(values = c(0.2, 0.6, 1),
                     breaks = c("light", "medium", "full")) +
  scale_fill_manual(values = c("#54F2F2", "#F76F8E", "purple"),
                    breaks = c("red", "blue", "yellow")) +
  ylab("Counterfactual Distance\nto MTX-211 1uM") +
  theme(text = element_text(size = 7, color = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks = element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.1, "cm"),
        axis.title.y = element_blank()) +
  guides(alpha = "none",
         fill = "none") +
  monocle3:::monocle_theme_opts()
ggsave("MrVI_PI3Ki_EGFRi_adaptive_resistance_genes_single_cell_distance_violin_colorchange.png",
       dpi = 600,
       height = 0.75, width = 2.25)
