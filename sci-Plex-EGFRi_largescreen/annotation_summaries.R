library(tidyverse)
this.dir <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(this.dir))
rm(this.dir)

annotations <- read_csv("Annotations_final_RG.csv") %>%
  filter(!drug %in% c("Puromycin", "NSC228155")) # filter controls for EGFRi visualization

reversible <- annotations %>%
  dplyr::rename(reverse = reversible) %>%
  group_by(reverse) %>%
  summarise(count = n())
reversible <- reversible %>%
  arrange(desc(reverse)) %>%
  mutate(prop = count / sum(reversible$count) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

ggplot(reversible, aes(x = "", y = prop, fill = reverse)) +
  geom_bar(stat = "identity", alpha = 0.8, color = "white", width = 1, linewidth = 0.1) +
  coord_polar("y", start = 0) +
  geom_text(aes(label = count), position = position_stack(vjust = 0.5), size = 1.5,
            color = case_when(reversible$reverse == "Yes" ~ "white", 
                              TRUE ~ "black")) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_line(linewidth = 0),
        text = element_text(size = 6),
        legend.key.size = unit(0.2, "cm"),
        legend.box.margin = unit(c(-3,0,0,-3), "mm"),
        axis.title = element_blank(),
        axis.line = element_line(linetype=0),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_fill_manual(name = "",
                    breaks = c("No", "Yes", NA),
                    labels = c("Covalent", "Reversible", "Unknown"),
                    values = c("red", "black")) +
  monocle3:::monocle_theme_opts() 
ggsave("drug_annotations/reversible_pie_chart.png", dpi = 600, height = 0.8, width = 1.3)


generation <- annotations %>%
  dplyr::rename(reverse = reversible) %>%
  group_by(generation) %>%
  summarise(count = n())
generation <- generation %>%
  mutate(generation = as.character(generation)) %>%
  arrange(desc(generation)) %>%
  mutate(prop = count / sum(reversible$count) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

ggplot(generation, aes(x = "", y = prop, fill = generation)) +
  geom_bar(stat = "identity", alpha = 0.5, color = "white", width = 1, linewidth = 0.1) +
  coord_polar("y", start = 0) +
  geom_text(aes(label = count), position = ggpp::position_stacknudge(vjust = 0.5, x = 0.2), 
            size = 1.5) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_line(linewidth = 0),
        text = element_text(size = 6),
        legend.key.size = unit(0.2, "cm"),
        legend.box.margin = unit(c(-3,0,0,-3), "mm"),
        axis.title = element_blank(),
        axis.line = element_line(linetype=0),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_fill_manual(name = "",
                    breaks = c(1, 2, 3, 4, NA),
                    labels = c("1st", "2nd", "3rd", "4th", NA),
                    values = c(RColorBrewer::brewer.pal(n = 4, "Set2"))) +
  monocle3:::monocle_theme_opts() 
ggsave("drug_annotations/generation_pie_chart.png", dpi = 600, height = 0.8, width = 1.3)

class <- annotations %>%
  group_by(class_broad) %>%
  summarise(count = n()) %>%
  # mutate(class = case_when(
  #   count > 1 ~ class_broad,
  #   is.na(class_broad) == NA ~ class_broad,
  #   TRUE ~ "other")) %>%
  mutate(class = case_when(
    count > 1 ~ class_broad,
    class_broad %in% c("PROTAC", 
                       "monoclonal antibody") ~ class_broad,
    TRUE ~ "other")) %>%
  # dplyr::rename(class = class_broad) %>%
  group_by(class) %>%
  summarise(count = sum(count))
classes <- class %>% 
  filter(class != "other") %>% 
  select(class) %>%
  unique() %>% pull()
class <- class %>%
  mutate(class = factor(class, levels = c(classes, "other"))) %>%
  arrange(desc(class)) %>%
  mutate(prop = count / sum(reversible$count) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop ) %>%
  mutate(count_label = ifelse(count >= 2, as.character(count), NA)) %>%
  mutate(count_label_special = ifelse(class %in% c("PROTAC", "monoclonal antibody"), 1, NA))

ggplot(class, aes(x = "", y = prop, fill = class)) +
  geom_bar(stat = "identity", alpha = 0.5, color = "white", width = 1, linewidth = 0.1,
           show.legend = T) +
  coord_polar("y", start = 0) +
  geom_text(aes(label = count_label), position = ggpp::position_stacknudge(vjust = 0.5, x = 0.25), 
            size = 1.5) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_line(linewidth = 0),
        text = element_text(size = 6),
        legend.key.size = unit(0.2, "cm"),
        legend.box.margin = unit(c(-3,0,0,-3), "mm"),
        axis.title = element_blank(),
        axis.line = element_line(linetype=0),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_fill_manual(name = "",values = c(RColorBrewer::brewer.pal(n = 6, "Set1"), RColorBrewer::brewer.pal(n = 4, "Set3"))) +
  monocle3:::monocle_theme_opts() 
ggsave("drug_annotations/class_pie_chart_legend.png", dpi = 600, height = 2, width = 1.5)


specificity <- annotations %>%
  select(specificity_RG) %>%
  separate(col = specificity_RG, sep = "_", into = c("X1", "X2", "X3")) %>%
  pivot_longer(cols = X1:X3, names_to = "X1", values_to = "specificity") %>%
  group_by(specificity) %>%
  summarise(count = n()) %>% 
  drop_na() %>%
  filter(specificity %in% c("HER2", "HER3", "HER4", "ABL", "ALK", "BTK", "HDAC", "PI3K", "RAF", "VEGFR")) %>%
  mutate(specificity = fct_reorder(specificity, desc(count))) %>%
  mutate(count = ifelse(count >= 10, 10, count))

ggplot(specificity, aes(x = specificity, y = as.double(count))) +
  geom_bar(stat = "identity", aes(fill = specificity), show.legend = F) +
  scale_y_continuous(breaks = c(0,5,10),
                     labels = c("0", "5", ">10"),
                     expand = c(0,0)) +
  # geom_text(aes(label = count), hjust = 0) +
  theme(
        # axis.text.x = element_blank(),
        # axis.text.y = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        text = element_text(size = 6),
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_blank(),
        axis.title = element_blank()) +
  scale_fill_manual(values = c(RColorBrewer::brewer.pal(n = 6, "Accent"), 
                                RColorBrewer::brewer.pal(n = 4, "Pastel1"))) +
  coord_flip() +
  monocle3:::monocle_theme_opts()
ggsave("drug_annotations/specificity_bar_graph.png", dpi = 600, height = 0.9, width = 1.1)

FDA <- annotations %>%
  group_by(FDA_approved) %>%
  mutate(FDA_approved = case_when(
    is.na(FDA_approved) == TRUE ~ "No",
    TRUE ~ "Yes"
  )) %>%
  summarise(count = n())
FDA <- FDA %>%
  mutate(FDA_approved = as.character(FDA_approved)) %>%
  arrange(desc(FDA_approved)) %>%
  mutate(prop = count / sum(FDA$count) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

ggplot(FDA, aes(x = "", y = prop, fill = FDA_approved)) +
  geom_bar(stat = "identity", alpha = 0.5, color = "white", width = 1, linewidth = 0.1) +
  coord_polar("y", start = 0) +
  geom_text(aes(label = count), position = ggpp::position_stacknudge(vjust = 0.5, x = 0.2), 
            size = 1.5) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_line(linewidth = 0),
        text = element_text(size = 6),
        legend.key.size = unit(0.2, "cm"),
        legend.box.margin = unit(c(-3,0,0,-3), "mm"),
        axis.title = element_blank(),
        axis.line = element_line(linetype=0),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_fill_manual(name = "",
                    breaks = c("No", "Yes"),
                    labels = c("Not Approved", "Approved"),
                    values = c(RColorBrewer::brewer.pal(n = 6, "Dark2")[c(2,5)])) +
  monocle3:::monocle_theme_opts() 
ggsave("drug_annotations/FDA_approved_pie_chart.png", dpi = 600, height = 0.8, width = 1.3)

