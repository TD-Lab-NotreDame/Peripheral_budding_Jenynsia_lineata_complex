library(tidyverse)
library(readr)
library(patchwork)

#Import and clean raw data files
# need to manually remove the "#" that is in the first line of each .pestPG file
file_paths <- list.files(path = "data", pattern = "*.thetas.idx.pestPG", full.names = TRUE)
# Use map to read in the data and add population identifier from file name
combined_data <- file_paths %>%
  # Read each file
  map_df(~ read.table(.x, header = TRUE) %>% 
           # Extract the population identifier from the file name (remove file extension)
           mutate(Population = gsub(".*/|\\..*", "", .x)))

combined_data <- combined_data %>%
  filter(nSites > 0)
combined_data <- combined_data %>%
  mutate(Contig = as.numeric(gsub("^dDocent_Contig_", "", Chr))) 


df <- combined_data %>%
  extract(
    X.indexStart.indexStop..firstPos_withData.lastPos_withData..WinStart.WinStop.,
    into = c("indexStartStop", "firstLast", "winStartStop"),
    regex = "\\(([^)]+)\\)\\(([^)]+)\\)\\(([^)]+)\\)"
  )

df <- df %>%
  separate(indexStartStop, into = c("indexStart", "indexStop"), sep = ",") %>%
  separate(firstLast, into = c("firstPos", "lastPos"), sep = ",") %>%
  separate(winStartStop, into = c("winStart", "winStop"), sep = ",")


chr_lengths <- df %>%
  group_by(Contig) %>%
  summarise(chr_len = max(winStop)) %>%
  mutate(chr_start = lag(cumsum(chr_len), default = 0))

df <- df %>%
  left_join(chr_lengths, by = "Contig") %>%
  mutate(pos = WinCenter + chr_start)
df <- df %>% filter(nSites > 0)  # Only keep windows with data
df$pi <- df$tP / df$nSites  # Convert to per-site π


df <- df %>%
  mutate(Population = case_when(
    Population=="luxata_thetas_per_site"~"luxata",
    Population=="lineataEW_thetas_per_site"~"lineataEW",
    Population=="lineataN_thetas_per_site"~"lineataN",
    Population=="lineataW_thetas_per_site"~"lineataW",
    Population=="lineataE_thetas_per_site"~"lineataE",
    Population=="darwini_thetas_per_site"~"darwini",
    Population=="onca_thetas_per_site"~"onca",
    TRUE ~ Population))

summary_stats <- df %>%
  filter(!is.na(tP), nSites > 0) %>%
  mutate(pi_per_site = tP / nSites) %>%
  group_by(Population) %>%
  summarise(
    mean_pi = weighted.mean(pi_per_site, w = nSites),
    se_pi = sqrt(Hmisc::wtd.var(pi_per_site, weights = nSites)) / sqrt(n()),
    mean_tajima = weighted.mean(Tajima, w = nSites),
    se_tajima = sqrt(Hmisc::wtd.var(Tajima, weights = nSites)) / sqrt(n())
  )

summary_stats <- df %>%
  group_by(Population) %>%
  summarise(
    mean_pi = mean(pi, na.rm = TRUE),
    se_pi = sd(pi, na.rm = TRUE) / sqrt(n()))
summary_stats <- summary_stats %>%
  mutate(group = case_when(
    Population == "darwini" ~ "bud",
    Population == "lineataE" ~ "ancestral",
    Population == "lineataW" ~ "ancestral",
    Population == "onca" ~ "bud",
    Population == "luxata" ~ "bud",
    Population == "lineataN" ~ "bud",
    Population == "lineataEW" ~ "ancestral",
    TRUE ~ Population
  )) %>%
  mutate(Population = factor(Population,
                             levels = c(
                               "lineataEW", "lineataE", "lineataW",  "darwini", "onca", "luxata", "lineataN")))

#graphing----
##pi----
jenynsia_pi<-ggplot(summary_stats, aes(x = Population, y = mean_pi, color = Population)) +
  geom_point(size=1.5) +
  geom_errorbar(aes(ymin = mean_pi - se_pi, ymax = mean_pi + se_pi), width = 0.4) +
  scale_color_manual(
    values = c("green3","#6CC7AB", "#6a59cdff","#EE9A00","#4D4D4D","#B95251","#1c90b9ff"))+
  labs(
    x = "Population",
    y = expression("Mean " * pi)
  ) +ylim(0, 0.005) +   theme_minimal() + theme(axis.title = element_text(size = 12),
                                                plot.title = element_text(size = 12),
                                                legend.text = element_text(size = 12),
                                                legend.title = element_text(size = 12),
                                                axis.text.x = element_text(size = 12),
                                                axis.text.y = element_text(size = 12),
                                                legend.key = element_rect(fill = NA)) +
  theme(legend.position =  "bottom",  legend.direction = "horizontal") + theme(axis.text = element_text(size = 12)) 




##tajimasD----
jenynsia_tajimasD<-ggplot(summary_stats, aes(x = Population, y = mean_tajima, color = Population)) +
  geom_point(size=1.5) +
  geom_errorbar(aes(ymin = mean_tajima - se_tajima, ymax = mean_tajima + se_tajima), width = 0.4) +
  scale_color_manual(
    values = c("green3","#6CC7AB", "#6a59cdff","#EE9A00","#4D4D4D","#B95251","#1c90b9ff"))+
  labs(
    x = "Population",
    y = expression("Mean Tajima's D")) +
    ylim(0, 0.25) + theme_minimal() + theme(axis.title = element_text(size = 12),
                                         plot.title = element_text(size = 12),
                                         legend.text = element_text(size = 12),
                                         legend.title = element_text(size = 12),
                                         axis.text.x = element_text(size = 12),
                                         axis.text.y = element_text(size = 12),
                                         legend.key = element_rect(fill = NA)) +
    theme(legend.position =  "bottom",  legend.direction = "horizontal") + theme(axis.text = element_text(size = 12))


ggsave("Jenynsia_tajimasD.pdf", plot = jenynsia_tajimasD, width = 3.5, height = 4, dpi = 300, units = "in")
ggsave("Jenynsia_pi.pdf", plot = jenynsia_pi, width = 3.5, height = 4, dpi = 300, units = "in")

Jenynsia_combined_T_pi<-(jenynsia_pi/jenynsia_tajimasD) + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave("Jenynsia_combined_T_pi.pdf", plot = Jenynsia_combined_T_pi, width = 7, height = 8, dpi = 300, units = "in")
