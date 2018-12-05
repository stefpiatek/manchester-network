library(tidyverse)
library(here)
library(ggrepel)

### Set custom theme ----

theme_custom <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      panel.grid.minor.y = element_line(colour="white"),
      panel.border = element_blank(),
      legend.key = element_rect(fill = NA, colour = NA),
      axis.title.x=element_text(vjust=-1.5),
      plot.margin= unit(c(0.5, 0, 0.5, 0), "cm")# top right bottom left
    )
}

theme_set(theme_custom())

palette_18 <- c("#000000","#4286ff","#b3a80a","#D3D3D3","#6e0068","#01dca1","#f44575","#005aac","#f58037","#003888","#ff8868","#c492ff","#884600","#60004c","#d48267","#ff7ea9","#808080", "#5fe17a")

### Moduland analysis ----

data <- read_csv(here("data", "BIOGRID-ORGANISM-Homo_sapiens-3.4.158.tab2.txt-0-bridgeness.csv")) %>%
  inner_join(read_csv(here("data", "BIOGRID-ORGANISM-Homo_sapiens-3.4.158.tab2.txt-0-linkland.csv"))) %>%
  mutate(first_group = if_else(condition = `ModuLand community centrality` > 5000 & `ModuLand bridgeness` > 1,
                               true = nodeID,
                               false = "")) %>%
  mutate(second_group = if_else(condition = `ModuLand community centrality` < 700 & `ModuLand community centrality` > 200 & `ModuLand bridgeness` > 2,
                                true = nodeID,
                                false = "")) %>%
  mutate(gene_of_interest = paste(first_group, second_group, sep="")) %>%
  select(-ends_with("_group"))

ggplot(data, aes(x = `ModuLand community centrality`, y = `ModuLand bridgeness`)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  annotation_logticks(base = 10, sides = "bl") +
  geom_smooth(method='lm', formula= y ~ x) +
  geom_label_repel(aes(label = paste(gene_of_interest)),
                   box.padding = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') 


data %>%
  filter(`ModuLand community centrality` > 5000) %>%
  filter(`ModuLand bridgeness` > 1)

data %>%
  filter(`ModuLand community centrality` < 700) %>%
  filter(`ModuLand community centrality` > 200) %>%
  filter(`ModuLand bridgeness` > 2)

data %>% 
  filter(nodeID == "IGF2")
