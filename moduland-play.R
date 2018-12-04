library(tidyverse)
library(here)

data <- read_csv(here("data", "string_interactions.tsv-0-bridgeness.csv")) %>%
  inner_join(read_csv(here("data", "string_interactions.tsv-0-linkland.csv")))

ggplot(data, aes(x = `ModuLand community centrality`, y = `ModuLand bridgeness`)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method='lm', formula= y ~ x)

data %>%
  filter(`ModuLand bridgeness` > 1) %>%
  filter(`ModuLand community centrality` > 5000)


data %>%
  filter(`ModuLand bridgeness` < 0.3) %>%
  filter(`ModuLand community centrality` > 10000)

data %>%
  filter(`ModuLand bridgeness` < 1e-05)
