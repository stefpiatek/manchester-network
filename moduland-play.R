library(tidyverse)
library(here)
library(ggrepel)
library(VarfromPDB)
library(RISmed)
library(stringi)

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

palette_5 <- c( "#16b616", "#A3A3A3", "#B80000", "#FFA347", "#142952")

### Moduland analysis ----

panelapp_green_genes <- read_delim(here("data", "BWS_genes-panelapp.txt"), delim = "\t") %>%
  filter(status == "green")

panelapp_red_genes <- read_delim(here("data", "BWS_genes-panelapp.txt"), delim = "\t") %>%
  filter(status == "red")

data <- read_csv(here("data", "BIOGRID-ORGANISM-Homo_sapiens-3.4.158.tab2.txt-0-bridgeness.csv")) %>%
  inner_join(read_csv(here("data", "BIOGRID-ORGANISM-Homo_sapiens-3.4.158.tab2.txt-0-linkland.csv"))) %>%
  mutate(first_group = if_else(condition = `ModuLand community centrality` > 5000 & `ModuLand bridgeness` > 1,
                               true = nodeID,
                               false = "")) %>%
  mutate(second_group = if_else(condition = `ModuLand community centrality` < 700 & `ModuLand community centrality` > 200 & `ModuLand bridgeness` > 2,
                                true = nodeID,
                                false = "")) %>%
  mutate(third_group = if_else(condition = nodeID == "HRAS",
                                true = nodeID,
                                false = "")) %>%
  mutate(gene_of_interest = paste(first_group, second_group, third_group, sep="")) %>%
  mutate(gene_selected = if_else(condition = gene_of_interest != "",
                                 true = "selected",
                                 false = "")) %>%
  select(-ends_with("_group")) %>%
  mutate(green_genes = if_else(condition = nodeID %in% panelapp_green_genes$`Gene Symbol`,
                               true = "green",
                               false = "")) %>%
  mutate(red_genes = if_else(condition = nodeID %in% panelapp_red_genes$`Gene Symbol`,
                             true = "red",
                             false = "")) %>%
  mutate(panel_app_status = if_else(condition = green_genes == "" & red_genes == "",
                                    true = "not assigned",
                                    false = paste(red_genes, green_genes, sep=""))) %>%
  select(-ends_with("_genes"))

ggplot(data, aes(x = `ModuLand community centrality`, y = `ModuLand bridgeness`)) +
  geom_smooth(method='lm', formula= y ~ x) +
  geom_point(aes(colour=panel_app_status, shape = gene_selected), size = 2.2) +
  scale_x_log10() +
  scale_y_log10() +
  annotation_logticks(base = 10, sides = "bl") +
  geom_label_repel(aes(label = gene_of_interest, fill = panel_app_status),
                   box.padding = 0.35, 
                   point.padding = 0.5,
                   alpha=0.8,  
                   segment.color = 'grey50',
                   nudge_x = 0.2) +
  scale_color_manual(values = palette_5) +
  scale_fill_manual(values = palette_5) +
  guides(shape = FALSE)



#### gene phenotype correlation -----



gene_pheno_matrix <- read_csv(here("data", "gene_pheno_matrix.csv")) %>%
  as.data.frame()
row.names(gene_pheno_matrix) <- gene_pheno_matrix$X1

gene_pheno_matrix <- gene_pheno_matrix[, 2:ncol(gene_pheno_matrix)]
gene_pheno_matrix <- as.matrix(gene_pheno_matrix)

matrix <- (gene_pheno_matrix %*% t(gene_pheno_matrix)) -1

library(gplots)
heatmap.2(matrix,trace="none")

hc.rows<-hclust(dist(t(matrix)))
h<-4
plot(hc.rows,cex=0.5)
ct<-cutree(hc.rows,h)
rect.hclust(hc.rows,h)

genes<-names(ct[which(ct==4|ct==3)])
important_gene_pheno <- gene_pheno_matrix[match(genes,rownames(gene_pheno_matrix)),]
important_gene_pheno <- important_gene_pheno[, colSums(important_gene_pheno) != 0]

rowSums(important_gene_pheno)

phenotypes <- read_delim(here("data", "ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt"), delim= "\t", comment = "#",
                         col_names = c("HPO", "term_name", "gene_id", "gene_name"))

phenotypes %>%
  select(HPO, term_name) %>%
  distinct() %>%
  filter(HPO %in% colnames(important_gene_pheno)) %>%
  as.data.frame()
