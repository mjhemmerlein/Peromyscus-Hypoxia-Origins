# Libraries
library(tidyverse)
library(readr)
library(readxl)
library(pheatmap)
library(BiocParallel)
library(dplyr)
library(ggplot2)
library(readxl)
library(variancePartition)
library(edgeR)
library(Matrix)
library(reshape2)
library(viridis)

# Early pregnancy raw data
EP_Sample_Info <- read_xlsx("RNA_Seq_RawData/MetaData_EP.xlsx")
EP_Sample_Info <- as.data.frame(EP_Sample_Info)
rownames(EP_Sample_Info) <- EP_Sample_Info$Sample_ID
EP_Sample_Info <- EP_Sample_Info %>% filter(Strain == "BW" | Strain == "ME")
EP_Sample_Info$Group <- as.factor(EP_Sample_Info$Group)
EP_Sample_Info$Group_Timepoint = paste0("EP_", EP_Sample_Info$Group)

Pman_rawreads_EP <- read_xlsx("RNA_Seq_RawData/EP_Pman_ExtMMFrac_readcounts_Exon.xlsx")
Pman_rawreads_EP <- as.data.frame(Pman_rawreads_EP)
Pman_rawreads_EP <- Pman_rawreads_EP %>% filter(!is.na(Geneid))
row.names(Pman_rawreads_EP) <- Pman_rawreads_EP$Geneid
Pman_rawreads_EP <- Pman_rawreads_EP[,-c(1:6)]

# Check and match columns
Check_EP <- EP_Sample_Info$Seq_Name
colnames(Pman_rawreads_EP) == Check_EP
colnames(Pman_rawreads_EP) <- rownames(EP_Sample_Info)

# Filter data
Pman_readcounts_EP <- as.matrix(Pman_rawreads_EP)
dPman_0_EP <- DGEList(Pman_readcounts_EP)
dPman_0_EP <- calcNormFactors(dPman_0_EP)
keep_EP <- rowSums(cpm(dPman_0_EP) > 0.5) >= 60
dPman_EP <- dPman_0_EP[keep_EP,]
dim(dPman_EP)

# Log-CPM values
logcpm <- cpm(dPman_EP, log = TRUE)

# Read in genes of interest
EP_Summary = read_xlsx("RNA_Seq_Output/EP_ISO_Ortho_Summary.xlsx")

# All hypoxia sensitive genes (404) --------

gene_class <- EP_Summary %>%
  filter(O2_SIG == "SIG" | IXN_SIG == "SIG") %>%
  mutate(Category = case_when(
    IXN_SIG == "SIG" ~ "IXN_SIG",
    O2_SIG == "SIG" ~ "O2_SIG")) %>%
  distinct(Pman_GeneID, .keep_all = TRUE)

count(gene_class %>%
        filter(Category == "O2_SIG"))

count(gene_class %>%
        filter(Category == "IXN_SIG"))

gene_class$Category <- factor(gene_class$Category, levels = c("O2_SIG", "IXN_SIG"))

logcpm_subset <- logcpm[gene_class$Pman_GeneID, ]
z_scores <- t(scale(t(logcpm_subset)))

annotation_row <- data.frame(Category = gene_class$Category)
row.names(annotation_row) <- gene_class$Pman_GeneID

# Add Group and Sex to column annotations
annotation_col <- data.frame(
  Group = EP_Sample_Info$Group,
  Sex = EP_Sample_Info$Sex
)
rownames(annotation_col) <- colnames(z_scores)

# Order columns by Group first, then Sex within each Group
sample_order <- order(annotation_col$Group, annotation_col$Sex)

# Add visual gaps between Groups
gaps_col <- cumsum(table(annotation_col$Group[sample_order]))

# (Optionally) also show gaps between sexes within each group:
# gaps_col <- cumsum(table(paste(annotation_col$Group, annotation_col$Sex)[sample_order]))

# Cluster genes by category (as before)
clustered_order <- c()
for (cat in levels(gene_class$Category)) {
  subset_genes <- gene_class$Pman_GeneID[gene_class$Category == cat]
  if (length(subset_genes) > 1) {
    hc <- hclust(dist(z_scores[subset_genes, sample_order]))
    clustered_order <- c(clustered_order, subset_genes[hc$order])
  } else {
    clustered_order <- c(clustered_order, subset_genes)
  }
}

gaps_row <- cumsum(table(gene_class$Category))

# Draw heatmap
pheatmap(
  z_scores[clustered_order, sample_order],
  scale = "none",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = annotation_row,
  annotation_col = annotation_col,
  gaps_row = gaps_row,
  gaps_col = gaps_col,
  main = "EP Z-scores of Unique O2/IXN-Significant Genes",
  show_rownames = ifelse(nrow(z_scores) < 50, TRUE, FALSE),
  color = viridis(100, option = "A")
)

save_pheatmap_pdf <- function(pheatmap_result, filename, width = 10, height = 7) {
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(pheatmap_result$gtable)
  dev.off()
}

hm = pheatmap(
  z_scores[clustered_order, sample_order],
  scale = "none",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = annotation_row,
  annotation_col = annotation_col,
  gaps_row = gaps_row,
  gaps_col = gaps_col,
  main = "EP Z-scores of Unique O2/IXN-Significant Genes",
  show_rownames = ifelse(nrow(z_scores) < 50, TRUE, FALSE),
  color = viridis(100, option = "A")
)

save_pheatmap_pdf(hm, "Plots/Heatmaps/EP_Zscore_heatmap.pdf")


