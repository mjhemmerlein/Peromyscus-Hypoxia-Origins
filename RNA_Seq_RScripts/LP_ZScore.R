
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

# Late Pregnancy ------------
LP_Sample_Info <- read_xlsx("RNA_Seq_RawData/MetaData_LP.xlsx")
LP_Sample_Info = as.data.frame(LP_Sample_Info)
LP_Sample_Info = LP_Sample_Info[-c(80,81),]
rownames(LP_Sample_Info) = LP_Sample_Info$Sample_ID_LZ
LP_Sample_Info = LP_Sample_Info %>% filter(Strain == "BW" | Strain == "ME")
LP_Sample_Info = LP_Sample_Info[-which(rownames(LP_Sample_Info) == "LZ089"),]

# Read in Files + QC
Pman_rawreads = read_xlsx("RNA_Seq_RawData/LP_Pman_ExtMMFrac_readcounts_Exon.xlsx")
Pman_rawreads = as.data.frame(Pman_rawreads)
Pman_rawreads = `row.names<-`(Pman_rawreads, Pman_rawreads$Geneid)
Pman_rawreads <- Pman_rawreads[,-c(1:6)]
Pman_rawreads <- subset(Pman_rawreads, select = -c(RNA201216ZC_LZ089_S16_L001_fastp_pman_Halign_liberal.bam))

# Check read table vs sample info
Check = LP_Sample_Info$Seq_Name
colnames(Pman_rawreads) == Check

colnames(Pman_rawreads) = rownames(LP_Sample_Info)

Pman_readcounts <- as.matrix(Pman_rawreads)
dPman_0 <- DGEList(Pman_readcounts)

dPman_0 <- calcNormFactors(dPman_0)
dim(dPman_0)
keep <- rowSums(cpm(dPman_0) > 0.5 ) >= 60
dPman <- dPman_0[keep,]
dim(dPman)


# Log-CPM values
logcpm <- cpm(dPman, log = TRUE)

# Read in genes of interest
LP_Summary = read_xlsx("RNA_Seq_Output/LP_ISO_Ortho_Summary.xlsx")


# All hypoxia sensitive genes (404) -------

gene_class <- LP_Summary %>%
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

annotation_col <- data.frame(Group = LP_Sample_Info$Group)
row.names(annotation_col) <- colnames(z_scores)

sample_order <- order(annotation_col$Group)
gaps_col <- cumsum(table(annotation_col$Group[sample_order]))

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

pheatmap(z_scores[clustered_order, sample_order],
         scale = "none",
         cluster_rows = FALSE,   # we already clustered manually
         cluster_cols = FALSE,   # we’re ordering by group
         annotation_row = annotation_row,
         annotation_col = annotation_col,
         gaps_row = gaps_row,
         gaps_col = gaps_col,
         main = "LP Z-scores of Unique O2/IXN-Significant Genes",
         show_rownames = ifelse(nrow(z_scores) < 50, TRUE, FALSE),
         color = viridis(100, option = "A"))
