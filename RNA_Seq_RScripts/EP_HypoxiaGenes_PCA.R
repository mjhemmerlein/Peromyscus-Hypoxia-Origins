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

# Read in genes of interest
EP_Summary = read_xlsx("RNA_Seq_Output/EP_ISO_Ortho_Summary.xlsx")

# Hypoxia
hyp_genes = EP_Summary %>%
  filter(O2_SIG == "SIG" | IXN_SIG == "SIG") %>%
  pull(Pman_GeneID) %>%
  unique()  

hyp_genes = as.vector(hyp_genes)

dPman_EP_hyp = dPman_EP[hyp_genes, ]

plotMDS(dPman_EP_hyp, col = as.numeric(EP_Sample_Info$Group), labels = EP_Sample_Info$Group)

# Get MDS coordinates from your object
mds <- plotMDS(dPman_EP_hyp, plot = FALSE)

# Create a data frame with coordinates and sample info
mds_df <- data.frame(
  Sample = colnames(dPman_EP_hyp),
  Dim1 = mds$x,
  Dim2 = mds$y,
  Group = EP_Sample_Info$Group,
  Treatment = EP_Sample_Info$O2, 
  Pop = EP_Sample_Info$Strain
)

# Plot with ggplot2
ggplot(mds_df, aes(x = Dim1, y = Dim2, color = Pop, shape = Treatment)) +
  geom_point(size = 3.5) +
  labs(
    x = paste("Dimension 1 (", round(mds$var.explained[1] * 100, 1), "%)", sep = ""),
    y = paste("Dimension 2 (", round(mds$var.explained[2] * 100, 1), "%)", sep = ""),
    title = "EP MDS Plot Hypoxia Genes (404)"
  ) +
  theme_classic() +
  scale_shape_manual(values = c(16,5)) +
  scale_color_manual(values = c(
    "BW" = "#F4C552",
    "ME" = "#247BB1")) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5))
