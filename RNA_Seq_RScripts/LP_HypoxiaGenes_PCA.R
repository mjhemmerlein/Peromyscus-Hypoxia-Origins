
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
dPman_LP <- dPman_0[keep,]
dim(dPman_LP)

# Log-CPM values
logcpm <- cpm(dPman_LP, log = TRUE)

# Read in genes of interest
LP_Summary = read_xlsx("RNA_Seq_Output/LP_ISO_Ortho_Summary.xlsx")

# Hypoxia

hyp_genes = LP_Summary %>%
  filter(O2_SIG == "SIG" | IXN_SIG == "SIG") %>%
  pull(Pman_GeneID) %>%
  unique()  

hyp_genes = as.vector(hyp_genes)

dPman_LP_hyp = dPman_LP[hyp_genes, ]

plotMDS(dPman_LP_hyp, col = as.numeric(LP_Sample_Info$Group), labels = LP_Sample_Info$Group)

# Get MDS coordinates from your object
mds <- plotMDS(dPman_LP_hyp, plot = FALSE)

# Create a data frame with coordinates and sample info
mds_df <- data.frame(
  Sample = colnames(dPman_LP_hyp),
  Dim1 = mds$x,
  Dim2 = mds$y,
  Group = LP_Sample_Info$Group,
  Treatment = LP_Sample_Info$O2, 
  Pop = LP_Sample_Info$Strain
)

# Plot with ggplot2
ggplot(mds_df, aes(x = Dim1, y = Dim2, color = Pop, shape = Treatment)) +
  geom_point(size = 3.5) +
  labs(
    x = paste("Dimension 1 (", round(mds$var.explained[1] * 100, 1), "%)", sep = ""),
    y = paste("Dimension 2 (", round(mds$var.explained[2] * 100, 1), "%)", sep = ""),
    title = "LP MDS Plot Hypoxia Genes (424)"
  ) +
  theme_classic() +
  scale_shape_manual(values = c(16,5)) +
  scale_color_manual(values = c(
    "BW" = "#F4C552",
    "ME" = "#247BB1")) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5))


