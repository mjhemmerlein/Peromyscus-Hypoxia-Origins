
# Libraries
library(tidyverse)
library(readr)
library(readxl)
library(pheatmap)
library('BiocParallel')
library('dplyr')
library('ggplot2')
library('lme4')
library('readxl')
library('variancePartition')
library('edgeR')
library('Matrix')

# Combined Analysis  ------------------
Sample_Info <- read_xlsx("RNA_Seq_RawData/MetaData_EP.xlsx")
Sample_Info = as.data.frame(Sample_Info)
rownames(Sample_Info) = Sample_Info$Sample_ID
Sample_Info = Sample_Info %>% filter(Strain == "BW" | Strain == "ME")

# Read in Files + QC
Pman_rawreads = read_xlsx("RNA_Seq_RawData/EP_Pman_ExtMMFrac_readcounts.xlsx")
Pman_rawreads = as.data.frame(Pman_rawreads)
Pman_rawreads = Pman_rawreads %>%
  filter(!is.na(Geneid))
Pman_rawreads = `row.names<-`(Pman_rawreads, Pman_rawreads$Geneid)
Pman_rawreads <- Pman_rawreads[,-c(1:6)]

# Check read table vs sample info
Check = Sample_Info$Seq_Name
colnames(Pman_rawreads) == Check

colnames(Pman_rawreads) = rownames(Sample_Info)

Pman_readcounts <- as.matrix(Pman_rawreads)
dPman_0 <- DGEList(Pman_readcounts)

dPman_0 <- calcNormFactors(dPman_0)
dim(dPman_0)
keep <- rowSums(cpm(dPman_0) > 0.5 ) >= 60
dPman <- dPman_0[keep,]
dim(dPman)

# List of glycolysis-related genes
Gly_genes <- c("Gpi", "Pfkm", "Aldoa", "LOC102904208", "LOC102916082", 
               "Pgam1", "Eno1", "LOC102923285", "LOC102928417", "Ldha", "LOC102927525")

Gly_genes <- c("Etfdh", "Etfa", "Etfb", "Abcd4","LOC102923527", "LOC102923832", "Tgfbi")

Gly_genes <- c("LOC102927525", "Etfdh", "Etfa", "Etfb")

Gly_genes <- c("Gpi", "Pfkm", "Aldoa", "LOC102904208", "LOC102916082", 
               "Pgam1", "Eno1", "LOC102923285", "LOC102928417", "Ldha", "LOC102927525")

# Compute logCPM matrix from dPman
logCPM_matrix <- cpm(dPman, log = TRUE)

# Subset logCPM_matrix for Gly_genes
Gly_matrix <- logCPM_matrix[rownames(logCPM_matrix) %in% Gly_genes, ]

# Reorder the rows in Gly_matrix based on Gly_genes order
Gly_matrix <- Gly_matrix[Gly_genes, ]

# Initialize empty matrix to hold group averages
group_avg <- sapply(unique(Sample_Info$Group), function(g) {
  # Identify the samples that belong to this group
  samples <- rownames(Sample_Info)[Sample_Info$Group == g]
  # Compute row means for these samples
  rowMeans(Gly_matrix[, samples, drop = FALSE])
})

# Convert to matrix
group_avg <- as.matrix(group_avg)

# Add row names (gene names) and column names (group names)
rownames(group_avg) <- rownames(Gly_matrix);
rownames(group_avg)[rownames(group_avg) == "LOC102904208"] <- "GAPDH-Like_1"
rownames(group_avg)[rownames(group_avg) == "LOC102916082"] <- "GAPDH-Like_2"
rownames(group_avg)[rownames(group_avg) == "LOC102923285"] <- "PKM-Like"
rownames(group_avg)[rownames(group_avg) == "LOC102928417"] <- "LDHA-Like"

colnames(group_avg) <- unique(Sample_Info$Group);

# Reorder the columns in your group_avg matrix
desired_order <- c("BW1N", "BW2H", "ME1N", "ME2H")
group_avg <- group_avg[, desired_order]

colnames(group_avg) <- c("Lowland Normoxia", 
                         "Lowland Hypoxia", 
                         "Highland Normoxia",
                         "Highland Hypoxia")

# Plot the heatmap with the reordered matrix
pheatmap(group_avg, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE, 
         scale = "none", 
         show_rownames = TRUE)


