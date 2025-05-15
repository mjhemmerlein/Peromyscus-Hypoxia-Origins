# Packages
library(ggplot2)
library(readxl)
library(dplyr) 
library(edgeR)
library(ggrepel)
library(gprofiler2)
library(readr)
library(tidyverse)

# Load and clean data
EP_ISO_Strain <- read_xlsx("RNA_Seq_Output/EP_ISO_Ortho_Summary.xlsx")

EP_ISO_Strain = EP_ISO_Strain %>%
  select(c(Pman_GeneID, strain_logFC, strain_adj_P.Val))

# Set cutoffs
logFC_cutoff <- 2
pval_cutoff <- 0.05

# Consolidate diffexpressed_EP calculation
EP_ISO_Strain <- EP_ISO_Strain %>%
  mutate(diffexpressed_EP = case_when(
    strain_logFC > logFC_cutoff & strain_adj_P.Val < pval_cutoff ~ "UP",
    strain_logFC < -logFC_cutoff & strain_adj_P.Val < pval_cutoff ~ "DOWN",
    TRUE ~ "NA"
  ))

# Create a new column "delabel" for the top 30 differentially expressed genes
# First, create a temporary data frame with the ordered genes
top_genes <- EP_ISO_Strain %>%
  arrange(strain_adj_P.Val) %>%
  head(30) %>%
  pull(Pman_GeneID)

# Now add the labels
EP_ISO_Strain$delabel <- ifelse(EP_ISO_Strain$Pman_GeneID %in% top_genes, 
                                EP_ISO_Strain$Pman_GeneID, 
                                NA)

# Volcano plot
ggplot(data = EP_ISO_Strain, 
       aes(x = strain_logFC, 
           y = -log10(strain_adj_P.Val), 
           label = delabel)) +
  geom_point() +
  geom_text_repel()

# For EP_ISO_Strain dataset
# Enhanced volcano plot with colored points and cutoff lines
ggplot(data = EP_ISO_Strain, 
       aes(x = strain_logFC, 
           y = -log10(strain_adj_P.Val), 
           color = diffexpressed_EP,
           label = delabel)) +
  geom_point(size = 1.5, alpha = 0.7) +
  geom_text_repel(max.overlaps = 15, 
                  box.padding = 0.5,
                  segment.color = "grey50",
                  show.legend = FALSE) +
  scale_color_manual(values = c("DOWN" = "blue", "UP" = "red", "NA" = "black")) +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "solid", color = "darkgrey") +
  geom_hline(yintercept = -log10(pval_cutoff), linetype = "solid", color = "darkgrey") +
  labs(x = "Log2 Fold Change",
       y = "-Log10 Adjusted p-value",
       color = "Differential\nExpression") +
  theme_bw()


# Load and clean data
LP_ISO_Strain <- read_xlsx("RNA_Seq_Output/LP_ISO_Ortho_Summary.xlsx")

LP_ISO_Strain = LP_ISO_Strain %>%
  select(c(Pman_GeneID, strain_logFC, strain_adj_P.Val))

# Set cutoffs
logFC_cutoff <- 2
pval_cutoff <- 0.05

# Consolidate diffexpressed_LP calculation
LP_ISO_Strain <- LP_ISO_Strain %>%
  mutate(diffexpressed_LP = case_when(
    strain_logFC > logFC_cutoff & strain_adj_P.Val < pval_cutoff ~ "UP",
    strain_logFC < -logFC_cutoff & strain_adj_P.Val < pval_cutoff ~ "DOWN",
    TRUE ~ "NA"
  ))

# Create a new column "delabel" for the top 30 differentially expressed genes
# First, create a temporary data frame with the ordered genes
top_genes <- LP_ISO_Strain %>%
  arrange(strain_adj_P.Val) %>%
  head(30) %>%
  pull(Pman_GeneID)

# Now add the labels
LP_ISO_Strain$delabel <- ifelse(LP_ISO_Strain$Pman_GeneID %in% top_genes, 
                                LP_ISO_Strain$Pman_GeneID, 
                                NA)

# Enhanced volcano plot with colored points and cutoff lines
ggplot(data = LP_ISO_Strain, 
       aes(x = strain_logFC, 
           y = -log10(strain_adj_P.Val), 
           color = diffexpressed_LP,
           label = delabel)) +
  geom_point(size = 1.5, alpha = 0.7) +
  geom_text_repel(max.overlaps = 15, 
                  box.padding = 0.5,
                  segment.color = "grey50",
                  show.legend = FALSE) +
  scale_color_manual(values = c("DOWN" = "blue", "UP" = "red", "NA" = "black")) +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "solid", color = "darkgrey") +
  geom_hline(yintercept = -log10(pval_cutoff), linetype = "solid", color = "darkgrey") +
  labs(x = "Log2 Fold Change",
       y = "-Log10 Adjusted p-value",
       color = "Differential\nExpression") +
  ylim(0,250) +
  theme_bw()
