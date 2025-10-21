
# Packages
library(ggplot2)
library(readxl)
library(dplyr) 
library(edgeR)
library(ggrepel)
library(gprofiler2)
library(readr)
library(tidyverse)

# Input and munge 

# Set cutoffs
logFC_cutoff <- 0

# Load and clean data
EP_ISO_Strain <- read_xlsx("RNA_Seq_Output/EP_ISO_Ortho_Summary.xlsx")

EP_ISO_Strain = EP_ISO_Strain %>%
  select(Pman_GeneID, 
         Mus_GeneID, 
         JLamb_GeneID, 
         Jlamb_PmanAttribute, 
         BW_O2_logFC, 
         BW_O2_adj_P.Val, 
         BW_O2_SIG,
         ME_O2_logFC, 
         ME_O2_adj_P.Val, 
         ME_O2_SIG)

# Consolidate diffexpressed_EP calculation
EP_ISO_Strain <- EP_ISO_Strain %>%
  mutate(diffexpressed_BW = case_when(
    BW_O2_logFC > logFC_cutoff ~ "Up",
    BW_O2_logFC < -logFC_cutoff ~ "Down"
  ))

EP_ISO_Strain <- EP_ISO_Strain %>%
  mutate(diffexpressed_ME = case_when(
    ME_O2_logFC > logFC_cutoff ~ "Up",
    ME_O2_logFC < -logFC_cutoff ~ "Down"
  ))

# Add a 'Quadrant' column based on the combination of statuses
EP_ISO_Strain <- EP_ISO_Strain %>%
  mutate(Quadrant =
           paste0(diffexpressed_BW, "-", diffexpressed_ME))


Top_genes = rbind(
  
  EP_ISO_Strain %>%
    filter(Quadrant == "Up-Up") %>%
    arrange(desc(BW_O2_logFC)) %>%
    head(10), 
  
  EP_ISO_Strain %>%
    filter(Quadrant == "Up-Up") %>%
    arrange(desc(ME_O2_logFC)) %>%
    head(10),
  
  EP_ISO_Strain %>%
    filter(Quadrant == "Up-Down") %>%
    arrange(ME_O2_logFC) %>%
    head(4), 
  
  EP_ISO_Strain %>%
    filter(Quadrant == "Up-Down") %>%
    arrange(desc(BW_O2_logFC)) %>%
    head(4),
  
  EP_ISO_Strain %>%
    filter(Quadrant == "Down-Down") %>%
    arrange(ME_O2_logFC) %>%
    head(10), 
  
  EP_ISO_Strain %>%
    filter(Quadrant == "Down-Down") %>%
    arrange(BW_O2_logFC) %>%
    head(10),
  
  EP_ISO_Strain %>%
    filter(Quadrant == "Down-Up") %>%
    arrange(desc(ME_O2_logFC)) %>%
    head(4), 
  
  EP_ISO_Strain %>%
    filter(Quadrant == "Down-Up") %>%
    arrange(BW_O2_logFC) %>%
    head(4))



# Define the logFC cutoffs
lower_cutoff <- -0
upper_cutoff <- 0

# Create a 4-quadrant scatter plot
ggplot(EP_ISO_Strain, aes(x = BW_O2_logFC, y = ME_O2_logFC)) +
  geom_point(size = 3.5) +
  xlim(-3,3) +
  ylim(-3,3) +
  labs(y = "ME\nHypoxia logFC", x = "BW\nHypoxia logFC") +
  geom_hline(yintercept = upper_cutoff, linetype = "solid") +
  geom_hline(yintercept = lower_cutoff, linetype = "solid") +
  geom_vline(xintercept = upper_cutoff, linetype = "solid") +
  geom_vline(xintercept = lower_cutoff, linetype = "solid") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none") +
  theme(
    axis.title.y = element_text(face = "bold", size = 20),
    axis.title.x = element_text(face = "bold", size = 20),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 12))


EP_ISO_Strain <- EP_ISO_Strain %>%
  mutate(Label = ifelse(JLamb_GeneID %in% Top_genes$JLamb_GeneID, JLamb_GeneID, NA))

ggplot(EP_ISO_Strain, aes(x = BW_O2_logFC, y = ME_O2_logFC)) +
  geom_point(alpha = 0.6, size = 3.5) +
  xlim(-3,3) +
  ylim(-3,3) +
  labs(y = NULL, x = NULL) +
  geom_hline(yintercept = upper_cutoff, linetype = "solid") +
  geom_hline(yintercept = lower_cutoff, linetype = "solid") +
  geom_vline(xintercept = upper_cutoff, linetype = "solid") +
  geom_vline(xintercept = lower_cutoff, linetype = "solid") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none") +
  theme(
    axis.text.y = element_text(face = "bold", size = 20),
    axis.title.y = element_text(face = "bold", size = 24),
    axis.title.x = element_text(face = "bold", size = 24),
    axis.text.x = element_text(face = "bold", size = 20)) +
  geom_text_repel(aes(label = Label), 
                  size = 5, 
                  na.rm = TRUE, 
                  color = "black",
                  segment.color = "black", 
                  segment.size = 0.5,     
                  segment.alpha = 0.5,
                  max.overlaps =  Inf) +
  labs(y = "ME\nHypoxia logFC", x = "BW\nHypoxia logFC")


