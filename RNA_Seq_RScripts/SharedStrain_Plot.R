
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
pval_cutoff <- 0.05


# Load and clean data
EP_ISO_Strain <- read_xlsx("RNA_Seq_Output/EP_ISO_Ortho_Summary.xlsx")

EP_ISO_Strain = EP_ISO_Strain %>%
  filter(strain_SIG == "SIG") %>%
  rename(EP_strain_logFC = strain_logFC) %>%
  rename(EP_strain_adj_P.Val = strain_adj_P.Val) %>%
  rename(EP_strain_SIG = strain_SIG) %>%
  select(Pman_GeneID, Mus_GeneID, JLamb_GeneID, Jlamb_PmanAttribute, EP_strain_logFC, EP_strain_adj_P.Val, EP_strain_SIG)

# Consolidate diffexpressed_EP calculation
EP_ISO_Strain <- EP_ISO_Strain %>%
  mutate(diffexpressed_EP = case_when(
    EP_strain_logFC > logFC_cutoff & EP_strain_adj_P.Val < pval_cutoff ~ "UP",
    EP_strain_logFC < -logFC_cutoff & EP_strain_adj_P.Val < pval_cutoff ~ "DOWN",
    TRUE ~ "NA"
  ))

# Load and clean data
LP_ISO_Strain <- read_xlsx("RNA_Seq_Output/LP_ISO_Ortho_Summary.xlsx")

LP_ISO_Strain = LP_ISO_Strain %>%
  filter(strain_SIG == "SIG") %>%
  rename(LP_strain_logFC = strain_logFC) %>%
  rename(LP_strain_adj_P.Val = strain_adj_P.Val) %>%
  rename(LP_strain_SIG = strain_SIG) %>%
  select(Pman_GeneID, LP_strain_logFC, LP_strain_adj_P.Val, LP_strain_SIG)

# Consolidate diffexpressed_LP calculation
LP_ISO_Strain <- LP_ISO_Strain %>%
  mutate(diffexpressed_LP = case_when(
    LP_strain_logFC > logFC_cutoff & LP_strain_adj_P.Val < pval_cutoff ~ "UP",
    LP_strain_logFC < -logFC_cutoff & LP_strain_adj_P.Val < pval_cutoff ~ "DOWN",
    TRUE ~ "NA"
  ))

count(intersect(
  EP_ISO_Strain %>% filter(EP_strain_SIG == "SIG") %>% select(Pman_GeneID),
  LP_ISO_Strain %>% filter(LP_strain_SIG == "SIG") %>% select(Pman_GeneID)
))


EP_LP_strain = full_join(EP_ISO_Strain, LP_ISO_Strain, by = c("Pman_GeneID"))
EP_LP_strain = na.omit(EP_LP_strain)

# Add a 'Quadrant' column based on the combination of statuses
EP_LP_strain <- EP_LP_strain %>%
  mutate(Quadrant = paste(diffexpressed_EP, diffexpressed_LP, sep = "-"))

Strain_summary <- EP_LP_strain %>%
  group_by(Quadrant) %>%
  summarise(Freq = n()) %>%
  mutate(Percent = Freq / sum(Freq) * 100)


# write.csv(EP_LP_strain, file = "RNA_Seq_Output/EP_LP_SharedStrain.csv")


quadrant_colors <- c("UP-UP" = "#3b8132", 
                     "UP-DOWN" = "darkgray", 
                     "DOWN-DOWN" = "#8bca84", 
                     "DOWN-UP" = "#262626")

# Create a 4-quadrant scatter plot
ggplot(EP_LP_strain, aes(x = EP_strain_logFC, y = LP_strain_logFC, color = Quadrant)) +
  geom_point(size = 5, alpha = 0.4, stroke = 0) +
  geom_hline(yintercept = 0, linetype = "solid") +
  geom_vline(xintercept = 0, linetype = "solid") +
  xlim(-6.5,6.5) +
  ylim(-6.5,6.5) +
  scale_color_manual(values = quadrant_colors) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none") +
  theme(
    axis.text.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank())

ggsave("Plots/Quadrant_Plots/Unlabeled_Quad.pdf", width = 6, height = 6, units = "in", dpi = 300)




# Top 3 genes per quadrant for labeling

Top_genes = rbind(
  
  EP_LP_strain %>%
    filter(Quadrant == "UP-UP") %>%
    arrange(desc(LP_strain_logFC)) %>%
    head(10), 
  
  EP_LP_strain %>%
    filter(Quadrant == "UP-UP") %>%
    arrange(desc(EP_strain_logFC)) %>%
    head(10),
  
  EP_LP_strain %>%
    filter(Quadrant == "UP-DOWN") %>%
    arrange(LP_strain_logFC) %>%
    head(4), 
  
  EP_LP_strain %>%
    filter(Quadrant == "UP-DOWN") %>%
    arrange(desc(EP_strain_logFC)) %>%
    head(4),
  
  EP_LP_strain %>%
    filter(Quadrant == "DOWN-DOWN") %>%
    arrange(LP_strain_logFC) %>%
    head(10), 
  
  EP_LP_strain %>%
    filter(Quadrant == "DOWN-DOWN") %>%
    arrange(EP_strain_logFC) %>%
    head(10),
  
  EP_LP_strain %>%
    filter(Quadrant == "DOWN-UP") %>%
    arrange(desc(LP_strain_logFC)) %>%
    head(4), 
  
  EP_LP_strain %>%
    filter(Quadrant == "DOWN-UP") %>%
    arrange(EP_strain_logFC) %>%
    head(4))


# Quadrant w/ labels 

EP_LP_strain <- EP_LP_strain %>%
  mutate(Label = ifelse(JLamb_GeneID %in% Top_genes$JLamb_GeneID, JLamb_GeneID, NA))

ggplot(EP_LP_strain, aes(x = EP_strain_logFC, y = LP_strain_logFC, color = Quadrant)) +
  geom_point(alpha = 0.6, size = 3.5) +
  geom_hline(yintercept = 0, linetype = "solid") +
  geom_vline(xintercept = 0, linetype = "solid") +
  xlim(-6.5,6.5) +
  ylim(-6.5,6.5) +
  labs(y = NULL, x = NULL) +
  scale_color_manual(values = quadrant_colors) +
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
  labs(y = "Late Pregnancy\npopulation logFC", x = "Early Pregnancy\npopulation logFC")

# ggsave("Dream_Output/Plots/Quadrant_Plots/AllgenesQuadrant.png", width = 10, height = 10, units = "in", dpi = 300)