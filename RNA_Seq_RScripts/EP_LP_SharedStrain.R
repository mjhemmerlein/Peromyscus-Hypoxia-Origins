
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
  select(Gene_ID, Pman_GeneID, P_maniculatus_attr, M_musculus_attr, EP_strain_logFC, EP_strain_adj_P.Val, EP_strain_SIG)

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
  select(Gene_ID, LP_strain_logFC, LP_strain_adj_P.Val, LP_strain_SIG)

# Consolidate diffexpressed_LP calculation
LP_ISO_Strain <- LP_ISO_Strain %>%
  mutate(diffexpressed_LP = case_when(
    LP_strain_logFC > logFC_cutoff & LP_strain_adj_P.Val < pval_cutoff ~ "UP",
    LP_strain_logFC < -logFC_cutoff & LP_strain_adj_P.Val < pval_cutoff ~ "DOWN",
    TRUE ~ "NA"
  ))

count(intersect(
  EP_ISO_Strain %>% filter(EP_strain_SIG == "SIG") %>% select(Gene_ID),
  LP_ISO_Strain %>% filter(LP_strain_SIG == "SIG") %>% select(Gene_ID)
))


EP_LP_strain = full_join(EP_ISO_Strain, LP_ISO_Strain, by = c("Gene_ID"))
EP_LP_strain = na.omit(EP_LP_strain)

Strain_summary = EP_LP_strain %>%
  group_by(diffexpressed_EP, diffexpressed_LP) %>%
  summarise(Freq = n())


# Directions

# Up_Up
Up_Up = EP_LP_strain %>% filter(diffexpressed_EP == "UP", diffexpressed_LP == "UP") %>%
  arrange(desc(EP_strain_logFC)) 

# Down_Down
Down_Down = EP_LP_strain %>% filter(diffexpressed_EP == "DOWN", diffexpressed_LP == "DOWN") %>%
  arrange(EP_strain_logFC) 


# Quadrant

# Define the logFC cutoffs
lower_cutoff <- -0.5
upper_cutoff <- 0.5

# Define color mapping for the quadrants, including gray for the cutoff range
quadrant_colors <- c("Up-Up" = "tomato2", 
                     "Up-Down" = "darkolivegreen", 
                     "Down-Down" = "deepskyblue3", 
                     "Down-Up" = "blueviolet")

# Add a 'Quadrant' column based on the combination of statuses
EP_LP_strain <- EP_LP_strain %>%
  mutate(Quadrant = case_when(
    EP_strain_logFC < lower_cutoff & LP_strain_logFC < lower_cutoff ~ "Down-Down",
    EP_strain_logFC < lower_cutoff & LP_strain_logFC > upper_cutoff ~ "Down-Up",
    EP_strain_logFC > upper_cutoff & LP_strain_logFC < lower_cutoff ~ "Up-Down",
    EP_strain_logFC > upper_cutoff & LP_strain_logFC > upper_cutoff ~ "Up-Up",
    EP_strain_logFC >= lower_cutoff & EP_strain_logFC <= upper_cutoff & 
    LP_strain_logFC >= lower_cutoff & LP_strain_logFC <= upper_cutoff ~ "InCutoff",  # Gray for points in cutoff range
    TRUE ~ "NA"
  ))

# Output excel

Strain_summary = EP_LP_strain %>%
  group_by(Quadrant) %>%
  summarise(Freq = n())


# write.csv(EP_LP_strain, file = "RNA_Seq_Output/EP_LP_SharedStrain.csv")


# Create a 4-quadrant scatter plot
ggplot(EP_LP_strain, aes(x = EP_strain_logFC, y = LP_strain_logFC, color = Quadrant)) +
  geom_point(size = 3.5) +
  xlim(-6.5,6.5) +
  ylim(-6.5,6.5) +
  labs(y = "Late Pregnancy\npopulation logFC", x = "Early Pregnancy\npopulation logFC") +
  scale_color_manual(values = quadrant_colors) +
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

ggsave("Plots/Quadrant_Plots/Unlabeled_Quad.pdf", width = 7, height = 7, units = "in", dpi = 300)

# Top 3 genes per quadrant for labeling

Top_genes = rbind(EP_LP_strain %>%
  filter(Quadrant == "Up-Up") %>%
  arrange(desc(LP_strain_logFC)) %>%
  head(5), 
  
  EP_LP_strain %>%
  filter(Quadrant == "Up-Up") %>%
  arrange(desc(EP_strain_logFC)) %>%
  head(5),
  
  EP_LP_strain %>%
    filter(Quadrant == "Up-Down") %>%
    arrange(LP_strain_logFC) %>%
    head(3), 
  
  EP_LP_strain %>%
    filter(Quadrant == "Up-Down") %>%
    arrange(desc(EP_strain_logFC)) %>%
    head(3),
  
  EP_LP_strain %>%
    filter(Quadrant == "Down-Down") %>%
    arrange(LP_strain_logFC) %>%
    head(5), 
  
  EP_LP_strain %>%
    filter(Quadrant == "Down-Down") %>%
    arrange(EP_strain_logFC) %>%
    head(5),
  
  EP_LP_strain %>%
    filter(Quadrant == "Down-Up") %>%
    arrange(desc(LP_strain_logFC)) %>%
    head(3), 
  
  EP_LP_strain %>%
    filter(Quadrant == "Down-Up") %>%
    arrange(EP_strain_logFC) %>%
    head(3))


# Quadrant w/ labels 

EP_LP_strain <- EP_LP_strain %>%
  mutate(Label = ifelse(Pman_GeneID %in% Top_genes$Pman_GeneID, Pman_GeneID, NA))

ggplot(EP_LP_strain, aes(x = EP_strain_logFC, y = LP_strain_logFC, color = Quadrant)) +
  geom_point(alpha = 0.6, size = 3) +
  xlim(-6,6) +
  ylim(-6,6) +
  labs(y = NULL, x = NULL) +
  scale_color_manual(values = quadrant_colors) +
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
                  segment.alpha = 0.5  ) +
  labs(y = "Late Pregnancy\npopulation logFC", x = "Early Pregnancy\npopulation logFC")

# ggsave("Dream_Output/Hypoxia_Figs/AllgenesQuadrant.png", width = 10, height = 10, units = "in", dpi = 300)

# Highlight glycolysis genes ----

Gly_genes <- c("Aldoa",
               "Bpgm",
               "Eno1",
               "Eno2",
               "LOC102904208",
               "LOC102916082",
               "Gpi",
               "Hk1",
               "Hk2",
               "Hkdc1",
               "Pfkl",
               "Pfkm",
               "Pgam1",
               "Pgam2",
               "Pklr",
               "LOC102923285",
               "Pkm",
               "Tpi1",
               "LOC102928417",
               "Ldha")


EP_LP_strain = EP_LP_strain %>%
  mutate(Gly_flag = ifelse(Gene_ID %in% Gly_genes, "Glycolysis", "NA"))

ggplot(EP_LP_strain, aes(x = EP_strain_logFC, y = LP_strain_logFC, color = Quadrant)) +
  geom_point(data = EP_LP_strain %>% filter(Gly_flag == "NA"),
             color = "gray", 
             alpha = 0.6, size = 3) +
  geom_point(data = EP_LP_strain %>% filter(Gly_flag == "Glycolysis"), 
             color = "black", 
             alpha = 0.8, size = 3) + 
  xlim(-6,6) +
  ylim(-8,8) +
  labs(y = "Late Pregnancy\npopulation logFC", x = "Early Pregnancy\npopulation logFC") +
  scale_color_manual(values = quadrant_colors) +
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

#ggsave("Dream_Output/Hypoxia_Figs/GlycQuadrant.png", width = 12, height = 6, units = "in", dpi = 300)

Glycol_genes = EP_LP_strain %>%
              filter(Gly_flag == "Glycolysis")

EP_LP_strain <- EP_LP_strain %>%
  mutate(Glyc_Label = ifelse(Pman_GeneID %in% Glycol_genes$Pman_GeneID, Pman_GeneID, NA))


ggplot(EP_LP_strain, aes(x = EP_strain_logFC, y = LP_strain_logFC, color = Quadrant)) +
  geom_point(data = EP_LP_strain %>% filter(Gly_flag == "NA"),
             color = "gray", 
             alpha = 0.6, size = 3) +
  geom_point(data = EP_LP_strain %>% filter(Gly_flag == "Glycolysis"), 
             color = "black", 
             alpha = 0.8, size = 3) + 
  xlim(-6,6) +
  ylim(-6,6) +
  labs(y = "Late Pregnancy\npopulation logFC", x = "Early Pregnancy\npopulation logFC") +
  scale_color_manual(values = quadrant_colors) +
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
  geom_text_repel(aes(label = Glyc_Label), size = 5, na.rm = TRUE, color = "black")

# ggsave("Dream_Output/Hypoxia_Figs/GlycQuadrant.png", width = 10, height = 10, units = "in", dpi = 300)


# Combined Quadrant and Glycolysis Plot


















# Individual quadrants -----

# Up Up

# Get top 5 genes for each direction without duplicates
top_genes_x <- EP_LP_strain %>%
  filter(Quadrant == "Up-Up") %>%
  arrange(desc(EP_strain_logFC)) %>%
  head(20) %>%
  mutate(Direction = "Top 5 in logFC.x")

top_genes_y <- EP_LP_strain %>%
  filter(Quadrant == "Up-Up") %>%
  arrange(desc(LP_strain_logFC)) %>%
  head(20) %>%
  mutate(Direction = "Top 5 in logFC.y")

# Combine top genes for labeling, keeping unique genes only
top_genes <- bind_rows(top_genes_x, top_genes_y) %>%
  distinct(Pman_GeneID, .keep_all = TRUE)  # Ensure unique genes only

# Create the scatter plot
ggplot(EP_LP_strain, aes(x = EP_strain_logFC, y = LP_strain_logFC, color = Quadrant)) +
  geom_point(alpha = 0.6, size = 3) +
  scale_color_manual(values = quadrant_colors) +
  geom_hline(yintercept = upper_cutoff, linetype = "dashed") +
  geom_hline(yintercept = lower_cutoff, linetype = "dashed") +
  geom_vline(xintercept = upper_cutoff, linetype = "dashed") +
  geom_vline(xintercept = lower_cutoff, linetype = "dashed") +
  labs(
    x = "EP Strain Log10 Fold Change",
    y = "LP Strain Log10 Fold Change",
    color = "Expression Status"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top") +
  # Zoom into the Up-Up quadrant
  xlim(upper_cutoff, max(EP_LP_strain$EP_strain_logFC, na.rm = TRUE)) +
  ylim(upper_cutoff, max(EP_LP_strain$LP_strain_logFC, na.rm = TRUE)) +
  # Add labels for top genes with repulsion
  geom_text_repel(data = top_genes, aes(label = Pman_GeneID),
                  size = 4, color = "black",
                  box.padding = 0.5,
                  point.padding = 0.5)



# Down Down

# Get top 10 genes with lowest log fold changes in each direction
top_genes_x <- EP_LP_strain %>%
  filter(Quadrant == "Down-Down") %>%
  arrange(EP_strain_logFC) %>%
  head(6) %>%
  mutate(Direction = "Top 10 in logFC.x")

top_genes_y <- EP_LP_strain %>%
  filter(Quadrant == "Down-Down") %>%
  arrange(LP_strain_logFC) %>%
  head(6) %>%
  mutate(Direction = "Top 10 in logFC.y")

# Combine top genes for labeling, keeping unique genes only
top_genes <- bind_rows(top_genes_x, top_genes_y) %>%
  distinct(Pman_GeneID, .keep_all = TRUE)  # Ensure unique genes only

# Create the scatter plot
ggplot(EP_LP_strain, aes(x = EP_strain_logFC, y = LP_strain_logFC, color = Quadrant)) +
  geom_point(alpha = 0.6, size = 3) +
  scale_color_manual(values = quadrant_colors) +
  geom_hline(yintercept = lower_cutoff, linetype = "dashed") +
  geom_vline(xintercept = lower_cutoff, linetype = "dashed") +
  labs(
    x = "EP Strain Log10 Fold Change",
    y = "LP Strain Log10 Fold Change",
    color = "Expression Status"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top") +
  # Zoom into the Down-Down quadrant
  xlim(min(EP_LP_strain$EP_strain_logFC, na.rm = TRUE), lower_cutoff) +
  ylim(min(EP_LP_strain$LP_strain_logFC, na.rm = TRUE), lower_cutoff) +
  # Add labels for top 10 genes in both directions with repulsion
  geom_text_repel(data = top_genes, aes(label = Pman_GeneID),
                  size = 4, color = "black", 
                  box.padding = 0.5,
                  point.padding = 0.5)


