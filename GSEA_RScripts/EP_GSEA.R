
library(fgsea)
library(tidyverse)
library(RColorBrewer)
library(msigdbr)
library(readxl)

# Download gene sets
gene_sets_df <- msigdbr(species = 'mouse', category = 'H')

# Convert to named list format required by fgsea
gene_sets <- gene_sets_df %>%
  split(x = .$gene_symbol, f = .$gs_name)

# Read in data
EP_ISO_Strain <- read_xlsx("RNA_Seq_Output/EP_ISO_Ortho_Summary.xlsx")


# POPULATION -------
# Ranked genes
rankings <- sign(EP_ISO_Strain$strain_logFC)*(-log10(EP_ISO_Strain$strain_adj_P.Val)) 
names(rankings) <- EP_ISO_Strain$Mus_GeneID 

head(rankings)
rankings <- sort(rankings, decreasing = TRUE) 
plot(rankings)

max(rankings)
min(rankings)

plot(rankings)

ggplot(data.frame(gene_symbol = names(rankings)[1:50], ranks = rankings[1:50]), aes(gene_symbol, ranks)) + 
  geom_point() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# Run GSEA
GSEAres <- fgsea(pathways = gene_sets, 
                 stats = rankings,
                 scoreType = 'std', 
                 minSize = 10,
                 maxSize = 500,
                 nproc = 1) 

head(GSEAres)

sum(GSEAres[, pval < 0.05])
sum(GSEAres[, padj < 0.05])


# HYPOXIA -------
# Ranked genes
rankings <- sign(EP_ISO_Strain$O2_logFC)*(-log10(EP_ISO_Strain$O2_adj_P.Val)) 
names(rankings) <- EP_ISO_Strain$Mus_GeneID 

head(rankings)

rankings <- sort(rankings, decreasing = TRUE) 
plot(rankings)

max(rankings)
min(rankings)

ggplot(data.frame(gene_symbol = names(rankings)[1:50], ranks = rankings[1:50]), aes(gene_symbol, ranks)) + 
  geom_point() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Run GSEA
GSEAres <- fgsea(pathways = gene_sets, 
                 stats = rankings,
                 scoreType = 'std', 
                 minSize = 10,
                 maxSize = 500,
                 nproc = 1) 

head(GSEAres)

head(GSEAres[order(pval), ])

sum(GSEAres[, pval < 0.05])
sum(GSEAres[, padj < 0.05])


GSEA_sig <- GSEAres %>% 
  filter(padj < 0.05) %>%  
  arrange(desc(NES))       


GSEA_sig$pathway <- factor(GSEA_sig$pathway, levels = GSEA_sig$pathway)

# Dot plot
ggplot(GSEA_sig, aes(x = NES, y = pathway, size = -log10(padj), color = NES)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(x = "Normalized Enrichment Score (NES)",
       y = "Pathway",
       size = "-log10(adj p-value)",
       color = "NES") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))


# LOWLAND ONLY ------
# Ranked genes
rankings_low <- sign(EP_ISO_Strain$BW_O2_logFC)*(-log10(EP_ISO_Strain$BW_O2_adj_P.Val)) 
names(rankings_low) <- EP_ISO_Strain$Mus_GeneID 

head(rankings_low)

rankings_low <- sort(rankings_low, decreasing = TRUE) 
plot(rankings_low)

max(rankings_low)
min(rankings_low)

ggplot(data.frame(gene_symbol = names(rankings_low)[1:50], ranks = rankings_low[1:50]), aes(gene_symbol, ranks)) + 
  geom_point() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Run GSEA
GSEA_low <- fgsea(pathways = gene_sets, 
                  stats = rankings_low,
                  scoreType = 'std', 
                  minSize = 10,
                  maxSize = 500,
                  nproc = 1) 

head(GSEA_low)

head(GSEA_low[order(pval), ])

sum(GSEA_low[, pval < 0.05])
sum(GSEA_low[, padj < 0.05])


GSEA_sig_low <- GSEA_low %>% 
  filter(padj < 0.05) %>%  
  arrange(desc(NES))       



GSEA_sig_low$pathway <- factor(GSEA_sig_low$pathway, levels = GSEA_sig_low$pathway)

# Dot plot
ggplot(GSEA_sig_low, aes(x = NES, y = pathway, size = -log10(padj), color = NES)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(x = "Normalized Enrichment Score (NES)",
       y = "Pathway",
       size = "-log10(adj p-value)",
       color = "NES") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

leading_edge_long <- GSEA_low %>%
  unnest(leadingEdge) %>%
  filter(pathway == "HALLMARK_OXIDATIVE_PHOSPHORYLATION")


# HIGHLAND ONLY ------
# Ranked genes
rankings_high <- sign(EP_ISO_Strain$ME_O2_logFC)*(-log10(EP_ISO_Strain$ME_O2_adj_P.Val)) 
names(rankings_high) <- EP_ISO_Strain$Mus_GeneID 

head(rankings_high)

rankings_high <- sort(rankings_high, decreasing = TRUE) 
plot(rankings_high)

max(rankings_high)
min(rankings_high)

ggplot(data.frame(gene_symbol = names(rankings_high)[1:50], ranks = rankings_high[1:50]), aes(gene_symbol, ranks)) + 
  geom_point() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Run GSEA
GSEA_high <- fgsea(pathways = gene_sets, 
                   stats = rankings_high,
                   scoreType = 'std', 
                   maxSize = 500,
                   nproc = 1)

head(GSEA_high)

head(GSEA_high[order(pval), ])

sum(GSEA_high[, pval < 0.05])
sum(GSEA_high[, padj < 0.05])


GSEA_sig_high <- GSEA_high %>% 
  filter(padj < 0.05) %>%  
  arrange(desc(NES))       


GSEA_sig_high$pathway <- factor(GSEA_sig_high$pathway, levels = GSEA_sig_high$pathway)

# Dot plot
ggplot(GSEA_sig_high, aes(x = NES, y = pathway, size = -log10(padj), color = NES)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(x = "Normalized Enrichment Score (NES)",
       y = "Pathway",
       size = "-log10(adj p-value)",
       color = "NES") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))


plotEnrichment(gene_sets[['HALLMARK_DNA_REPAIR']],
               rankings_high) + 
  labs(title = 'Hallmark: DNA Repair') + 
  theme_classic() +
  scale_x_continuous('Rank', breaks = seq(0, 32000, 5000)) +
  scale_y_continuous('Enrichment score (ES)') +
  geom_line(col = 'purple', linewidth = 2)

plotEnrichment(gene_sets[['HALLMARK_HEME_METABOLISM']],
               rankings_high) + 
  labs(title = 'Hallmark: Heme Metabolism') + 
  theme_classic() +
  scale_x_continuous('Rank', breaks = seq(0, 32000, 5000)) +
  scale_y_continuous('Enrichment score (ES)') +
  geom_line(col = 'purple', linewidth = 2) 

plotEnrichment(gene_sets[['HALLMARK_TNFA_SIGNALING_VIA_NFKB']],
               rankings_high) + 
  labs(title = 'Hallmark: TNFA Signalling via NFKB') + 
  theme_classic() +
  scale_x_continuous('Rank', breaks = seq(0, 32000, 5000)) +
  scale_y_continuous('Enrichment score (ES)') +
  geom_line(col = 'purple', linewidth = 2) 

plotEnrichment(gene_sets[['HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION']],
               rankings_high) + 
  labs(title = 'Hallmark: Epithelial Mesenchymal Transition') + 
  theme_classic() +
  scale_x_continuous('Rank', breaks = seq(0, 32000, 5000)) +
  scale_y_continuous('Enrichment score (ES)') +
  geom_line(col = 'purple', linewidth = 2) 


# IXN ------
# Ranked genes
rankings <- sign(EP_ISO_Strain$IXN_logFC)*(-log10(EP_ISO_Strain$IXN_adj_P.Val)) 
names(rankings) <- EP_ISO_Strain$Mus_GeneID 

head(rankings)
rankings <- sort(rankings, decreasing = TRUE) 
plot(rankings)

max(rankings)
min(rankings)

plot(rankings)

ggplot(data.frame(gene_symbol = names(rankings)[1:50], ranks = rankings[1:50]), aes(gene_symbol, ranks)) + 
  geom_point() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Run GSEA
GSEAres <- fgsea(pathways = gene_sets, 
                 stats = rankings,
                 scoreType = 'std', 
                 minSize = 10,
                 maxSize = 500,
                 nproc = 1) 

head(GSEAres)

sum(GSEAres[, pval < 0.05])
sum(GSEAres[, padj < 0.05])

GSEA_sig <- GSEAres %>% 
  filter(padj < 0.05) %>%  
  arrange(desc(NES))       


GSEA_sig$pathway <- factor(GSEA_sig$pathway, levels = GSEA_sig$pathway)

# Dot plot
ggplot(GSEA_sig, aes(x = NES, y = pathway, size = -log10(padj), color = NES)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(x = "Normalized Enrichment Score (NES)",
       y = "Pathway",
       size = "-log10(adj p-value)",
       color = "NES") +
  xlim(-3,3) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

# ggsave("Dream_Output/Plots/GSEA_Plots/Ixn_GSEA.png", width = 6, height = 7, units = "in", dpi = 300)

plotEnrichment(gene_sets[['HALLMARK_OXIDATIVE_PHOSPHORYLATION']],
               rankings) + 
  labs(title = 'Hallmark: Oxidative Phosphorylation') + 
  theme_classic() +
  scale_x_continuous('Rank', breaks = seq(0, 32000, 5000)) +
  scale_y_continuous('Enrichment score (ES)') +
  geom_line(col = 'purple', linewidth = 2) 

plotEnrichment(gene_sets[['HALLMARK_DNA_REPAIR']],
               rankings) + 
  labs(title = 'Hallmark: DNA Repair') + 
  theme_classic() +
  scale_x_continuous('Rank', breaks = seq(0, 32000, 5000)) +
  scale_y_continuous('Enrichment score (ES)') +
  geom_line(col = 'purple', linewidth = 2) 


leading_edge_long <- GSEAres %>%
  unnest(leadingEdge) %>%
  filter(pathway == "HALLMARK_OXIDATIVE_PHOSPHORYLATION") %>%
  select(leadingEdge)

leading_edge_long <- GSEAres %>%
  unnest(leadingEdge) %>%
  filter(pathway == "HALLMARK_DNA_REPAIR")

leading_edge_long <- GSEAres %>%
  unnest(leadingEdge) %>%
  filter(pathway == "HALLMARK_PROTEIN_SECRETION") %>%
  select(leadingEdge)

leading_edge_long <- GSEAres %>%
  unnest(leadingEdge) %>%
  filter(pathway == "HALLMARK_PEROXISOME") %>%
  select(leadingEdge)


