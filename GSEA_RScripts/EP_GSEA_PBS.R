
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

cat(paste("Loaded", length(gene_sets), "gene sets\n"))

# Read in data
EP_ISO_Strain <- read_xlsx("RNA_Seq_Output/EP_ISO_Ortho_Summary.xlsx")


# HIGHLAND ONLY ------
# Ranked genes
rankings_high <- sign(EP_ISO_Strain$ME_O2_logFC)*(-log10(EP_ISO_Strain$ME_O2_adj_P.Val)) # signed p values from spatial DGE as ranking
names(rankings_high) <- EP_ISO_Strain$Mus_GeneID # genes as names

head(rankings_high)

rankings_high <- sort(rankings_high, decreasing = TRUE) # sort genes by ranking
plot(rankings_high)

max(rankings_high)
min(rankings_high)

ggplot(data.frame(gene_symbol = names(rankings_high)[1:50], ranks = rankings_high[1:50]), aes(gene_symbol, ranks)) + 
  geom_point() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Run GSEA
GSEA_high <- fgsea(pathways = gene_sets, # List of gene sets to check
                   stats = rankings_high,
                   scoreType = 'std', # in this case we have both pos and neg rankings_high. if only pos or neg, set to 'pos', 'neg'
                   minSize = 10,
                   maxSize = 500,
                   nproc = 1) # for parallelisation

head(GSEA_high)

head(GSEA_high[order(pval), ])

sum(GSEA_high[, pval < 0.05])
sum(GSEA_high[, padj < 0.05])

# Let's filter for significant pathways first (optional)
GSEA_sig_high <- GSEA_high %>% 
  filter(padj < 0.05) %>%  # or pval < 0.01
  arrange(desc(NES))       # NES = normalized enrichment score

# Make sure pathway names are factors so they plot nicely
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



# PBS Outliers

PBS = read_xlsx("RNA_Seq_RawData/PBS_RDA_outliers.xlsx")

# OXPHOS
leading_edge_long <- GSEA_high %>%
  unnest(leadingEdge) %>%
  filter(pathway == "HALLMARK_OXIDATIVE_PHOSPHORYLATION") %>%
  select(leadingEdge)

overlap = intersect(PBS$gene_mus, leading_edge_long$leadingEdge)

# DNA REPAIR
leading_edge_long <- GSEA_high %>%
  unnest(leadingEdge) %>%
  filter(pathway == "HALLMARK_DNA_REPAIR") %>%
  select(leadingEdge)

overlap = intersect(PBS$gene_mus, leading_edge_long$leadingEdge)


# MYC
leading_edge_long <- GSEA_high %>%
  unnest(leadingEdge) %>%
  filter(pathway == "HALLMARK_MYC_TARGETS_V1") %>%
  select(leadingEdge)

overlap = intersect(PBS$gene_mus, leading_edge_long$leadingEdge)


# HEME_METABOLISM
leading_edge_long <- GSEA_high %>%
  unnest(leadingEdge) %>%
  filter(pathway == "HALLMARK_HEME_METABOLISM") %>%
  select(leadingEdge)

overlap = intersect(PBS$gene_mus, leading_edge_long$leadingEdge)


# E2F
leading_edge_long <- GSEA_high %>%
  unnest(leadingEdge) %>%
  filter(pathway == "HALLMARK_E2F_TARGETS") %>%
  select(leadingEdge)

overlap = intersect(PBS$gene_mus, leading_edge_long$leadingEdge)




