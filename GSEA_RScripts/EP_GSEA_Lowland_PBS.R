
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


# LOWLAND ONLY ------
# Ranked genes
rankings_low <- sign(EP_ISO_Strain$BW_O2_logFC)*(-log10(EP_ISO_Strain$BW_O2_P.Val)) # signed p values from spatial DGE as ranking
names(rankings_low) <- EP_ISO_Strain$Mus_GeneID # genes as names

head(rankings_low)

rankings_low <- sort(rankings_low, decreasing = TRUE) # sort genes by ranking
plot(rankings_low)

max(rankings_low)
min(rankings_low)

ggplot(data.frame(gene_symbol = names(rankings_low)[1:50], ranks = rankings_low[1:50]), aes(gene_symbol, ranks)) + 
  geom_point() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Run GSEA
GSEA_low <- fgsea(pathways = gene_sets, # List of gene sets to check
                  stats = rankings_low,
                  scoreType = 'std', # in this case we have both pos and neg rankings_low. if only pos or neg, set to 'pos', 'neg'
                  minSize = 10,
                  maxSize = 500,
                  nproc = 1) # for parallelisation

head(GSEA_low)

head(GSEA_low[order(pval), ])

sum(GSEA_low[, pval < 0.05])
sum(GSEA_low[, padj < 0.05])

# Let's filter for significant pathways first (optional)
GSEA_sig_low <- GSEA_low %>% 
  filter(padj < 0.05) %>%  # or pval < 0.01
  arrange(desc(NES))       # NES = normalized enrichment score

# Make sure pathway names are factors so they plot nicely
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
  filter(pathway == "HALLMARK_OXIDATIVE_PHOSPHORYLATION") %>%
  select(leadingEdge)



# PBS Outliers

PBS = read_xlsx("RNA_Seq_RawData/PBS_RDA_outliers.xlsx")

# OXPHOS
leading_edge_long <- GSEA_low %>%
  unnest(leadingEdge) %>%
  filter(pathway == "HALLMARK_OXIDATIVE_PHOSPHORYLATION") %>%
  select(leadingEdge)

overlap = intersect(PBS$gene_mus, leading_edge_long$leadingEdge)

# DNA REPAIR
leading_edge_long <- GSEA_low %>%
  unnest(leadingEdge) %>%
  filter(pathway == "HALLMARK_DNA_REPAIR") %>%
  select(leadingEdge)

overlap = intersect(PBS$gene_mus, leading_edge_long$leadingEdge)


# MYC
leading_edge_long <- GSEA_low %>%
  unnest(leadingEdge) %>%
  filter(pathway == "HALLMARK_MYC_TARGETS_V1") %>%
  select(leadingEdge)

overlap = intersect(PBS$gene_mus, leading_edge_long$leadingEdge)



