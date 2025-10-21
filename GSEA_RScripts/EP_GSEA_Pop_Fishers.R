

library(fgsea)
library(tidyverse)
library(RColorBrewer)
library(msigdbr)
library(readxl)

# READ IN DATA -------
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
rankings_low <- sign(EP_ISO_Strain$BW_O2_logFC)*(-log10(EP_ISO_Strain$BW_O2_adj_P.Val)) # signed p values from spatial DGE as ranking
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


# Get significant pathways from lowland GSEA
sig_pathways_low <- GSEA_low %>% filter(padj < 0.05) %>% pull(pathway)

# Create summary for lowland pathways
lowland_fishers <- data.frame(Pathway = sig_pathways_low) %>%
  mutate(
    # Get leading edge genes for each pathway
    leadingEdge = map(Pathway, ~{
      GSEA_low %>% filter(pathway == .x) %>% pull(leadingEdge) %>% .[[1]]
    }),
    
    # Counts for Fisher's test
    Total_genes = 20733,  # adjust to your total
    PopDE_total = 7750,   # your pop DE count
    LeadingEdge_size = map_int(leadingEdge, length),
    PopDE_in_LeadingEdge = map_int(leadingEdge, ~{
      sum(.x %in% EP_ISO_Strain$Mus_GeneID[EP_ISO_Strain$strain_SIG == "SIG"])
    }),
    NotPopDE_in_LeadingEdge = LeadingEdge_size - PopDE_in_LeadingEdge
  ) %>%
  # Run Fisher's tests  
  mutate(
    fisher_table = pmap(list(PopDE_in_LeadingEdge, NotPopDE_in_LeadingEdge, PopDE_total, Total_genes),
                        ~ matrix(c(..1, ..2, ..3 - ..1, ..4 - ..2 - ..3 + ..1), nrow = 2)),
    Fisher_pval = map_dbl(fisher_table, ~ fisher.test(.)$p.value),
    Fisher_FDR = p.adjust(Fisher_pval, method = "fdr"),
    Fisher_SIG = if_else(Fisher_FDR < 0.05, "SIG", "NS"),
    Odds_Ratio = map_dbl(fisher_table, ~ fisher.test(.)$estimate)
  ) %>%
  select(Pathway, LeadingEdge_size, PopDE_in_LeadingEdge, Fisher_pval, Fisher_FDR, Fisher_SIG, Odds_Ratio) %>%
  arrange(Fisher_pval)




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

# Highland Fisher's tests
sig_pathways_high <- GSEA_high %>% filter(padj < 0.05) %>% pull(pathway)

highland_fishers <- data.frame(Pathway = sig_pathways_high) %>%
  mutate(
    # Get leading edge genes for each pathway
    leadingEdge = map(Pathway, ~{
      GSEA_high %>% filter(pathway == .x) %>% pull(leadingEdge) %>% .[[1]]
    }),
    
    # Counts for Fisher's test
    Total_genes = 20733,  # adjust to your total
    PopDE_total = 7750,   # your pop DE count
    LeadingEdge_size = map_int(leadingEdge, length),
    PopDE_in_LeadingEdge = map_int(leadingEdge, ~{
      sum(.x %in% EP_ISO_Strain$Mus_GeneID[EP_ISO_Strain$strain_SIG == "SIG"])
    }),
    NotPopDE_in_LeadingEdge = LeadingEdge_size - PopDE_in_LeadingEdge
  ) %>%
  # Run Fisher's tests  
  mutate(
    fisher_table = pmap(list(PopDE_in_LeadingEdge, NotPopDE_in_LeadingEdge, PopDE_total, Total_genes),
                        ~ matrix(c(..1, ..2, ..3 - ..1, ..4 - ..2 - ..3 + ..1), nrow = 2)),
    Fisher_pval = map_dbl(fisher_table, ~ fisher.test(.)$p.value),
    Fisher_FDR = p.adjust(Fisher_pval, method = "fdr"),
    Fisher_SIG = if_else(Fisher_FDR < 0.05, "SIG", "NS"),
    Odds_Ratio = map_dbl(fisher_table, ~ fisher.test(.)$estimate)
  ) %>%
  select(Pathway, LeadingEdge_size, PopDE_in_LeadingEdge, Fisher_pval, Fisher_FDR, Fisher_SIG, Odds_Ratio) %>%
  arrange(Fisher_pval)


# Look at leading edge sizes
lowland_fishers %>% 
  arrange(desc(LeadingEdge_size)) %>%
  select(Pathway, LeadingEdge_size, PopDE_in_LeadingEdge, Odds_Ratio, Fisher_pval)

highland_fishers %>% 
  arrange(desc(LeadingEdge_size)) %>%
  select(Pathway, LeadingEdge_size, PopDE_in_LeadingEdge, Odds_Ratio, Fisher_pval)

