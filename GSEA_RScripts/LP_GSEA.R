
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
LP_ISO_Strain <- read_xlsx("RNA_Seq_Output/LP_ISO_Ortho_Summary.xlsx")

# EARLY PREG -------

# POPULATION -------
# Ranked genes
rankings <- sign(LP_ISO_Strain$strain_logFC)*(-log10(LP_ISO_Strain$strain_adj_P.Val)) # signed p values from spatial DGE as ranking
names(rankings) <- LP_ISO_Strain$Mus_GeneID # genes as names

head(rankings)
rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
plot(rankings)

max(rankings)
min(rankings)

plot(rankings)

# Some genes have such low p values that the signed pval is +- inf, we need to change it to the maximum * constant to avoid problems with fgsea
max_ranking <- max(rankings[is.finite(rankings)])
min_ranking <- min(rankings[is.finite(rankings)])
rankings <- replace(rankings, rankings > max_ranking, max_ranking * 10)
rankings <- replace(rankings, rankings < min_ranking, min_ranking * 10)
rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking

plot(rankings)

ggplot(data.frame(gene_symbol = names(rankings)[1:50], ranks = rankings[1:50]), aes(gene_symbol, ranks)) + 
  geom_point() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# Run GSEA
GSEAres <- fgsea(pathways = gene_sets, # List of gene sets to check
                 stats = rankings,
                 scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                 minSize = 10,
                 maxSize = 500,
                 nproc = 1) # for parallelisation

head(GSEAres)

sum(GSEAres[, pval < 0.01])
sum(GSEAres[, padj < 0.01])


# HYPOXIA -------
# Ranked genes
rankings <- sign(LP_ISO_Strain$O2_logFC)*(-log10(LP_ISO_Strain$O2_adj_P.Val)) # signed p values from spatial DGE as ranking
names(rankings) <- LP_ISO_Strain$Mus_GeneID # genes as names

head(rankings)

rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
plot(rankings)

max(rankings)
min(rankings)

ggplot(data.frame(gene_symbol = names(rankings)[1:50], ranks = rankings[1:50]), aes(gene_symbol, ranks)) + 
  geom_point() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Run GSEA
GSEAres <- fgsea(pathways = gene_sets, # List of gene sets to check
                 stats = rankings,
                 scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                 minSize = 10,
                 maxSize = 500,
                 nproc = 1) # for parallelisation

head(GSEAres)

head(GSEAres[order(pval), ])

sum(GSEAres[, pval < 0.01])
sum(GSEAres[, padj < 0.01])

# Let's filter for significant pathways first (optional)
GSEA_sig <- GSEAres %>% 
  filter(padj < 0.05) %>%  # or pval < 0.01
  arrange(desc(NES))       # NES = normalized enrichment score

# Make sure pathway names are factors so they plot nicely
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
rankings <- sign(LP_ISO_Strain$BW_O2_logFC)*(-log10(LP_ISO_Strain$BW_O2_adj_P.Val)) # signed p values from spatial DGE as ranking
names(rankings) <- LP_ISO_Strain$Mus_GeneID # genes as names

head(rankings)

rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
plot(rankings)

max(rankings)
min(rankings)

ggplot(data.frame(gene_symbol = names(rankings)[1:50], ranks = rankings[1:50]), aes(gene_symbol, ranks)) + 
  geom_point() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Run GSEA
GSEAres <- fgsea(pathways = gene_sets, # List of gene sets to check
                 stats = rankings,
                 scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                 minSize = 10,
                 maxSize = 500,
                 nproc = 1) # for parallelisation

head(GSEAres)

head(GSEAres[order(pval), ])

sum(GSEAres[, pval < 0.05])
sum(GSEAres[, padj < 0.05])

# Let's filter for significant pathways first (optional)
GSEA_sig <- GSEAres %>% 
  filter(padj < 0.05) %>%  # or pval < 0.01
  arrange(desc(NES))       # NES = normalized enrichment score

# Make sure pathway names are factors so they plot nicely
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


# HIGHLAND ONLY ------
# Ranked genes
rankings <- sign(LP_ISO_Strain$ME_O2_logFC)*(-log10(LP_ISO_Strain$ME_O2_adj_P.Val)) # signed p values from spatial DGE as ranking
names(rankings) <- LP_ISO_Strain$Mus_GeneID # genes as names

head(rankings)

rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
plot(rankings)

max(rankings)
min(rankings)

ggplot(data.frame(gene_symbol = names(rankings)[1:50], ranks = rankings[1:50]), aes(gene_symbol, ranks)) + 
  geom_point() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Run GSEA
GSEAres <- fgsea(pathways = gene_sets, # List of gene sets to check
                 stats = rankings,
                 scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                 minSize = 10,
                 maxSize = 500,
                 nproc = 1) # for parallelisation

head(GSEAres)

head(GSEAres[order(pval), ])

sum(GSEAres[, pval < 0.05])
sum(GSEAres[, padj < 0.05])

# Let's filter for significant pathways first (optional)
GSEA_sig <- GSEAres %>% 
  filter(padj < 0.05) %>%  # or pval < 0.01
  arrange(desc(NES))       # NES = normalized enrichment score

# Make sure pathway names are factors so they plot nicely
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




# HYPOXIA ENSEMBL
rankings <- sign(LP_ISO_Strain$O2_logFC)*(-log10(LP_ISO_Strain$O2_adj_P.Val)) # signed p values from spatial DGE as ranking
names(rankings) <- LP_ISO_Strain$Ensembl_ID # genes as names 

# StLP 2: turn into dataframe to handle duplicates
rankings_df <- data.frame(
  gene  = names(rankings),
  score = as.numeric(rankings),
  stringsAsFactors = FALSE
)

# StLP 3: collapse duplicates (keLP max score per gene)
rankings_df_unique <- rankings_df %>%
  group_by(gene) %>%
  slice_max(order_by = abs(score), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(score))

rankings_df_unique <- rankings_df_unique %>%
  filter(gene != "None")


# StLP 4: back to named vector
rankings_unique <- rankings_df_unique$score
names(rankings_unique) <- rankings_df_unique$gene

# StLP 5: sanity checks
head(rankings_unique)
plot(rankings_unique, main = "Ranked gene statistics (unique IDs)",
     ylab = "Score", xlab = "Genes")
max(rankings_unique, na.rm = TRUE)
min(rankings_unique, na.rm = TRUE)

# Run GSEA
GSEAres <- fgsea(pathways = gene_sets, # List of gene sets to check
                 stats = rankings_unique,
                 scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                 minSize = 10,
                 maxSize = 500,
                 nproc = 1) # for parallelisation

head(GSEAres)

head(GSEAres[order(pval), ])

sum(GSEAres[, pval < 0.01])
sum(GSEAres[, padj < 0.01])

# Let's filter for significant pathways first (optional)
GSEA_sig <- GSEAres %>% 
  filter(padj < 0.05) %>%  # or pval < 0.01
  arrange(desc(NES))       # NES = normalized enrichment score

# Make sure pathway names are factors so they plot nicely
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




# IXN
# Ranked genes
rankings <- sign(LP_ISO_Strain$IXN_logFC)*(-log10(LP_ISO_Strain$IXN_adj_P.Val)) # signed p values from spatial DGE as ranking
names(rankings) <- LP_ISO_Strain$Pman_GeneID # genes as names

head(rankings)
rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
plot(rankings)

max(rankings)
min(rankings)

plot(rankings)

ggplot(data.frame(gene_symbol = names(rankings)[1:50], ranks = rankings[1:50]), aes(gene_symbol, ranks)) + 
  geom_point() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Run GSEA
GSEAres <- fgsea(pathways = gene_sets, # List of gene sets to check
                 stats = rankings,
                 scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                 minSize = 10,
                 maxSize = 500,
                 nproc = 1) # for parallelisation

head(GSEAres)

sum(GSEAres[, pval < 0.01])
sum(GSEAres[, padj < 0.01])

# Let's filter for significant pathways first (optional)
GSEA_sig <- GSEAres %>% 
  filter(padj < 0.05) %>%  # or pval < 0.01
  arrange(desc(NES))       # NES = normalized enrichment score

# Make sure pathway names are factors so they plot nicely
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







