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

# WGCNA MODULE ANALYSIS ------
# Get available modules
available_modules <- unique(EP_ISO_Strain$EP_BW_Modules[!is.na(EP_ISO_Strain$EP_BW_Modules)])
cat("Available WGCNA modules:", paste(available_modules, collapse = ", "), "\n")


# Get significant pathways from lowland GSEA - fix the filter
target_pathways <- c("HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                     "HALLMARK_MYC_TARGETS_V1", 
                     "HALLMARK_DNA_REPAIR")

sig_pathways_low <- GSEA_low %>%
  filter(padj < 0.05 & pathway %in% target_pathways) %>%
  pull(pathway)

# Function to analyze module enrichment for each module
analyze_module_enrichment <- function(module_name) {
  # Get genes in this module
  module_genes <- EP_ISO_Strain$Mus_GeneID[EP_ISO_Strain$EP_BW_Modules == module_name & !is.na(EP_ISO_Strain$EP_BW_Modules)]
  
  # Create summary for this module
  module_fishers <- data.frame(Pathway = sig_pathways_low) %>%
    mutate(
      Module = module_name,
      # Get leading edge genes for each pathway
      leadingEdge = map(Pathway, ~{
        GSEA_low %>% filter(pathway == .x) %>% pull(leadingEdge) %>% .[[1]]
      }),
      
      # Counts for Fisher's test
      Total_genes = nrow(EP_ISO_Strain),
      Module_total = length(module_genes),
      LeadingEdge_size = map_int(leadingEdge, length),
      Module_in_LeadingEdge = map_int(leadingEdge, ~{
        sum(.x %in% module_genes)
      }),
      NotModule_in_LeadingEdge = LeadingEdge_size - Module_in_LeadingEdge
    ) %>%
    # Run Fisher's tests  
    mutate(
      fisher_table = pmap(list(Module_in_LeadingEdge, NotModule_in_LeadingEdge, Module_total, Total_genes),
                          ~ matrix(c(..1, ..2, ..3 - ..1, ..4 - ..2 - ..3 + ..1), nrow = 2)),
      Fisher_pval = map_dbl(fisher_table, ~ fisher.test(.)$p.value),
      Fisher_FDR = p.adjust(Fisher_pval, method = "fdr"),
      Fisher_SIG = if_else(Fisher_FDR < 0.05, "SIG", "NS"),
      Odds_Ratio = map_dbl(fisher_table, ~ fisher.test(.)$estimate)
    ) %>%
    select(Module, Pathway, LeadingEdge_size, Module_total, Module_in_LeadingEdge, 
           Fisher_pval, Fisher_FDR, Fisher_SIG, Odds_Ratio) %>%
    arrange(Fisher_pval)
  
  return(module_fishers)
}

# Analyze all modules
all_module_results <- map_dfr(available_modules, analyze_module_enrichment)

# View results
print("Module enrichment results:")
print(all_module_results)

# Summary by module
module_summary <- all_module_results %>%
  group_by(Module) %>%
  summarise(
    Total_pathways = n(),
    Significant_pathways = sum(Fisher_SIG == "SIG"),
    Min_pvalue = min(Fisher_pval),
    .groups = "drop"
  ) %>%
  arrange(desc(Significant_pathways), Min_pvalue)

print("Summary by module:")
print(module_summary)

# Plot significant enrichments
sig_enrichments <- all_module_results %>%
  filter(Fisher_SIG == "SIG") %>%
  arrange(Fisher_pval)

if(nrow(sig_enrichments) > 0) {
  ggplot(sig_enrichments, aes(x = -log10(Fisher_FDR), y = reorder(paste(Module, Pathway, sep = " - "), -log10(Fisher_FDR)), 
                              fill = Module, size = Odds_Ratio)) +
    geom_point(shape = 21, alpha = 0.7) +
    labs(x = "-log10(Global Fisher FDR)",
         y = "Module - Pathway",
         title = "Significant WGCNA Module Enrichments (Global FDR < 0.05)",
         size = "Odds Ratio",
         fill = "WGCNA Module") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8))
} else {
  print("No significant module enrichments found after global FDR correction.")
}

