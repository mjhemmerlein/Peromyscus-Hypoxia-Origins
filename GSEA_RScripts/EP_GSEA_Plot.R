
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

leading_edge_long <- GSEA_low %>%
  unnest(leadingEdge) %>%
  filter(pathway == "HALLMARK_OXIDATIVE_PHOSPHORYLATION") %>%
  select(leadingEdge)


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

leading_edge_long <- GSEA_high %>%
  unnest(leadingEdge) %>%
  filter(pathway == "HALLMARK_OXIDATIVE_PHOSPHORYLATION") %>%
  select(leadingEdge)

leading_edge_long <- GSEA_high %>%
  unnest(leadingEdge) %>%
  filter(pathway == "HALLMARK_DNA_REPAIR") %>%
  select(leadingEdge)

leading_edge_long <- GSEA_high %>%
  unnest(leadingEdge) %>%
  filter(pathway == "HALLMARK_MYC_TARGETS_V1") %>%
  select(leadingEdge)


# FACETED POPULATION GSEA PLOT ------
# Facet by population specificity, rank by NES within groups

# Setup data
GSEA_high$condition <- "Highland"
GSEA_low$condition  <- "Lowland"
all_data <- bind_rows(GSEA_high, GSEA_low)

# Define significant pathways
sig_high <- GSEA_high %>% filter(padj < 0.05)
sig_low  <- GSEA_low  %>% filter(padj < 0.05)
shared_sig <- intersect(sig_high$pathway, sig_low$pathway)
highland_only <- setdiff(sig_high$pathway, shared_sig)
lowland_only <- setdiff(sig_low$pathway, shared_sig)

# Create the grouped dataset
plot_data <- all_data %>%
  # Classify pathways into groups
  mutate(
    pathway_group = case_when(
      pathway %in% shared_sig ~ "Shared",
      pathway %in% highland_only ~ "Highland-only", 
      pathway %in% lowland_only ~ "Lowland-only",
      TRUE ~ "Non-significant"
    )
  ) %>%
  # Filter to only significant pathways
  filter(pathway_group != "Non-significant") %>%
  # For ordering within groups:
  group_by(pathway_group, pathway) %>%
  mutate(
    # For shared: use absolute difference between populations (biggest differences FIRST)
    # For population-specific: use the NES from that population
    rank_value = case_when(
      pathway_group == "Shared" ~ abs(diff(NES)),  # Absolute difference between Highland and Lowland
      pathway_group == "Highland-only" ~ NES[condition == "Highland"],
      pathway_group == "Lowland-only" ~ NES[condition == "Lowland"]
    )
  ) %>%
  ungroup() %>%
  # Create ordering within each group - LARGEST differences first for shared
  arrange(pathway_group, ifelse(pathway_group == "Shared", -rank_value, rank_value)) %>%
  mutate(
    # Order factor levels for plotting - SHARED IN MIDDLE
    pathway_group = factor(pathway_group, levels = c("Lowland-only", "Shared", "Highland-only")),
    pathway_ordered = factor(pathway, levels = unique(pathway))
  )

# Create the plot with p-value-based point sizing
  
ggplot(plot_data, aes(x = NES, y = pathway_ordered, color = condition)) +
  # Dumbbell lines for shared pathways
  geom_line(
    data = plot_data %>% filter(pathway_group == "Shared"),
    aes(group = pathway_ordered),
    color = "gray60",
    linewidth = 0.8
  ) +
  # Points - only show significant ones, sized by p-value
  geom_point(
    data = plot_data %>% filter(padj < 0.05),
    aes(size = -log10(padj))) +
  geom_vline(xintercept = 0, linetype = "solid", color = "gray50") +
  
  # Facet by pathway group
  facet_grid(pathway_group ~ ., scales = "free_y", space = "free_y") +
  
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = "Pathway",
    color = "Population",
    size = "-log10(adj p-value)"
  ) +
  scale_x_continuous(limits = c(-3, 3)) +
  scale_color_manual(values = c("Highland" = "#99c6cc", "Lowland" = "#cba34e")) +
  scale_size_continuous(range = c(1.5, 4.5), guide = guide_legend(override.aes = list(color = "gray50"))) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    strip.text.y = element_text(size = 11, face = "bold", angle = 0),
    strip.background = element_rect(fill = "gray95", color = "white"),
    legend.position = "right",
    legend.box = "horizontal",  # Put legends side by side
    panel.spacing = unit(0.8, "lines"),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "gray40")
  )

ggsave("Plots/GSEA_Plots/SharedPathway_GSEA.pdf", width = 12, height = 8, units = "in", dpi = 300)

# PATHWAY PLOTS ------------

library(patchwork)

# Create both plots first
plot1 <- plotEnrichment(gene_sets[['HALLMARK_OXIDATIVE_PHOSPHORYLATION']],
                        rankings_low) + 
  labs(title = 'Population 1', subtitle = 'Hallmark: OXPHOS') + 
  theme_classic() +
  scale_x_continuous('Rank', breaks = seq(0, 32000, 5000)) +
  geom_line(col = 'purple', linewidth = 2)

plot2 <- plotEnrichment(gene_sets[['HALLMARK_OXIDATIVE_PHOSPHORYLATION']],
                        rankings_high) + 
  labs(title = 'Population 2', subtitle = 'Hallmark: OXPHOS') + 
  theme_classic() +
  scale_x_continuous('Rank', breaks = seq(0, 32000, 5000)) +
  geom_line(col = 'orange', linewidth = 2)

# Find the y-axis range that encompasses both plots
y_min <- min(c(layer_data(plot1)$y, layer_data(plot2)$y))
y_max <- max(c(layer_data(plot1)$y, layer_data(plot2)$y))

# Apply same y-axis limits to both
plot1_fixed <- plot1 + scale_y_continuous('Enrichment score (ES)', 
                                          limits = c(y_min, y_max))
plot2_fixed <- plot2 + scale_y_continuous('Enrichment score (ES)', 
                                          limits = c(y_min, y_max))

# Combine with shared y-axis
plot1_fixed | plot2_fixed



# Calculate enrichment scores for both populations
geneset <- gene_sets[['HALLMARK_OXIDATIVE_PHOSPHORYLATION']]

# Get enrichment data for both
enrich_data_low <- plotEnrichment(geneset, rankings_low)$data
enrich_data_high <- plotEnrichment(geneset, rankings_high)$data

# Add population labels
enrich_data_low$population <- "Low"
enrich_data_high$population <- "High"

# Combine data
combined_data <- rbind(enrich_data_low, enrich_data_high)

# Plot with taller gene markers and expanded y-axis
ggplot(combined_data, aes(x = x, y = y, color = population)) +
  # Add horizontal reference lines
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.7) +
  geom_hline(yintercept = max(combined_data$y), linetype = "dotted", color = "gray70", alpha = 0.5) +
  geom_hline(yintercept = min(combined_data$y), linetype = "dotted", color = "gray70", alpha = 0.5) +
  
  # Main enrichment curves
  geom_line(linewidth = 2) +
  
  # Add taller gene tick marks at bottom
  geom_rug(data = subset(combined_data, population == "Low"), 
           aes(x = x), 
           sides = "b", 
           color = "purple", 
           alpha = 0.6,
           length = unit(0.05, "npc")) +  # Made taller (was 0.02)
  geom_rug(data = subset(combined_data, population == "High"), 
           aes(x = x), 
           sides = "b", 
           color = "orange", 
           alpha = 0.6,
           length = unit(0.05, "npc")) +  # Made taller (was 0.02)
  
  scale_color_manual(values = c("Low" = "purple", 
                                "High" = "orange")) +
  labs(title = 'Hallmark: Oxidative Phosphorylation - Comparison',
       x = 'Rank', 
       y = 'Enrichment score (ES)',
       color = 'Population') +
  theme_classic() +
  scale_x_continuous(breaks = seq(0, 32000, 5000)) +
  # Expand y-axis to give more room for rug marks
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.05))) +  # 10% below, 5% above
  theme(legend.position = "bottom")






