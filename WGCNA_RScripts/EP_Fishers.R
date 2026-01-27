
library(ggplot2)
library(ggrepel)
library(dplyr)
library(purrr)
library(readxl)
library(tidyverse)

BW = read.csv("EP_WGCNA_Output/EP_BW_Summary.csv")

EP_BW = read_xlsx("RNA_Seq_Output/EP_ISO_Ortho_Summary.xlsx")

# Get unique list of colors
colors <- unique(BW$ModuleColor)

# ANCESTRY ------

strain_sig_counts <- lapply(colors, function(col) {
  count <- EP_BW %>%
    filter(EP_BW_Modules == col, strain_SIG == "SIG") %>%
    nrow()
  
  data.frame(ModuleColor = col, Strain_SIG_Count = count)}) %>%
  bind_rows()

BW_with_counts <- BW %>%
  left_join(strain_sig_counts, by = "ModuleColor")

BW_with_counts <- BW_with_counts %>%
  mutate(
    Total_genes = 20938,
    PopDE = 7596,
    PopDE_InModule = Strain_SIG_Count,
    NotPopDE_InModule = ModuleSize - PopDE_InModule,)

BW_with_counts <- BW_with_counts %>%
  mutate(
    fisher_table = pmap(list(PopDE_InModule, NotPopDE_InModule, PopDE, Total_genes),
                        ~ matrix(c(..1, ..2, ..3, ..4), nrow = 2,
                                 dimnames = list(PopDE = c("Yes", "No"), InModule = c("Yes", "No")))),
    Fisher_pval = map_dbl(fisher_table, ~ fisher.test(.)$p.value),
    Fisher_FDR = p.adjust(Fisher_pval, method = "fdr"),
    Fisher_SIG = if_else(Fisher_FDR < 0.05, "SIG", "NS"),
    Odds_Ratio = map_dbl(fisher_table, ~ fisher.test(.)$estimate))


fisher_summary <- BW_with_counts %>%
  select(ModuleColor,
         ModuleSize,
         ModuleLetter,
         PopDE_InModule,
         NotPopDE_InModule,
         PopDE,
         Total_genes,
         Zsummary.pres,
         Fisher_pval,
         Fisher_FDR,
         Fisher_SIG,
         Odds_Ratio) 

fisher_summary_sig = fisher_summary %>%
  filter(Fisher_SIG == "SIG")

# write.csv(fisher_summary, "EP_WGCNA_Output/EP_BW_Fishers.csv")

# PIE CHARTS

BW_long <- BW_with_counts %>%
  mutate(
    PopDE_prop = PopDE_InModule / ModuleSize,
    NotPopDE_prop = 1 - PopDE_prop) %>%
  select(ModuleColor, ModuleLetter, PopDE_prop, NotPopDE_prop) %>%
  pivot_longer(cols = c(PopDE_prop, NotPopDE_prop),
               names_to = "Category", values_to = "Prop")


ggplot(BW_long, aes(x = "", y = Prop, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y") +
  facet_wrap(~ ModuleLetter) +
  theme_void() +
  scale_fill_manual(values = c("PopDE_prop" = "steelblue", "NotPopDE_prop" = "grey90"))


# ggsave("Plots/WGCNA_Plots/WGCNA_Fishers.pdf", width = 12, height = 8, units = "in", dpi = 300)


BW_with_counts %>%
  ggplot(aes(x = ModuleLetter, y = Odds_Ratio)) +
  geom_point(aes(size = -log10(Fisher_FDR), color = Odds_Ratio > 1)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_color_manual(values = c("TRUE" = "steelblue", "FALSE" = "grey60")) +
  theme_minimal(base_size = 12) +
  labs(
    y = "Odds Ratio (PopDE Enrichment)",
    x = "Module",
    color = "Direction",
    size = "-log10(Fisher p-value)"
  )







# INTERACTION --------

# Build a summary dataframe with counts for each color
IXN_sig_counts <- lapply(colors, function(col) {
  count <- EP_BW %>%
    filter(EP_BW_Modules == col, IXN_SIG == "SIG") %>%
    nrow()
  
  data.frame(ModuleColor = col, IXN_SIG_Count = count)}) %>%
  bind_rows()

# Now join with BW
BW_with_counts_Ixn <- BW %>%
  left_join(IXN_sig_counts, by = "ModuleColor")

BW_with_counts_Ixn <- BW_with_counts_Ixn %>%
  mutate(
    Total_genes = 20938,
    IxnDE = 28,
    IxnDE_InModule = IXN_SIG_Count,
    NotIxnDE_InModule = ModuleSize - IxnDE_InModule,)


BW_with_counts_Ixn <- BW_with_counts_Ixn %>%
  mutate(
    fisher_table = pmap(list(IxnDE_InModule, NotIxnDE_InModule, IxnDE, Total_genes),
                        ~ matrix(c(..1, ..2, ..3, ..4), nrow = 2,
                                 dimnames = list(PopDE = c("Yes", "No"), InModule = c("Yes", "No")))),
    Fisher_pval = map_dbl(fisher_table, ~ fisher.test(.)$p.value),
    Fisher_FDR = p.adjust(Fisher_pval, method = "fdr"),
    Fisher_SIG = if_else(Fisher_FDR < 0.05, "SIG", "NS"),
    Odds_Ratio = map_dbl(fisher_table, ~ fisher.test(.)$estimate))


fisher_summary <- BW_with_counts_Ixn %>%
  select(ModuleColor,
         ModuleLetter,
         ModuleSize,
         IxnDE_InModule,
         NotIxnDE_InModule,
         IxnDE,
         Total_genes,
         Zsummary.pres,
         Fisher_pval,
         Fisher_FDR,
         Fisher_SIG,
         Odds_Ratio) 

fisher_summary_sig = fisher_summary %>%
  filter(Fisher_SIG == "SIG")

# write.csv(fisher_summary, "EP_WGCNA_Output/EP_BW_IXN_Fishers.csv")

BW_long <- BW_with_counts_Ixn %>%
  mutate(
    IxnDE_prop = IxnDE_InModule / ModuleSize,
    NotIxnDE_prop = 1 - IxnDE_prop) %>%
  select(ModuleColor, ModuleLetter, IxnDE_prop, NotIxnDE_prop) %>%
  pivot_longer(cols = c(IxnDE_prop, NotIxnDE_prop),
               names_to = "Category", values_to = "Prop")


ggplot(BW_long, aes(x = "", y = Prop, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y") +
  facet_wrap(~ ModuleLetter) +
  theme_void() +
  scale_fill_manual(values = c("IxnDE_prop" = "steelblue", "IxnPopDE_prop" = "grey90"))



# Combine
fisher_combined <- bind_rows(
  BW_with_counts %>%
    select(ModuleLetter, Odds_Ratio, Fisher_FDR) %>%
    mutate(Category = "PopDE"),
  BW_with_counts_Ixn %>%
    select(ModuleLetter, Odds_Ratio, Fisher_FDR) %>%
    mutate(Category = "IxnDE")
)

library(ggplot2)
ggplot(fisher_combined, aes(x = ModuleLetter, y = Category, fill = log2(Odds_Ratio))) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato", midpoint = 0) +
  geom_text(aes(label = ifelse(Fisher_FDR < 0.05, "*", "")), color = "black") +
  theme_minimal(base_size = 12) +
  labs(
    x = "Module",
    y = "",
    fill = "log2(OR)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# UNESTED FISHERS ------
# POP PURPLE

Sig_genes <- nrow(POP)
total_genes <- 20938
group_name <- "Strain"

# Count significant genes in this module
InModule <- POP %>% 
  filter(EP_BW_Modules == "purple") %>% 
  nrow()

# Get module info from BW
module_info <- BW %>% 
  filter(ModuleColor == "purple") %>% 
  select(-X)

# Build the Fisher table (2x2)
fisher_table <- matrix(
  c(
    InModule,
    module_info$ModuleSize,
    Sig_genes,
    total_genes),
  nrow = 2,
  dimnames = list(
    Table1 = c("DiffExp", "Total"),
    Table2 = c("Module", "Transcriptome")
  )
)

# Run Fisher's test
fisher_result <- fisher.test(fisher_table)

# Collect results
Fisher_summary <- tibble(
  Factor = group_name,
  Fisher_pval = fisher_result$p.value,
  Odds_Ratio = fisher_result$estimate,
  CI_lower = fisher_result$conf.int[1],
  CI_upper = fisher_result$conf.int[2]
)




