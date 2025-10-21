
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
    Total_genes = 20733,
    PopDE = 7750,
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
         ModuleLetter,
         ModuleSize,
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
  facet_wrap(~ ModuleColor) +
  theme_void() +
  scale_fill_manual(values = c("PopDE_prop" = "steelblue", "NotPopDE_prop" = "grey90"))


# ggsave("Plots/WGCNA_Plots/WGCNA_Fishers.pdf", width = 12, height = 8, units = "in", dpi = 300)


# INTERACTION --------

# Build a summary dataframe with counts for each color
IXN_sig_counts <- lapply(colors, function(col) {
  count <- EP_BW %>%
    filter(EP_BW_Modules == col, IXN_SIG == "SIG") %>%
    nrow()
  
data.frame(ModuleColor = col, IXN_SIG_Count = count)}) %>%
  bind_rows()

# Now join with BW
BW_with_counts <- BW %>%
  left_join(IXN_sig_counts, by = "ModuleColor")

BW_with_counts <- BW_with_counts %>%
  mutate(
    Total_genes = 20733,
    IxnDE = 43,
    IxnDE_InModule = IXN_SIG_Count,
    NotIxnDE_InModule = ModuleSize - IxnDE_InModule,)


BW_with_counts <- BW_with_counts %>%
  mutate(
    fisher_table = pmap(list(IxnDE_InModule, NotIxnDE_InModule, IxnDE, Total_genes),
                        ~ matrix(c(..1, ..2, ..3, ..4), nrow = 2,
                                 dimnames = list(PopDE = c("Yes", "No"), InModule = c("Yes", "No")))),
    Fisher_pval = map_dbl(fisher_table, ~ fisher.test(.)$p.value),
    Fisher_FDR = p.adjust(Fisher_pval, method = "fdr"),
    Fisher_SIG = if_else(Fisher_FDR < 0.05, "SIG", "NS"),
    Odds_Ratio = map_dbl(fisher_table, ~ fisher.test(.)$estimate))


fisher_summary <- BW_with_counts %>%
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

# write.csv(fisher_summary, "EP_WGCNA_Output/EP_BW_Fishers.csv")






