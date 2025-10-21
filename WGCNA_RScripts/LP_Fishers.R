
library(ggplot2)
library(ggrepel)
library(dplyr)
library(purrr)
library(readxl)

BW = read.csv("LP_WGCNA_Output/LP_BW_Summary.csv")

LP_BW = read_xlsx("RNA_Seq_Output/LP_ISO_Ortho_Summary.xlsx")

# Get unique list of colors
colors <- unique(BW$ModuleColor)

## All strain signficiant
# Build a summary dataframe with counts for each color
strain_sig_counts <- lapply(colors, function(col) {
  count <- LP_BW %>%
    filter(LP_BW_Modules == col, strain_SIG == "SIG") %>%
    nrow()
  
  data.frame(ModuleColor = col, Strain_SIG_Count = count)
}) %>%
  bind_rows()

# Now join with BW
BW_with_counts <- BW %>%
  left_join(strain_sig_counts, by = "ModuleColor")


BW_with_counts <- BW_with_counts %>%
  mutate(
    Total_genes = 20733,
    PopDE = 7750,
    PopDE_InModule = Strain_SIG_Count,
    NotPopDE_InModule = ModuleSize - PopDE_InModule,
  )


BW_with_counts <- BW_with_counts %>%
  mutate(
    fisher_table = pmap(list(PopDE_InModule, NotPopDE_InModule, PopDE, Total_genes),
                        ~ matrix(c(..1, ..2, ..3, ..4), nrow = 2,
                                 dimnames = list(PopDE = c("Yes", "No"), InModule = c("Yes", "No")))),
    Fisher_pval = map_dbl(fisher_table, ~ fisher.test(.)$p.value),
    Fisher_FDR = p.adjust(Fisher_pval, method = "fdr"),
    Fisher_SIG = if_else(Fisher_FDR < 0.05, "SIG", "NS"),
    Odds_Ratio = map_dbl(fisher_table, ~ fisher.test(.)$estimate)
  )


fisher_summary <- BW_with_counts %>%
  select(ModuleColor, 
         ModuleSize,
         PopDE_InModule,
         Zsummary.pres,
         Fisher_pval,
         Fisher_FDR,
         Fisher_SIG,
         Odds_Ratio) 

# write.csv(fisher_summary, "LP_WGCNA_Output/LP_BW_Fishers.csv")

BW_with_counts %>%
  ggplot(aes(x = log2(Odds_Ratio), y = -log10(Fisher_FDR), label = ModuleColor)) +
  geom_point(aes(size = ModuleSize), color = "steelblue") +
  geom_text_repel(data = ~filter(., Odds_Ratio > 2 | Fisher_FDR < 1e-5), size = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(x = "Log2(Odds Ratio)", y = "-Log10(FDR)", title = "Enrichment of PopDE genes in Modules") +
  theme_minimal()

library(ggrepel)

BW %>%
  ggplot(aes(x = Module_Size, y = Zsummary.pres, size = PopDE, fill = Module_Color)) + 
  geom_point(shape = 21, color = "black", stroke = 0.5, alpha = 0.8) +
  geom_text_repel(aes(label = PopDE), size = 4, nudge_y = 2) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 2, color = "blue", linetype = "dashed") +
  geom_hline(yintercept = 10, color = "red", linetype = "dashed") +
  scale_fill_identity(guide = "legend", name = "Module Color") +
  scale_y_continuous(limits = c(0, 100)) +
  theme_bw() +
  labs(x = "Module Size",
       y = "Zsummary Preservation")

BW %>%
  ggplot(aes(x = Module_Size, y = Zsummary.pres, size = PopDE, fill = Module_Color)) + 
  geom_point(shape = 21, color = "black", stroke = 0.5, alpha = 0.8) +
  geom_text(aes(label = PopDE), size = 4, color = "black") +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 2, color = "blue", linetype = "dashed") +
  geom_hline(yintercept = 10, color = "red", linetype = "dashed") +
  scale_fill_identity(guide = "legend", name = "Module Color") +
  scale_size_continuous(
    name = "Pop DE Genes",
    breaks = c(10, 100, 500, 1000, 2000, 2500),
    range = c(2, 20),
    trans = "sqrt"
  ) +
  scale_y_continuous(limits = c(0, 100)) +
  theme_bw() +
  labs(x = "Module Size", y = "Zsummary Preservation")









