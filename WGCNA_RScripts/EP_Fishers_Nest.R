
library(ggplot2)
library(ggrepel)
library(dplyr)
library(purrr)
library(readxl)
library(tidyverse)

BW = read.csv("EP_WGCNA_Output/EP_BW_Summary.csv")

EP_BW = read_xlsx("RNA_Seq_Output/EP_ISO_Ortho_Summary.xlsx")

# Define your subsets
IXN <- EP_BW %>%
  filter(IXN_SIG == "SIG")

POP <- EP_BW %>%
  filter(strain_SIG == "SIG", IXN_SIG != "SIG")

HYP <- EP_BW %>%
  filter(O2_SIG == "SIG", IXN_SIG != "SIG")

# Get unique module colors
colors <- unique(BW$ModuleColor)


# NESTED FISHERS --------
# Define compute_fisher function
compute_fisher <- function(df, module_col, total_genes, Sig_genes, group_name = "Group") {
  
  sig_counts <- df %>%
    count(.data[[module_col]], name = "SIG_Count") %>%
    rename(ModuleColor = .data[[module_col]])
  
  # Join counts with module info from BW
  df_with_counts <- BW %>%
    left_join(sig_counts, by = "ModuleColor") %>%
    mutate(
      SIG_Count = replace_na(SIG_Count, 0),
      Total_genes = total_genes,
      Sig_genes = Sig_genes,
      InModule = SIG_Count
    )
  
  # Perform Fisher’s exact test per module
  df_with_counts <- df_with_counts %>%
    mutate(
      fisher_table = pmap(list(InModule, ModuleSize, Sig_genes, Total_genes),
                          ~ matrix(c(..1, ..2, ..3, ..4), nrow = 2,
                                   dimnames = list( Table1= c("DiffExp", "Total"),
                                                    Table2 = c("Module", "Transcriptome")))),
      fisher_result = map(fisher_table, ~ fisher.test(.)),
      Fisher_pval = map_dbl(fisher_result, ~ .$p.value),
      Odds_Ratio = map_dbl(fisher_result, ~ .$estimate),
      CI_lower = map_dbl(fisher_result, ~ .$conf.int[1]),
      CI_upper = map_dbl(fisher_result, ~ .$conf.int[2]),
      Fisher_FDR = p.adjust(Fisher_pval, method = "fdr"),
      Fisher_SIG = if_else(Fisher_FDR < 0.05, "SIG", "NS"),
      Factor = group_name
    )
  
  # Return summary table
  df_with_counts %>%
    select(ModuleColor, 
           ModuleLetter,
           InModule,
           ModuleSize,
           Sig_genes, 
           Total_genes,
           fisher_table,
           Zsummary.pres, 
           Fisher_pval, 
           Fisher_FDR, 
           Fisher_SIG, 
           Odds_Ratio, 
           CI_lower, 
           CI_upper, 
           Factor)
}

# RESULTS -----------
# POPULATION
fisher_POP <- compute_fisher(POP, "EP_BW_Modules",
                             total_genes = 20938,
                             Sig_genes = nrow(POP),
                             group_name = "POP")

print(fisher_POP$fisher_table[[1]])

# fisher_POP_filtered <- fisher_POP %>%
#  select(-fisher_table)
# write.csv(fisher_POP_filtered, "EP_WGCNA_Output/EP_Fishers_POP.csv", row.names = FALSE)

# HYPOXIA
fisher_HYP <- compute_fisher(HYP, "EP_BW_Modules",
                             total_genes = 20938,
                             Sig_genes = nrow(HYP),
                             group_name = "HYP")

print(fisher_HYP$fisher_table[[1]])

# fisher_HYP_filtered <- fisher_HYP %>%
#  select(-fisher_table)
# write.csv(fisher_HYP_filtered, "EP_WGCNA_Output/EP_Fishers_HYP.csv", row.names = FALSE)

fisher_IXN <- compute_fisher(IXN, "EP_BW_Modules",
                             total_genes = 20938,
                             Sig_genes = nrow(IXN),
                             group_name = "IXN")

print(fisher_IXN$fisher_table[[1]])

# fisher_IXN_filtered <- fisher_IXN %>%
#  select(-fisher_table)
# write.csv(fisher_IXN_filtered, "EP_WGCNA_Output/EP_Fishers_IXN.csv", row.names = FALSE)




# PLOTTING ------

module_size_order <- fisher_HYP %>%
  arrange(desc(ModuleSize)) %>%
  pull(ModuleLetter) %>%
  unique()  # ensure unique values

fisher_all <- bind_rows(fisher_POP, fisher_HYP, fisher_IXN) %>%
  # Remove rows with invalid Odds Ratios
  filter(Odds_Ratio > 0, !is.infinite(Odds_Ratio), !is.nan(Odds_Ratio)) %>%
  mutate(
    # Calculate log2 values (now all valid)
    log2OR = log2(Odds_Ratio),
    log2CI_low = log2(CI_lower),
    log2CI_high = log2(CI_upper),
    
    # Set factor levels based on module size (largest to smallest)
    ModuleLetter = factor(ModuleLetter, levels = module_size_order),
    Factor = factor(Factor, levels = c("Strain", "HYP", "IXN"))
  )


ggplot(fisher_all, aes(x = log2OR, y = ModuleLetter)) +
  geom_vline(xintercept = 0, color = "gray40", linetype = "solid", linewidth = 0.8) +
  
  geom_errorbarh(
    aes(xmin = log2CI_low, xmax = log2CI_high, group = Factor, color = Fisher_SIG),
    height = 0,
    linewidth = 1,
    position = position_dodge(width = 0.7)
  ) +

  geom_point(
    aes(shape = Factor, color = Fisher_SIG, group = Factor),
    size = 4.5,
    position = position_dodge(width = 0.7)
  ) +
  
  scale_color_manual(
    values = c("NS" = "gray40", "SIG" = "#DF4B68"),
    name = "Significance"
  ) +
  
  scale_shape_manual(
    values = c("Strain" = 16, "HYP" = 17, "IXN" = 15),
    name = "Group"
  ) +

  labs(
    x = expression(bold(log[2]*"(Odds Ratio)")),
    y = "Module",
    title = "Module Enrichment by Group (Fisher Test)"
  ) +
  
  scale_x_continuous(
    limits = c(-8, 8),
    breaks = seq(-8, 8, 1),     
    labels = seq(-8, 8, 1),     
    minor_breaks = seq(-8, 8, 0.5)
  ) +
  
  theme_bw(base_size = 12) +
  theme(
    # Grid
    panel.grid.major.y = element_line(color = "gray95", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
    
    # Axes
    axis.title = element_text(face = "bold", size = 12),
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 10),
    
    # Title and legend
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    legend.position = "top",
    legend.box = "horizontal",
    legend.text = element_text(size = 10)
  )

# ggsave("Plots/WGCNA_Plots/WGCNA_Fishers_Forest.pdf", width = 6, height = 6, units = "in", dpi = 300)



# Raw Odds Ratios

fisher_all <- bind_rows(fisher_POP, fisher_HYP, fisher_IXN)

manual_levels <- c(
  "A_POP","A_HYP","A_IXN",
  "B_POP","B_HYP","B_IXN",
  "C_POP","C_HYP","C_IXN",
  "D_POP","D_HYP","D_IXN",
  "E_POP","E_HYP","E_IXN",
  "F_POP","F_HYP","F_IXN",
  "G_POP","G_HYP","G_IXN",
  "H_POP","H_HYP","H_IXN",
  "I_POP","I_HYP","I_IXN",
  "J_POP","J_HYP","J_IXN",
  "K_POP","K_HYP","K_IXN",
  "L_POP","L_HYP","L_IXN"
)

# Assign Module_Type and set factor levels
fisher_all <- fisher_all %>%
  mutate(
    Module_Type = paste(ModuleLetter, Factor, sep="_"),
    Module_Type = factor(Module_Type, levels = rev(manual_levels)),
    
    # Placeholder for NS points
    OR_plot      = ifelse(Fisher_SIG == "NS", 0.01, Odds_Ratio),
    CI_low_plot  = ifelse(Fisher_SIG == "NS", 0.01, CI_lower),
    CI_high_plot = ifelse(Fisher_SIG == "NS", 0.01, CI_upper),
    color_plot   = ifelse(Fisher_SIG == "NS", "gray", "black")
  )


ggplot(fisher_all, aes(x = OR_plot, y = Module_Type, shape = Factor)) +
  geom_errorbarh(aes(xmin = CI_low_plot, xmax = CI_high_plot, color = color_plot),
                 height = 0.2, linewidth = 0.5) +
  geom_point(aes(color = color_plot), size = 4) +
  geom_vline(xintercept = 1, color = "gray40") +
  scale_color_identity() +
  scale_x_continuous(
    limits = c(0.01, 100),
    trans = "log10",
    breaks = c(0.01,0.1,1,10,50,100),
    labels = c("0.01","0.1","1","10","50","100")
  ) +
  labs(x = "Odds Ratio (log10 scale)", y = "Module") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(face = "bold"),
    panel.grid.major.y = element_line(color = "gray95", linewidth = 0.3),
    panel.grid.major.x = element_line(color = "gray95", linewidth = 0.3),
  ) +
  scale_shape_manual(
    values = c("POP" = 1, "HYP" = 0, "IXN" = 2),
    name = "Factor")


# ggsave("Plots/WGCNA_Plots/WGCNA_Fishers_Forestv2.pdf", width = 7, height = 6, units = "in", dpi = 300)





