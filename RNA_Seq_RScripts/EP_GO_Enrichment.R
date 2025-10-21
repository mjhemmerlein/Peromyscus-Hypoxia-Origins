
# Libraries
library(tidyverse)
library(readr)
library(readxl)
library(gprofiler2)


# EP ISO interesting genes
# Gene Counts

EP_ISO_Summary = read_xlsx("RNA_Seq_Output/EP_ISO_Ortho_Summary.xlsx")

run_gost_analysis <- function(gene_list, background, organism = "mmusculus") {
  
gost_results = gost(
    gene_list, 
    organism = organism,
    custom_bg = background,
    ordered_query = F,
    correction_method = "fdr",
    evcodes = T)

results_df = data.frame(
    Cluster = gost_results$result$query,
    Term.ID = gost_results$result$term_id,
    Term.Name = gost_results$result$term_name,
    geneid = gost_results$result$intersection,
    P.value = gost_results$result$p_value,
    Source = gost_results$result$source,
    Term.Size = gost_results$result$term_size,
    Precision = gost_results$result$precision,
    intersection_size = gost_results$result$intersection_size,
    query_size = gost_results$result$query_size,
    Gene_ratio = as.numeric((gost_results$result$intersection_size/gost_results$result$query_size)))

}


background = EP_ISO_Summary$Mus_GeneID


# Population -----
EP_POP = EP_ISO_Summary %>%
  filter(strain_SIG == "SIG") %>%
  pull(Mus_GeneID)

EP_POP_GOResults = run_gost_analysis(EP_POP, background)

# Pop Up
EP_POP_UP = EP_ISO_Summary %>%
  filter(strain_SIG == "SIG") %>%
  filter(strain_logFC > 0) %>%
  pull(Mus_GeneID)

EP_POP_UP_GOResults = run_gost_analysis(EP_POP_UP, background)


EP_POP_DOWN = EP_ISO_Summary %>%
  filter(strain_SIG == "SIG") %>%
  filter(strain_logFC < 0) %>%
  pull(Mus_GeneID)

EP_POP_DOWN_GOResults = run_gost_analysis(EP_POP_DOWN, background)


# Absolute median of pop
median_effect_POP <- EP_ISO_Summary %>%
  filter(strain_SIG == "SIG")

median(abs(median_effect_POP$strain_logFC))

PBScount = EP_ISO_Summary %>%
  filter(strain_SIG == "SIG") %>%
  filter(PBS == "TRUE")

EP_POP_PBS_UP = PBScount %>%
  filter(strain_logFC > 0) %>%
  pull(Mus_GeneID)

EP_POP_PBS_UP_GOResults = run_gost_analysis(EP_POP_PBS_UP, background)

EP_POP_PBS_DOWN = PBScount %>%
  filter(strain_logFC < 0) %>%
  pull(Mus_GeneID)

EP_POP_PBS_DOWN_GOResults = run_gost_analysis(EP_POP_PBS_DOWN, background)


# Hypoxia --------

EP_O2 = EP_ISO_Summary %>%
  filter(O2_SIG == "SIG") %>%
  pull(Mus_GeneID)

EP_O2_GOResults = run_gost_analysis(EP_O2, background)

# Hypoxia Up
EP_O2_UP = EP_ISO_Summary %>%
  filter(O2_SIG == "SIG") %>%
  filter(O2_logFC > 0) %>%
  pull(Mus_GeneID)

EP_O2_UP_GOResults = run_gost_analysis(EP_O2_UP, background)

# write.csv(EP_O2_UP_GOResults, "EP_GOEnrichment_Ouput/EP_Hypoxia_UP.csv")

# Hypoxia Down
EP_O2_DOWN = EP_ISO_Summary %>%
  filter(O2_SIG == "SIG") %>%
  filter(O2_logFC < 0) %>%
  pull(Mus_GeneID)

EP_O2_DOWN_GOResults = run_gost_analysis(EP_O2_DOWN, background)

# write.csv(EP_O2_DOWN_GOResults, "EP_GOEnrichment_Output/EP_Hypoxia_Down.csv")

# Absolute median of hypoxia
median_effect_O2 <- EP_ISO_Summary %>%
  filter(O2_SIG == "SIG")

median(abs(median_effect_O2$O2_logFC))

# Interaction -------

EP_IXN = EP_ISO_Summary %>%
  filter(IXN_SIG == "SIG")

median_effect_size <- median(abs(EP_IXN$IXN_logFC))

EP_IXN = EP_ISO_Summary %>%
  filter(IXN_SIG == "SIG") %>%
  pull(Mus_GeneID)

EP_IXN_GOResults = run_gost_analysis(EP_IXN, background)

# write.csv(EP_IXN_GOResults, "EP_GOEnrichment_Output/EP_IXN.csv")

# Interaction Up

count(EP_ISO_Summary %>%
  filter(IXN_SIG == "SIG") %>%
  filter(IXN_logFC > 0))

EP_IXN_UP = EP_ISO_Summary %>%
  filter(IXN_SIG == "SIG") %>%
  filter(IXN_logFC > 0) %>%
  pull(Mus_GeneID)

EP_IXN_UP_GOResults = run_gost_analysis(EP_IXN_UP, background)

# Interaction Down

count(EP_ISO_Summary %>%
        filter(IXN_SIG == "SIG") %>%
        filter(IXN_logFC < 0))

EP_IXN_DOWN = EP_ISO_Summary %>%
  filter(IXN_SIG == "SIG") %>%
  filter(IXN_logFC < 0) %>%
  pull(Mus_GeneID)

EP_IXN_DOWN_GOResults = run_gost_analysis(EP_IXN_DOWN, background)


# BW
EP_BW = EP_ISO_Summary %>%
  filter(BW_O2_SIG == "SIG") %>%
  pull(Mus_GeneID)

EP_BW_GOResults = run_gost_analysis(EP_BW, background)

median(abs(EP_ISO_Summary$BW_O2_logFC))

# ME
EP_ME = EP_ISO_Summary %>%
  filter(ME_O2_SIG == "SIG") %>%
  pull(Mus_GeneID)

EP_ME_GOResults = run_gost_analysis(EP_ME, background)

median(abs(EP_ISO_Summary$ME_O2_logFC))
