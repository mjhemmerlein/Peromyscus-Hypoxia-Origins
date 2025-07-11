# Gene Ontology for Early Pregnancy
# Derived from EP_ISO_Summary file

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


background = EP_ISO_Summary$Pman_GeneID

# Early Pregnancy
# Population
EP_POP = EP_ISO_Summary %>%
  filter(strain_SIG == "SIG") %>%
  pull(Pman_GeneID)

EP_POP_GOResults = run_gost_analysis(EP_POP, background)

# Pop Up
EP_POP_UP = EP_ISO_Summary %>%
  filter(strain_SIG == "SIG") %>%
  filter(O2_logFC > 1) %>%
  pull(Pman_GeneID)

EP_POP_UP_GOResults = run_gost_analysis(EP_POP_UP, background)


EP_POP_DOWN = EP_ISO_Summary %>%
  filter(strain_SIG == "SIG") %>%
  filter(O2_logFC < -1) %>%
  pull(Pman_GeneID)

EP_POP_DOWN_GOResults = run_gost_analysis(EP_POP_DOWN, background)



# Hypoxia
EP_O2 = EP_ISO_Summary %>%
  filter(O2_SIG == "SIG") %>%
  pull(Pman_GeneID)

EP_O2_GOResults = run_gost_analysis(EP_O2, background)

# Hypoxia Up
EP_O2_UP = EP_ISO_Summary %>%
  filter(O2_SIG == "SIG") %>%
  filter(O2_logFC > 0.5) %>%
  pull(Pman_GeneID)

EP_O2_UP_GOResults = run_gost_analysis(EP_O2_UP, background)

# Hypoxia Down
EP_O2_DOWN = EP_ISO_Summary %>%
  filter(O2_SIG == "SIG") %>%
  filter(O2_logFC < -0.5) %>%
  pull(Pman_GeneID)

EP_O2_DOWN_GOResults = run_gost_analysis(EP_O2_DOWN, background)





# Interaction
EP_IXN = EP_ISO_Summary %>%
  filter(IXN_SIG == "SIG") %>%
  pull(Pman_GeneID)

EP_IXN_GOResults = run_gost_analysis(EP_IXN, background)

# Interaction Up

count(EP_ISO_Summary %>%
  filter(IXN_SIG == "SIG") %>%
  filter(IXN_logFC > 0.5))

EP_IXN_UP = EP_ISO_Summary %>%
  filter(IXN_SIG == "SIG") %>%
  filter(IXN_logFC > 0.5) %>%
  pull(Pman_GeneID)

EP_IXN_UP_GOResults = run_gost_analysis(EP_IXN_UP, background)

# Interaction Down
EP_IXN_DOWN = EP_ISO_Summary %>%
  filter(IXN_SIG == "SIG") %>%
  filter(IXN_logFC < -0.5) %>%
  pull(Pman_GeneID)

EP_IXN_DOWN_GOResults = run_gost_analysis(EP_IXN_DOWN, background)


# BW
EP_BW = EP_ISO_Summary %>%
  filter(BW_O2_SIG == "SIG") %>%
  pull(Pman_GeneID)

EP_BW_GOResults = run_gost_analysis(EP_BW, background)


# ME Up
EP_ME = EP_ISO_Summary %>%
  filter(ME_O2_SIG == "SIG") %>%
  pull(Pman_GeneID)

EP_ME_GOResults = run_gost_analysis(EP_ME, background)


