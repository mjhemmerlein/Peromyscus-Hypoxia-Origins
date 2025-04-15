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

EP_POP_PBS = EP_ISO_Summary %>%
  filter(strain_SIG == "SIG") %>%
  filter(PBS == "TRUE") %>%
  pull(Pman_GeneID)

EP_POP_PBS_GOResults = run_gost_analysis(EP_POP_PBS, background)

# Hypoxia
EP_O2 = EP_ISO_Summary %>%
  filter(O2_SIG == "SIG") %>%
  pull(Pman_GeneID)

EP_O2_GOResults = run_gost_analysis(EP_O2, background)

# Interaction
EP_IXN = EP_ISO_Summary %>%
  filter(IXN_SIG == "SIG") %>%
  pull(Pman_GeneID)

EP_IXN_GOResults = run_gost_analysis(EP_IXN, background)




# Shared across EP and LP
EP_LP_strain = read_csv("RNA_Seq_Output/EP_LP_SharedStrain.csv")

# EP Population Upregulated
EP_POP_UP = EP_LP_strain %>%
  filter(diffexpressed_EP == "UP") %>%
  pull(Pman_GeneID)

EP_POP_UP_GOResults = run_gost_analysis(EP_POP_UP, background)

# EP Population Downregulated
EP_POP_DOWN = EP_LP_strain %>%
  filter(diffexpressed_EP == "DOWN") %>%
  filter(diffexpressed_LP == "DOWN") %>%
  pull(Pman_GeneID)

EP_POP_DOWN_GOResults = run_gost_analysis(EP_POP_DOWN, background)

# Up-Up
UP_UP = EP_LP_strain %>%
  filter(Quadrant == "Up-Up") %>%
  pull(Pman_GeneID)

UP_UP_GOResults = run_gost_analysis(UP_UP, background)

# Down-Down
DOWN_DOWN = EP_LP_strain %>%
  filter(Quadrant == "Down-Down") %>%
  pull(Pman_GeneID)

DOWN_DOWN_GOResults = run_gost_analysis(DOWN_DOWN, background)

