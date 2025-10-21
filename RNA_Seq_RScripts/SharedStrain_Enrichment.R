
# Gene Ontology for Early Pregnancy

# Libraries
library(tidyverse)
library(readr)
library(readxl)
library(gprofiler2)

# Shared across EP and LP
EP_LP_strain = read_csv("RNA_Seq_Output/EP_LP_SharedStrain.csv")

background = EP_LP_strain$Mus_GeneID

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


# Up-Up
EP_POP_UP = EP_LP_strain %>%
  filter(Quadrant == "Up-Up") %>%
  pull(Mus_GeneID)

EP_POP_UP_GOResults = run_gost_analysis(EP_POP_UP, background)

# Down-Down

EP_POP_DOWN = EP_LP_strain %>%
  filter(Quadrant == "Down-Down") %>%
  pull(Mus_GeneID)

EP_POP_DOWN_GOResults = run_gost_analysis(EP_POP_DOWN, background)

# Up-Down
EP_POP_UPDOWN = EP_LP_strain %>%
  filter(Quadrant == "Up-Down") %>%
  pull(Mus_GeneID)

EP_POP_UPDOWN_GOResults = run_gost_analysis(EP_POP_UPDOWN, background)

# Down-Up
EP_POP_DOWNUP = EP_LP_strain %>%
  filter(Quadrant == "Down-Up") %>%
  pull(Mus_GeneID)

EP_POP_DOWNUP_GOResults = run_gost_analysis(EP_POP_DOWNUP, background)


