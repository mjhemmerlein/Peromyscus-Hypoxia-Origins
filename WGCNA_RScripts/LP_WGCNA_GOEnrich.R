
library(WGCNA)
library(dplyr)
library(edgeR)
library(lme4)
library(lmerTest)
library(readxl)
library(ggplot2)
library(gprofiler2)

LP_ISO_Summary = read_xlsx("RNA_Seq_Output/LP_ISO_Ortho_Summary.xlsx")

# Define the function
run_gost_analysis <- function(gene_list, background, organism = "mmusculus") {
  gost_results <- gost(
    gene_list, 
    organism = organism,
    custom_bg = background,
    ordered_query = FALSE,
    correction_method = "fdr",
    evcodes = TRUE
  )
  
  if (is.null(gost_results)) {
    return(NULL)
  }
  
  data.frame(
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
    Gene_ratio = gost_results$result$intersection_size / gost_results$result$query_size
  )
}

# Background gene list
background <- LP_ISO_Summary$Mus_GeneID

# Unique module colors
colors <- unique(LP_ISO_Summary$LP_BW_Modules)

# Loop over module colors
for (col in colors) {
  gene_list <- LP_ISO_Summary %>%
    filter(LP_BW_Modules == col) %>%
    pull(Mus_GeneID)
  
  go_result <- run_gost_analysis(gene_list, background)
  
  if (!is.null(go_result)) {
    write.csv(go_result, paste0("LP_WGCNA_Output/GO_Results/", col, ".csv"), row.names = FALSE)
  }
}







