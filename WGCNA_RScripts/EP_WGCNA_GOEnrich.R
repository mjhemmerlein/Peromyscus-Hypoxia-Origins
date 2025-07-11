
library(WGCNA)
library(dplyr)
library(edgeR)
library(lme4)
library(lmerTest)
library(readxl)
library(ggplot2)

gene_infoBWEP = read.csv("EP_WGCNA_Output/EP_BW_Modules.csv")


# GO Analysis for BW modules ####

musGeneinfo = read_xlsx("RNA_Seq_Output/EP_ISO_Ortho_Summary.xlsx")

## Match and add gene names from Mus from Kates musGeneInfo
## Filter genes from module colors that are in musGeneInfo

# Module colors
colors <- c("green", "pink", "red", "magenta", "blue", "brown", "purple", "black", "turquoise", "tan", "greenyellow", "yellow")
  
# Loop through each color
  for (color in colors) {
    print(paste("Processing color:", color))
    

colorME <- gene_infoBWEP %>% filter(moduleColor == color)
colorME <- as.data.frame(colorME)
colorME <- colorME[, 1:2]
rownames(colorME) <- colorME$Gene
intersection2 <- musGeneinfo$Pman_GeneID %in% rownames(colorME)
colorME0 <- musGeneinfo[which(intersection2), ] 
colorMEGO <- as.vector(colorME0$Pman_GeneID)
    
colorGO <- gost(colorMEGO,
                organism = "mmusculus",
                user_threshold = 0.05,
                custom_bg = background,
                ordered_query = TRUE,
                correction_method = "fdr")
    
# Create results dataframe
colorGOResults <- data.frame(
  Cluster = colorGO$result$query,
  Term.ID = colorGO$result$term_id,
  Term.Name = colorGO$result$term_name,
  geneid = colorGO$result$intersection,
  P.value = colorGO$result$p_value,
  Source = colorGO$result$source,
  Term.Size = colorGO$result$term_size,
  Precision = colorGO$result$precision,
  intersection_size = colorGO$result$intersection_size,
  query_size = colorGO$result$query_size,
  Gene_ratio = as.numeric((colorGO$result$intersection_size/colorGO$result$query_size)))
    
# Save to CSV
    filename <- paste0("EP_WGCNA_Output/GO_Results/", color, "_GO_results.csv")
    write.csv(colorGOResults, filename, row.names = FALSE)
    print(paste("Saved", nrow(colorGOResults), "GO terms to", filename))
  }