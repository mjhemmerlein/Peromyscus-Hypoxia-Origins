# Libraries
library(tidyverse)
library(readr)
library(readxl)
library(pheatmap)
library(BiocParallel)
library(dplyr)
library(ggplot2)
library(readxl)
library(variancePartition)
library(edgeR)
library(Matrix)

# Process Early Pregnancy (EP) Data
EP_Sample_Info <- read_xlsx("RNA_Seq_RawData/MetaData_EP.xlsx")
EP_Sample_Info <- as.data.frame(EP_Sample_Info)
rownames(EP_Sample_Info) <- EP_Sample_Info$Sample_ID
EP_Sample_Info <- EP_Sample_Info %>% filter(Strain == "BW" | Strain == "ME")
EP_Sample_Info$Group <- as.factor(EP_Sample_Info$Group)
EP_Sample_Info$Group_Timepoint = paste0("EP_", EP_Sample_Info$Group)

# Read EP Raw Data
Pman_rawreads_EP <- read_xlsx("RNA_Seq_RawData/EP_Pman_ExtMMFrac_readcounts.xlsx")
Pman_rawreads_EP <- as.data.frame(Pman_rawreads_EP)
Pman_rawreads_EP <- Pman_rawreads_EP %>% filter(!is.na(Geneid))
row.names(Pman_rawreads_EP) <- Pman_rawreads_EP$Geneid
Pman_rawreads_EP <- Pman_rawreads_EP[,-c(1:6)]

# Check and match columns
Check_EP <- EP_Sample_Info$Seq_Name
colnames(Pman_rawreads_EP) == Check_EP
colnames(Pman_rawreads_EP) <- rownames(EP_Sample_Info)

# Prepare EP data for edgeR
Pman_readcounts_EP <- as.matrix(Pman_rawreads_EP)
dPman_0_EP <- DGEList(Pman_readcounts_EP)
dPman_0_EP <- calcNormFactors(dPman_0_EP)
keep_EP <- rowSums(cpm(dPman_0_EP) > 0.5) >= 60
dPman_EP <- dPman_0_EP[keep_EP,]
dim(dPman_EP)

# Compute logCPM matrix for EP data
logCPM_matrix_EP <- cpm(dPman_EP, log = TRUE)

# Process Late Pregnancy (LP) Data
LP_Sample_Info <- read_xlsx("RNA_Seq_RawData/MetaData_LP.xlsx")
LP_Sample_Info <- as.data.frame(LP_Sample_Info)
LP_Sample_Info <- LP_Sample_Info[-c(80,81),]
rownames(LP_Sample_Info) <- LP_Sample_Info$Sample_ID_LZ
LP_Sample_Info <- LP_Sample_Info %>% filter(Strain == "BW" | Strain == "ME")
LP_Sample_Info <- LP_Sample_Info[-which(rownames(LP_Sample_Info) == "LZ089"),]
LP_Sample_Info$Group_Timepoint = paste0("LP_", LP_Sample_Info$Group)
LP_Sample_Info <- LP_Sample_Info %>% rename(Sample_ID = Sample_ID_LZ)

# Read LP Raw Data
Pman_rawreads_LP <- read_xlsx("RNA_Seq_RawData/LP_Pman_ExtMMFrac_readcounts.xlsx")
Pman_rawreads_LP <- as.data.frame(Pman_rawreads_LP)
Pman_rawreads_LP <- Pman_rawreads_LP %>% filter(!is.na(Geneid))
row.names(Pman_rawreads_LP) <- Pman_rawreads_LP$Geneid
Pman_rawreads_LP <- Pman_rawreads_LP[,-c(1:6)]
Pman_rawreads_LP <- subset(Pman_rawreads_LP, select = -c(RNA201216ZC_LZ089_S16_L001_fastp_pman_Halign_liberal.bam))

# Check and match columns
Check_LP <- LP_Sample_Info$Seq_Name
colnames(Pman_rawreads_LP) == Check_LP
colnames(Pman_rawreads_LP) <- rownames(LP_Sample_Info)

# Prepare LP data for edgeR
Pman_readcounts_LP <- as.matrix(Pman_rawreads_LP)
dPman_0_LP <- DGEList(Pman_readcounts_LP)
dPman_0_LP <- calcNormFactors(dPman_0_LP)
keep_LP <- rowSums(cpm(dPman_0_LP) > 0.5) >= 60
dPman_LP <- dPman_0_LP[keep_LP,]
dim(dPman_LP)

# Compute logCPM matrix for LP data
logCPM_matrix_LP <- cpm(dPman_LP, log = TRUE)

# Combine EP and LP
common_genes <- intersect(rownames(logCPM_matrix_EP), rownames(logCPM_matrix_LP))

# 2. Subset both matrices to the common genes and combine
logCPM_matrix_combined <- cbind(
  logCPM_matrix_EP[common_genes, ],
  logCPM_matrix_LP[common_genes, ]
)


Sample_Info_Combined <- rbind(
  EP_Sample_Info[, c("Sample_ID", "Group_Timepoint")],
  LP_Sample_Info[, c("Sample_ID", "Group_Timepoint")]
)


rm(logCPM_matrix_EP, 
   logCPM_matrix_LP, 
   EP_Sample_Info, 
   LP_Sample_Info,
   Pman_rawreads_EP,
   Pman_rawreads_LP,
   Pman_readcounts_EP,
   Pman_readcounts_LP,
   dPman_0_EP,
   dPman_0_LP,
   dPman_EP,
   dPman_LP,
   keep_EP,
   keep_LP)  # Remove specific objects



# Step 1: Universal LOC renaming dictionary
loc_renames <- c("LOC102904208" = "GAPDH-Like_1",
                 "LOC102916082" = "GAPDH-Like_2",
                 "LOC102923285" = "PKM-Like",
                 "LOC102928417" = "LDHA-Like")



# COMBINED ---------

glycolysis_genes_raw <- c("Aldoa",
                          "Bpgm",
                          "Eno1",
                          "Eno2",
                          "LOC102904208",
                          "LOC102916082",
                          "Gpi",
                          "Hk1",
                          "Hk2",
                          "Hkdc1",
                          "Pfkl",
                          "Pfkm",
                          "Pgam1",
                          "Pgam2",
                          "Pklr",
                          "LOC102923285",
                          "Pkm",
                          "Tpi1",
                          "LOC102928417",
                          "Ldha")

FA_genes_raw = c("Acadl",
                 "Acadm",
                 "Acads",
                 "Acadvl",
                 "Acat1",
                 "Acsl1",
                 "Acsl3",
                 "Acsl4",
                 "Acsl5",
                 "Acss2",
                 "Chkb",
                 "Cpt1a",
                 "Cpt1b",
                 "Cpt2",
                 "Crat",
                 "Decr1",
                 "Dld",
                 "Echs1",
                 "Gcdh",
                 "Gk",
                 "Gpd2",
                 "Hadh",
                 "Hadha",
                 "Hadhb",
                 "Lipc",
                 "Lipe",
                 "Lpl",
                 "Pnpla2",
                 "Slc25a20",
                 "Tpi1")

Solute_genes_raw = c("Mtch1",
                     "Mtch2",
                     "Slc25a1",
                     "Slc25a10",
                     "Slc25a11",
                     "Slc25a12",
                     "Slc25a13",
                     "Slc25a14",
                     "Slc25a15",
                     "Slc25a16",
                     "Slc25a17",
                     "Slc25a18",
                     "Slc25a19",
                     "Slc25a20",
                     "Slc25a21",
                     "Slc25a22",
                     "Slc25a23",
                     "Slc25a25",
                     "Slc25a26",
                     "Slc25a27",
                     "Slc25a28",
                     "Slc25a29",
                     "Slc25a3",
                     "LOC102918924",
                     "LOC121832364",
                     "SLC25A30",
                     "Slc25a32",
                     "Slc25a33",
                     "Slc25a35",
                     "Slc25a36",
                     "Slc25a37",
                     "Slc25a38",
                     "Slc25a39",
                     "Slc25a4",
                     "Slc25a40",
                     "Slc25a42",
                     "Slc25a43",
                     "Slc25a44",
                     "Slc25a45",
                     "Slc25a46",
                     "Slc25a47",
                     "Slc25a5",
                     "LOC102912071",
                     "Slc25a51",
                     "Ucp2")

AA_genes_raw = c("Pdpn",
                 "Sh3bp4",
                 "Slc16a10",
                 "Slc1a4",
                 "Slc1a5",
                 "Slc25a29",
                 "Slc36a1",
                 "Slc36a2",
                 "LOC121830744",
                 "Slc36a4",
                 "Slc38a1",
                 "Slc38a2",
                 "Slc38a3",
                 "Slc38a4",
                 "Slc3a1",
                 "Slc3a2",
                 "Slc43a1",
                 "Slc43a2",
                 "Slc6a1",
                 "Slc6a13",
                 "Slc6a17",
                 "Slc6a2",
                 "Slc6a4",
                 "Slc6a6",
                 "Slc6a8",
                 "LOC102917669",
                 "Slc7a1",
                 "Slc7a14",
                 "Slc7a2",
                 "Slc7a5",
                 "Slc7a6",
                 "Slc7a8",
                 "Slc7a9")

PPP_genes_raw = c("Aldoa",
                  "Dera",
                  "G6pd",
                  "Gpi",
                  "Pfkl",
                  "Pfkm",
                  "Pgd",
                  "Pgls",
                  "Pgm1",
                  "Pgm3",
                  "Prpsap2",
                  "Prpsap1",
                  "Prps2",
                  "Prps1",
                  "Rbks",
                  "Rpe",
                  "Rpia",
                  "Taldo1",
                  "Tkt")


# Step 3: Rename LOCs in gene lists
rename_genes <- function(gene_list, rename_dict) {
  ifelse(gene_list %in% names(rename_dict), rename_dict[gene_list], gene_list)
}

glycolysis_genes <- rename_genes(glycolysis_genes_raw, loc_renames)
FA_genes <- rename_genes(FA_genes_raw, loc_renames)
Solute_genes <- rename_genes(Solute_genes_raw, loc_renames)
AA_genes <- rename_genes(AA_genes_raw, loc_renames)
PPP_genes <- rename_genes(PPP_genes_raw, loc_renames)


# Step 4: Combine all genes and subset expression matrix
all_genes_raw <- unique(c(glycolysis_genes_raw, 
                          FA_genes_raw,
                          Solute_genes_raw,
                          AA_genes_raw,
                          PPP_genes_raw))


subset_expr_matrix <- logCPM_matrix_combined[rownames(logCPM_matrix_combined) %in% all_genes_raw, ]

# Step 5: Rename LOCs in the matrix
rownames(subset_expr_matrix) <- rename_genes(rownames(subset_expr_matrix), loc_renames)

# Step 6: Pathway annotation with renamed gene lists
gene_pathway_annotation <- c(
  setNames(rep("Glycolysis", length(glycolysis_genes)), glycolysis_genes),
  setNames(rep("Fatty Acid Metabolism", length(FA_genes)), FA_genes),
  setNames(rep("Solute Carriers", length(Solute_genes)), Solute_genes),
  setNames(rep("Amino Acid Transport", length(AA_genes)), AA_genes),
  setNames(rep("Pentose Phosphate Pathway", length(PPP_genes)), PPP_genes)
)

gene_annotation_df <- data.frame(
  Pathway = gene_pathway_annotation[rownames(subset_expr_matrix)],
  row.names = rownames(subset_expr_matrix)
)


# Define annotation colors
annotation_colors <- list(
  Pathway = c(
    "Glycolysis" = "#E69F00",
    "Fatty Acid Metabolism" = "#009E73",
    "Solute Carriers" = "lightpink",
    "Amino Acid Transport" = "gray",
    "Pentose Phosphate Pathway" = "goldenrod1"
  )
)

# Ensure sample info has correct rownames
rownames(Sample_Info_Combined) <- Sample_Info_Combined$Sample_ID

# Compute group averages
group_avg_expr <- sapply(unique(Sample_Info_Combined$Group_Timepoint), function(group) {
  samples_in_group <- rownames(Sample_Info_Combined)[Sample_Info_Combined$Group_Timepoint == group]
  rowMeans(subset_expr_matrix[, samples_in_group, drop = FALSE])
})

# Format matrix
group_avg_expr <- as.matrix(group_avg_expr)
rownames(group_avg_expr) <- rownames(subset_expr_matrix)

# Reorder columns
desired_order <- c("EP_BW1N", "EP_BW2H", "EP_ME1N", "EP_ME2H", "LP_BW1N", "LP_BW2H", "LP_ME1N", "LP_ME2H")
colnames(group_avg_expr) <- unique(Sample_Info_Combined$Group_Timepoint)
group_avg_expr <- group_avg_expr[, desired_order]

# Plot the heatmap
pheatmap(
  group_avg_expr,
  annotation_row = gene_annotation_df,
  annotation_colors = annotation_colors,
  cluster_cols = FALSE,
  fontsize_row = 5
)


pheatmap(
  mat = group_avg_expr,
  annotation_row = gene_annotation_df,
  annotation_colors = annotation_colors,
  cluster_cols = FALSE,
  cluster_rows = TRUE,  # Enable row clustering if genes are many
  fontsize_row = 8,     # Smaller gene name font
  fontsize_col = 10,    # Column label font size
  cellwidth = 20,       # Adjust based on plot width
  cellheight = 10,      # Adjust to reduce overlap in gene labels
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  show_rownames = TRUE,
  show_colnames = TRUE,
  border_color = NA,    # Cleaner look without borders
  main = "Pathway Gene Expression Heatmap"
)





# GLYCOLYSIS ONLY ----  

# Step 1: Universal LOC renaming dictionary
loc_renames <- c("LOC102904208" = "GAPDH",
                 "LOC102916082" = "GAPDH_",
                 "LOC102923285" = "PKM",
                 "LOC102928417" = "LDHA",
                 "LOC102918924" = "Slc25a30_",
                 "LOC121832364" = "Slc25a30__",
                 "LOC102912071" = "Slc25a51_",
                 "LOC121830744" = "Slc36a4_",
                 "LOC102917669" = "Slc6a8")

# Step 2: Define raw gene lists (unrenamed)
glycolysis_genes_raw <- c("Aldoa",
                          "Bpgm",
                          "Eno1",
                          "Eno2",
                          "LOC102904208",
                          "LOC102916082",
                          "Gpi",
                          "Hk1",
                          "Hk2",
                          "Hkdc1",
                          "Pfkl",
                          "Pfkm",
                          "Pgam1",
                          "Pgam2",
                          "Pklr",
                          "LOC102923285",
                          "Pkm",
                          "Tpi1",
                          "LOC102928417",
                          "Ldha")

# Step 3: Rename LOCs in gene lists
rename_genes <- function(gene_list, rename_dict) {
  ifelse(gene_list %in% names(rename_dict), rename_dict[gene_list], gene_list)
}

glycolysis_genes <- rename_genes(glycolysis_genes_raw, loc_renames)

# Step 4: Combine all genes and subset expression matrix
all_genes_raw <- unique(c(glycolysis_genes_raw))
subset_expr_matrix <- logCPM_matrix_combined[rownames(logCPM_matrix_combined) %in% all_genes_raw, ]

# Step 5: Rename LOCs in the matrix
rownames(subset_expr_matrix) <- rename_genes(rownames(subset_expr_matrix), loc_renames)

# Step 6: Pathway annotation with renamed gene lists
gene_pathway_annotation <- c(
  setNames(rep("Glycolysis", length(glycolysis_genes)), glycolysis_genes)
)

gene_annotation_df <- data.frame(
  Pathway = gene_pathway_annotation[rownames(subset_expr_matrix)],
  row.names = rownames(subset_expr_matrix)
)


# Define annotation colors
annotation_colors <- list(
  Pathway = c(
    "Glycolysis" = "#E69F00"
  )
)

# Ensure sample info has correct rownames
rownames(Sample_Info_Combined) <- Sample_Info_Combined$Sample_ID

# Compute group averages
group_avg_expr <- sapply(unique(Sample_Info_Combined$Group_Timepoint), function(group) {
  samples_in_group <- rownames(Sample_Info_Combined)[Sample_Info_Combined$Group_Timepoint == group]
  rowMeans(subset_expr_matrix[, samples_in_group, drop = FALSE])
})

# Format matrix
group_avg_expr <- as.matrix(group_avg_expr)
rownames(group_avg_expr) <- rownames(subset_expr_matrix)

# Reorder columns
desired_order <- c("EP_BW1N", "EP_BW2H", "EP_ME1N", "EP_ME2H", "LP_BW1N", "LP_BW2H", "LP_ME1N", "LP_ME2H")
colnames(group_avg_expr) <- unique(Sample_Info_Combined$Group_Timepoint)
group_avg_expr <- group_avg_expr[, desired_order]

# Plot the heatmap
pheatmap(
  group_avg_expr,
  annotation_row = gene_annotation_df,
  annotation_colors = annotation_colors,
  cluster_cols = FALSE
)


pheatmap(
  group_avg_expr,
  annotation_row = gene_annotation_df,
  annotation_colors = annotation_colors,
  cluster_cols = FALSE,
  filename = "Dream_Output/Heatmaps/Glycolysis.png",   # Save as PNG
  width = 15,                        
  height = 12,
  fontsize = 16,                  
  fontsize_row = 25,             
  fontsize_col = 14)



# FATTY ACID METABOLISM ONLY -------

FA_genes_raw = c("Acadl",
                 "Acadm",
                 "Acads",
                 "Acadvl",
                 "Acat1",
                 "Acsl1",
                 "Acsl3",
                 "Acsl4",
                 "Acsl5",
                 "Acss2",
                 "Chkb",
                 "Cpt1a",
                 "Cpt1b",
                 "Cpt2",
                 "Crat",
                 "Decr1",
                 "Dld",
                 "Echs1",
                 "Gcdh",
                 "Gk",
                 "Gpd2",
                 "Hadh",
                 "Hadha",
                 "Hadhb",
                 "Lipc",
                 "Lipe",
                 "Lpl",
                 "Pnpla2",
                 "Slc25a20",
                 "Tpi1")


# Step 3: Rename LOCs in gene lists
rename_genes <- function(gene_list, rename_dict) {
  ifelse(gene_list %in% names(rename_dict), rename_dict[gene_list], gene_list)
}

FA_genes <- rename_genes(FA_genes_raw, loc_renames)

# Step 4: Combine all genes and subset expression matrix
all_genes_raw <- unique(c(FA_genes_raw))
subset_expr_matrix <- logCPM_matrix_combined[rownames(logCPM_matrix_combined) %in% all_genes_raw, ]

# Step 5: Rename LOCs in the matrix
rownames(subset_expr_matrix) <- rename_genes(rownames(subset_expr_matrix), loc_renames)

# Step 6: Pathway annotation with renamed gene lists
gene_pathway_annotation <- c(
  setNames(rep("Fatty Acid Metabolism", length(FA_genes)), FA_genes)
)

gene_annotation_df <- data.frame(
  Pathway = gene_pathway_annotation[rownames(subset_expr_matrix)],
  row.names = rownames(subset_expr_matrix)
)


# Define annotation colors
annotation_colors <- list(
  Pathway = c(
    "Fatty Acid Metabolism" = "#009E73"
  )
)

# Ensure sample info has correct rownames
rownames(Sample_Info_Combined) <- Sample_Info_Combined$Sample_ID

# Compute group averages
group_avg_expr <- sapply(unique(Sample_Info_Combined$Group_Timepoint), function(group) {
  samples_in_group <- rownames(Sample_Info_Combined)[Sample_Info_Combined$Group_Timepoint == group]
  rowMeans(subset_expr_matrix[, samples_in_group, drop = FALSE])
})

# Format matrix
group_avg_expr <- as.matrix(group_avg_expr)
rownames(group_avg_expr) <- rownames(subset_expr_matrix)

# Reorder columns
desired_order <- c("EP_BW1N", "EP_BW2H", "EP_ME1N", "EP_ME2H", "LP_BW1N", "LP_BW2H", "LP_ME1N", "LP_ME2H")
colnames(group_avg_expr) <- unique(Sample_Info_Combined$Group_Timepoint)
group_avg_expr <- group_avg_expr[, desired_order]

# Plot the heatmap
pheatmap(
  group_avg_expr,
  annotation_row = gene_annotation_df,
  annotation_colors = annotation_colors,
  cluster_cols = FALSE
)


pheatmap(
  group_avg_expr,
  annotation_row = gene_annotation_df,
  annotation_colors = annotation_colors,
  cluster_cols = FALSE,
  filename = "Dream_Output/Heatmaps/FattyAcid.png",   # Save as PNG
  width = 15,                        
  height = 12,
  fontsize = 16,                  
  fontsize_row = 25,             
  fontsize_col = 14)


# SOLUTE CARRIERS ---------

Solute_genes_raw = c("Mtch1",
                     "Mtch2",
                     "Slc25a1",
                     "Slc25a10",
                     "Slc25a11",
                     "Slc25a12",
                     "Slc25a13",
                     "Slc25a14",
                     "Slc25a15",
                     "Slc25a16",
                     "Slc25a17",
                     "Slc25a18",
                     "Slc25a19",
                     "Slc25a20",
                     "Slc25a21",
                     "Slc25a22",
                     "Slc25a23",
                     "Slc25a25",
                     "Slc25a26",
                     "Slc25a27",
                     "Slc25a28",
                     "Slc25a29",
                     "Slc25a3",
                     "LOC102918924",
                     "LOC121832364",
                     "SLC25A30",
                     "Slc25a32",
                     "Slc25a33",
                     "Slc25a35",
                     "Slc25a36",
                     "Slc25a37",
                     "Slc25a38",
                     "Slc25a39",
                     "Slc25a4",
                     "Slc25a40",
                     "Slc25a42",
                     "Slc25a43",
                     "Slc25a44",
                     "Slc25a45",
                     "Slc25a46",
                     "Slc25a47",
                     "Slc25a5",
                     "LOC102912071",
                     "Slc25a51",
                     "Ucp2")


# Step 3: Rename LOCs in gene lists
rename_genes <- function(gene_list, rename_dict) {
  ifelse(gene_list %in% names(rename_dict), rename_dict[gene_list], gene_list)
}

Solute_genes <- rename_genes(Solute_genes_raw, loc_renames)

# Step 4: Combine all genes and subset expression matrix
all_genes_raw <- unique(c(Solute_genes_raw))
subset_expr_matrix <- logCPM_matrix_combined[rownames(logCPM_matrix_combined) %in% all_genes_raw, ]

# Step 5: Rename LOCs in the matrix
rownames(subset_expr_matrix) <- rename_genes(rownames(subset_expr_matrix), loc_renames)

# Step 6: Pathway annotation with renamed gene lists
gene_pathway_annotation <- c(
  setNames(rep("Solute Carrier", length(Solute_genes)), Solute_genes)
)

gene_annotation_df <- data.frame(
  Pathway = gene_pathway_annotation[rownames(subset_expr_matrix)],
  row.names = rownames(subset_expr_matrix)
)


# Define annotation colors
annotation_colors <- list(
  Pathway = c(
    "Solute Carrier" = "lightpink"
  )
)

# Ensure sample info has correct rownames
rownames(Sample_Info_Combined) <- Sample_Info_Combined$Sample_ID

# Compute group averages
group_avg_expr <- sapply(unique(Sample_Info_Combined$Group_Timepoint), function(group) {
  samples_in_group <- rownames(Sample_Info_Combined)[Sample_Info_Combined$Group_Timepoint == group]
  rowMeans(subset_expr_matrix[, samples_in_group, drop = FALSE])
})

# Format matrix
group_avg_expr <- as.matrix(group_avg_expr)
rownames(group_avg_expr) <- rownames(subset_expr_matrix)

# Reorder columns
desired_order <- c("EP_BW1N", "EP_BW2H", "EP_ME1N", "EP_ME2H", "LP_BW1N", "LP_BW2H", "LP_ME1N", "LP_ME2H")
colnames(group_avg_expr) <- unique(Sample_Info_Combined$Group_Timepoint)
group_avg_expr <- group_avg_expr[, desired_order]

# Plot the heatmap
pheatmap(
  group_avg_expr,
  annotation_row = gene_annotation_df,
  annotation_colors = annotation_colors,
  cluster_cols = FALSE
)

pheatmap(
  group_avg_expr,
  annotation_row = gene_annotation_df,
  annotation_colors = annotation_colors,
  cluster_cols = FALSE,
  filename = "Dream_Output/Heatmaps/Solute.png",   # Save as PNG
  width = 15,                        
  height = 12,
  fontsize = 16,                  
  fontsize_row = 25,             
  fontsize_col = 14)




# AMINO ACID TRANSPORT ---------

AA_genes_raw = c("Pdpn",
                     "Sh3bp4",
                     "Slc16a10",
                     "Slc1a4",
                     "Slc1a5",
                     "Slc25a29",
                     "Slc36a1",
                     "Slc36a2",
                     "LOC121830744",
                     "Slc36a4",
                     "Slc38a1",
                     "Slc38a2",
                     "Slc38a3",
                     "Slc38a4",
                     "Slc3a1",
                     "Slc3a2",
                     "Slc43a1",
                     "Slc43a2",
                     "Slc6a1",
                     "Slc6a13",
                     "Slc6a17",
                     "Slc6a2",
                     "Slc6a4",
                     "Slc6a6",
                     "Slc6a8",
                     "LOC102917669",
                     "Slc7a1",
                     "Slc7a14",
                     "Slc7a2",
                     "Slc7a5",
                     "Slc7a6",
                     "Slc7a8",
                     "Slc7a9")


# Step 3: Rename LOCs in gene lists
rename_genes <- function(gene_list, rename_dict) {
  ifelse(gene_list %in% names(rename_dict), rename_dict[gene_list], gene_list)
}

AA_genes <- rename_genes(AA_genes_raw, loc_renames)

# Step 4: Combine all genes and subset expression matrix
all_genes_raw <- unique(c(AA_genes_raw))
subset_expr_matrix <- logCPM_matrix_combined[rownames(logCPM_matrix_combined) %in% all_genes_raw, ]

# Step 5: Rename LOCs in the matrix
rownames(subset_expr_matrix) <- rename_genes(rownames(subset_expr_matrix), loc_renames)

# Step 6: Pathway annotation with renamed gene lists
gene_pathway_annotation <- c(
  setNames(rep("Amino Acid Transport", length(AA_genes)), AA_genes)
)

gene_annotation_df <- data.frame(
  Pathway = gene_pathway_annotation[rownames(subset_expr_matrix)],
  row.names = rownames(subset_expr_matrix)
)


# Define annotation colors
annotation_colors <- list(
  Pathway = c(
    "Amino Acid Transport" = "gray"
  )
)

# Ensure sample info has correct rownames
rownames(Sample_Info_Combined) <- Sample_Info_Combined$Sample_ID

# Compute group averages
group_avg_expr <- sapply(unique(Sample_Info_Combined$Group_Timepoint), function(group) {
  samples_in_group <- rownames(Sample_Info_Combined)[Sample_Info_Combined$Group_Timepoint == group]
  rowMeans(subset_expr_matrix[, samples_in_group, drop = FALSE])
})

# Format matrix
group_avg_expr <- as.matrix(group_avg_expr)
rownames(group_avg_expr) <- rownames(subset_expr_matrix)

# Reorder columns
desired_order <- c("EP_BW1N", "EP_BW2H", "EP_ME1N", "EP_ME2H", "LP_BW1N", "LP_BW2H", "LP_ME1N", "LP_ME2H")
colnames(group_avg_expr) <- unique(Sample_Info_Combined$Group_Timepoint)
group_avg_expr <- group_avg_expr[, desired_order]

# Plot the heatmap
pheatmap(
  group_avg_expr,
  annotation_row = gene_annotation_df,
  annotation_colors = annotation_colors,
  cluster_cols = FALSE
)

pheatmap(
  group_avg_expr,
  annotation_row = gene_annotation_df,
  annotation_colors = annotation_colors,
  cluster_cols = FALSE,
  filename = "Dream_Output/Heatmaps/AminoAcid.png",   # Save as PNG
  width = 15,                        
  height = 12,
  fontsize = 16,                  
  fontsize_row = 25,             
  fontsize_col = 14)


# PENTOSE PHOSPHATE PATHWAY ---------

PPP_genes_raw = c("Aldoa",
                 "Dera",
                 "G6pd",
                 "Gpi",
                 "Pfkl",
                 "Pfkm",
                 "Pgd",
                 "Pgls",
                 "Pgm1",
                 "Pgm3",
                 "Prpsap2",
                 "Prpsap1",
                 "Prps2",
                 "Prps1",
                 "Rbks",
                 "Rpe",
                 "Rpia",
                 "Taldo1",
                 "Tkt")

# Step 3: Rename LOCs in gene lists
rename_genes <- function(gene_list, rename_dict) {
  ifelse(gene_list %in% names(rename_dict), rename_dict[gene_list], gene_list)
}

PPP_genes <- rename_genes(PPP_genes_raw, loc_renames)

# Step 4: Combine all genes and subset expression matrix
all_genes_raw <- unique(c(PPP_genes_raw))
subset_expr_matrix <- logCPM_matrix_combined[rownames(logCPM_matrix_combined) %in% all_genes_raw, ]

# Step 5: Rename LOCs in the matrix
rownames(subset_expr_matrix) <- rename_genes(rownames(subset_expr_matrix), loc_renames)

# Step 6: Pathway annotation with renamed gene lists
gene_pathway_annotation <- c(
  setNames(rep("Pentose Phosphate Pathway", length(PPP_genes)), PPP_genes)
)

gene_annotation_df <- data.frame(
  Pathway = gene_pathway_annotation[rownames(subset_expr_matrix)],
  row.names = rownames(subset_expr_matrix)
)


# Define annotation colors
annotation_colors <- list(
  Pathway = c(
    "Pentose Phosphate Pathway" = "goldenrod1"
  )
)

# Ensure sample info has correct rownames
rownames(Sample_Info_Combined) <- Sample_Info_Combined$Sample_ID

# Compute group averages
group_avg_expr <- sapply(unique(Sample_Info_Combined$Group_Timepoint), function(group) {
  samples_in_group <- rownames(Sample_Info_Combined)[Sample_Info_Combined$Group_Timepoint == group]
  rowMeans(subset_expr_matrix[, samples_in_group, drop = FALSE])
})

# Format matrix
group_avg_expr <- as.matrix(group_avg_expr)
rownames(group_avg_expr) <- rownames(subset_expr_matrix)

# Reorder columns
desired_order <- c("EP_BW1N", "EP_BW2H", "EP_ME1N", "EP_ME2H", "LP_BW1N", "LP_BW2H", "LP_ME1N", "LP_ME2H")
colnames(group_avg_expr) <- unique(Sample_Info_Combined$Group_Timepoint)
group_avg_expr <- group_avg_expr[, desired_order]

# Plot the heatmap
pheatmap(
  group_avg_expr,
  annotation_row = gene_annotation_df,
  annotation_colors = annotation_colors,
  cluster_cols = FALSE
)

pheatmap(
  group_avg_expr,
  annotation_row = gene_annotation_df,
  annotation_colors = annotation_colors,
  cluster_cols = FALSE,
  filename = "Dream_Output/Heatmaps/PentosePathway.png",   # Save as PNG
  width = 15,                        
  height = 12,
  fontsize = 16,                  
  fontsize_row = 25,             
  fontsize_col = 14)
