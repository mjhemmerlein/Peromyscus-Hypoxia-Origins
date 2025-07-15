
library(WGCNA, exclude = "GO.db")
library(dplyr)
library(edgeR)
library(lme4)
library(lmerTest)
library(readxl)
library(ggplot2)

# Read in sample info
Traits = read_xlsx("RNA_Seq_RawData/MetaData_LP.xlsx")
Traits = as.data.frame(Traits)
Traits = Traits %>%
  filter(Strain != "LL")
Traits = Traits %>%
  select(-c("Sample_ID_JZ", "JZ_Lane"))
Traits = Traits[-c(80,81),]

rownames(Traits) = Traits$Sample_ID_LZ

Traits$Sample_ID_LZ = as.factor(Traits$Sample_ID_LZ)
Traits$LZ_Lane = as.factor(Traits$LZ_Lane)
Traits$TubeID = as.factor(Traits$TubeID)
Traits$Mom = as.factor(Traits$Mom)
Traits$Strain = as.factor(Traits$Strain)
Traits$O2 = as.factor(Traits$O2)
Traits$Sex = as.factor(Traits$Sex)
Traits$Group = as.factor(Traits$Group)

Traits = Traits %>%
  select(c(-Sex_Group, -JZ_Seq_Name))

Traits = Traits[-which(rownames(Traits) == "LZ089"),]

TraitsBW_LP = Traits %>%
  filter(Strain == "BW")

TraitsME_LP = Traits %>%
  filter(Strain == "ME")


# Read in Files + QC
LPcounts = read_xlsx("RNA_Seq_RawData/LP_Pman_ExtMMFrac_readcounts_Exon.xlsx")
LPcounts = as.data.frame(LPcounts)
LPcounts = `row.names<-`(LPcounts, LPcounts$Geneid)
LPcounts = LPcounts[,-c(1:6)]
LPcounts = subset(LPcounts, select = -c(RNA201216ZC_LZ089_S16_L001_fastp_pman_Halign_liberal.bam))


# Check read table vs sample info
Check = Traits$Seq_Name
colnames(LPcounts) == Check

colnames(LPcounts) = rownames(Traits)

dPman_0 <- DGEList(LPcounts)
dPman_0 <- calcNormFactors(dPman_0)
dim(dPman_0)

keep = rowSums(cpm(dPman_0) > 0.5 ) >= 60
dPman = dPman_0[keep,]
dim(dPman)

rownames(Traits) == colnames(dPman)

# WGCNA Analysis
ExprData_LP = cpm(dPman, log=TRUE, prior.count=2, normalized.lib.sizes=TRUE) #cpm normalized and log transformed expression dataset

# Filter ME and BW
ExprData_MELP = ExprData_LP[,colnames(ExprData_LP) %in% rownames(TraitsME_LP)]
ExprData_BWLP = ExprData_LP[,colnames(ExprData_LP) %in% rownames(TraitsBW_LP)]

# Transpose expression data for further analysis
# ME
ExprData_MELP = as.data.frame(t(ExprData_MELP))
ExprData_MELP$treatment = TraitsME_LP$Strain
ExprData_MELP = ExprData_MELP %>% select(-treatment)

# BW
ExprData_BWLP = as.data.frame(t(ExprData_BWLP))
ExprData_BWLP$treatment = TraitsBW_LP$Strain
ExprData_BWLP = ExprData_BWLP %>% select(-treatment)

match(rownames(TraitsME_LP), rownames(ExprData_MELP))# should be in order if datasets are matched
match(rownames(TraitsBW_LP), rownames(ExprData_BWLP))

# WGCNA

# Lowlanders, i.e., BW ####
# Cluster the samples to check for outliers
sampleTreeBW = hclust(dist(ExprData_BWLP), method = "average")
# Plot the sample trees
plot(sampleTreeBW, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2) 

## Set power threshold
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sftBW = pickSoftThreshold(ExprData_BWLP, powerVector = powers, verbose = 5, 
                          networkType = "signed hybrid") 
# Plot the results:
#pdf('wgcna/rv_beta_plot.pdf', h=4, w=7)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sftBW$fitIndices[,1], 
     -sign(sftBW$fitIndices[,3]) * sftBW$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sftBW$fitIndices[,1], -sign(sftBW$fitIndices[,3])*sftBW$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

# this line corresponds to using an R^2 cut-off of h
abline(h=0.8,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sftBW$fitIndices[,1], 
     sftBW$fitIndices[,5],
     xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sftBW$fitIndices[,1], sftBW$fitIndices[,5], labels=powers, cex=cex1,col="red")


# Call the network topology analysis function

# POWER = 4
NetBWLP = blockwiseModules(ExprData_BWLP, 
                           power = 4, 
                           maxBlockSize = 14000,
                           TOMType = "signed", networkType = "signed hybrid",
                           minModuleSize = 30,
                           reassignThreshold = 0, 
                           mergeCutHeight = 0.25,
                           numericLabels = TRUE, 
                           pamRespectsDendro = FALSE,
                           saveTOMs = TRUE,
                           saveTOMFileBase = "LP_WGCNA_Output_Test/LPBWExprTOM",
                           verbose = 3)

load("LP_WGCNA_Output_Test/LPBWExprTOM-block.1.RData"); load("LP_WGCNA_Output_Test/LPBWExprTOM-block.2.RData")
table(NetBWLP$colors)
moduleLabelsBWLP = NetBWLP$colors
moduleColorsBWLP = labels2colors(NetBWLP$colors)

MEsBW_LP = NetBWLP$MEs # These are module eigengenes:first principal component of the expression matrix, 
# i.e., weighted average expression profile of a sample in a module
geneTreeBWLP = NetBWLP$dendrograms[[1]]
table(moduleColorsBWLP)
dim(table(moduleColorsBWLP))

# POWER = 6
NetBWLP = blockwiseModules(ExprData_BWLP, 
                           power = 6, 
                           maxBlockSize = 14000,
                           TOMType = "signed", networkType = "signed hybrid",
                           minModuleSize = 30,
                           reassignThreshold = 0, 
                           mergeCutHeight = 0.25,
                           numericLabels = TRUE, 
                           pamRespectsDendro = FALSE,
                           saveTOMs = TRUE,
                           saveTOMFileBase = "LP_WGCNA_Output_Test/LPBWExprTOM",
                           verbose = 3)

load("LP_WGCNA_Output_Test/LPBWExprTOM-block.1.RData"); load("LP_WGCNA_Output_Test/LPBWExprTOM-block.2.RData")
table(NetBWLP$colors)
moduleLabelsBWLP = NetBWLP$colors
moduleColorsBWLP = labels2colors(NetBWLP$colors)

MEsBW_LP = NetBWLP$MEs # These are module eigengenes:first principal component of the expression matrix, 
# i.e., weighted average expression profile of a sample in a module
geneTreeBWLP = NetBWLP$dendrograms[[1]]
table(moduleColorsBWLP)
dim(table(moduleColorsBWLP))

# POWER = 7
NetBWLP = blockwiseModules(ExprData_BWLP, 
                           power = 7, 
                           maxBlockSize = 14000,
                           TOMType = "signed", networkType = "signed hybrid",
                           minModuleSize = 30,
                           reassignThreshold = 0, 
                           mergeCutHeight = 0.25,
                           numericLabels = TRUE, 
                           pamRespectsDendro = FALSE,
                           saveTOMs = TRUE,
                           saveTOMFileBase = "LP_WGCNA_Output_Test/LPBWExprTOM",
                           verbose = 3)

load("LP_WGCNA_Output_Test/LPBWExprTOM-block.1.RData"); load("LP_WGCNA_Output_Test/LPBWExprTOM-block.2.RData")
table(NetBWLP$colors)
moduleLabelsBWLP = NetBWLP$colors
moduleColorsBWLP = labels2colors(NetBWLP$colors)

MEsBW_LP = NetBWLP$MEs # These are module eigengenes:first principal component of the expression matrix, 
# i.e., weighted average expression profile of a sample in a module
geneTreeBWLP = NetBWLP$dendrograms[[1]]
table(moduleColorsBWLP)
dim(table(moduleColorsBWLP))


# POWER = 8
NetBWLP = blockwiseModules(ExprData_BWLP, 
                           power = 8, 
                           maxBlockSize = 14000,
                           TOMType = "signed", networkType = "signed hybrid",
                           minModuleSize = 30,
                           reassignThreshold = 0, 
                           mergeCutHeight = 0.25,
                           numericLabels = TRUE, 
                           pamRespectsDendro = FALSE,
                           saveTOMs = TRUE,
                           saveTOMFileBase = "LP_WGCNA_Output_Test/LPBWExprTOM",
                           verbose = 3)

load("LP_WGCNA_Output_Test/LPBWExprTOM-block.1.RData"); load("LP_WGCNA_Output_Test/LPBWExprTOM-block.2.RData")
table(NetBWLP$colors)
moduleLabelsBWLP = NetBWLP$colors
moduleColorsBWLP = labels2colors(NetBWLP$colors)

MEsBW_LP = NetBWLP$MEs # These are module eigengenes:first principal component of the expression matrix, 
# i.e., weighted average expression profile of a sample in a module
geneTreeBWLP = NetBWLP$dendrograms[[1]]
table(moduleColorsBWLP)
dim(table(moduleColorsBWLP))






## Highlanders WGCNA network build, i.e., ME ####

#cluster the samples to check for outliers
sampleTreeMELP = hclust(dist(ExprData_MELP), method = "average")
# Plot the sample trees
plot(sampleTreeMELP, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2) 

## Set power threshold
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sftME = pickSoftThreshold(ExprData_MELP, powerVector = powers, verbose = 5, 
                          networkType = "signed hybrid") 
# Plot the results:
#pdf('wgcna/rv_beta_plot.pdf', h=4, w=7)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sftME$fitIndices[,1], -sign(sftME$fitIndices[,3])*sftME$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sftME$fitIndices[,1], -sign(sftME$fitIndices[,3])*sftME$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

# this line corresponds to using an R^2 cut-off of h
abline(h=0.9,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sftME$fitIndices[,1], sftME$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sftME$fitIndices[,1], sftME$fitIndices[,5], labels=powers, cex=cex1,col="red")

# Call the network topology analysis function

# POWER = 4
NetMELP = blockwiseModules(ExprData_MELP, 
                           power = 4, 
                           maxBlockSize = 14000,
                           TOMType = "signed", 
                           networkType = "signed hybrid",
                           minModuleSize = 30,
                           reassignThreshold = 0, 
                           mergeCutHeight = 0.25,
                           numericLabels = TRUE, 
                           pamRespectsDendro = FALSE,
                           saveTOMs = TRUE,
                           saveTOMFileBase = "LP_WGCNA_Output_Test/LPMEExprTOM",
                           verbose = 3)

load("LP_WGCNA_Output_Test/LPMEExprTOM-block.1.RData"); load("LP_WGCNA_Output_Test/LPMEExprTOM-block.2.RData")
table(NetMELP$colors)
moduleLabelsMELP = NetMELP$colors
moduleColorsMELP = labels2colors(NetMELP$colors)
MEsMELP = NetMELP$MEs # These are module eigengenes:first principal component of the expression matrix, 
# i.e., weighted average expression profile of a sample in a module

geneTreeMELP = NetMELP$dendrograms[[1]]
table(moduleColorsMELP)
dim(table(moduleColorsMELP))

# POWER = 6
NetMELP = blockwiseModules(ExprData_MELP, 
                           power = 6, 
                           maxBlockSize = 14000,
                           TOMType = "signed", 
                           networkType = "signed hybrid",
                           minModuleSize = 30,
                           reassignThreshold = 0, 
                           mergeCutHeight = 0.25,
                           numericLabels = TRUE, 
                           pamRespectsDendro = FALSE,
                           saveTOMs = TRUE,
                           saveTOMFileBase = "LP_WGCNA_Output_Test/LPMEExprTOM",
                           verbose = 3)

load("LP_WGCNA_Output_Test/LPMEExprTOM-block.1.RData"); load("LP_WGCNA_Output_Test/LPMEExprTOM-block.2.RData")
table(NetMELP$colors)
moduleLabelsMELP = NetMELP$colors
moduleColorsMELP = labels2colors(NetMELP$colors)
MEsMELP = NetMELP$MEs # These are module eigengenes:first principal component of the expression matrix, 
# i.e., weighted average expression profile of a sample in a module

geneTreeMELP = NetMELP$dendrograms[[1]]
table(moduleColorsMELP)
dim(table(moduleColorsMELP))

# POWER = 7
NetMELP = blockwiseModules(ExprData_MELP, 
                           power = 7, 
                           maxBlockSize = 14000,
                           TOMType = "signed", 
                           networkType = "signed hybrid",
                           minModuleSize = 30,
                           reassignThreshold = 0, 
                           mergeCutHeight = 0.25,
                           numericLabels = TRUE, 
                           pamRespectsDendro = FALSE,
                           saveTOMs = TRUE,
                           saveTOMFileBase = "LP_WGCNA_Output_Test/LPMEExprTOM",
                           verbose = 3)

load("LP_WGCNA_Output_Test/LPMEExprTOM-block.1.RData"); load("LP_WGCNA_Output_Test/LPMEExprTOM-block.2.RData")
table(NetMELP$colors)
moduleLabelsMELP = NetMELP$colors
moduleColorsMELP = labels2colors(NetMELP$colors)
MEsMELP = NetMELP$MEs # These are module eigengenes:first principal component of the expression matrix, 
# i.e., weighted average expression profile of a sample in a module

geneTreeMELP = NetMELP$dendrograms[[1]]
table(moduleColorsMELP)
dim(table(moduleColorsMELP))

# POWER = 8
NetMELP = blockwiseModules(ExprData_MELP, 
                           power = 8, 
                           maxBlockSize = 14000,
                           TOMType = "signed", 
                           networkType = "signed hybrid",
                           minModuleSize = 30,
                           reassignThreshold = 0, 
                           mergeCutHeight = 0.25,
                           numericLabels = TRUE, 
                           pamRespectsDendro = FALSE,
                           saveTOMs = TRUE,
                           saveTOMFileBase = "LP_WGCNA_Output_Test/LPMEExprTOM",
                           verbose = 3)

load("LP_WGCNA_Output_Test/LPMEExprTOM-block.1.RData"); load("LP_WGCNA_Output_Test/LPMEExprTOM-block.2.RData")
table(NetMELP$colors)
moduleLabelsMELP = NetMELP$colors
moduleColorsMELP = labels2colors(NetMELP$colors)
MEsMELP = NetMELP$MEs # These are module eigengenes:first principal component of the expression matrix, 
# i.e., weighted average expression profile of a sample in a module

geneTreeMELP = NetMELP$dendrograms[[1]]
table(moduleColorsMELP)
dim(table(moduleColorsMELP))
