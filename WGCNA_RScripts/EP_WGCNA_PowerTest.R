
library(WGCNA, exclude = "GO.db")
library(dplyr)
library(edgeR)
library(lme4)
library(lmerTest)
library(readxl)
library(ggplot2)

# Read in sample info
Traits = read_xlsx("RNA_Seq_RawData/MetaData_EP.xlsx")
Traits = as.data.frame(Traits)
Traits = Traits %>%
  filter(Strain != "LL")

rownames(Traits) = Traits$Sample_ID

Traits$Sample_ID = as.factor(Traits$Sample_ID)
Traits$Lane = as.factor(Traits$Lane)
Traits$TubeID = as.factor(Traits$TubeID)
Traits$Mom = as.factor(Traits$Mom)
Traits$Strain = as.factor(Traits$Strain)
Traits$O2 = as.factor(Traits$O2)
Traits$Sex = as.factor(Traits$Sex)
Traits$Group = as.factor(Traits$Group)

Traits = Traits %>%
  select(c(-Order, -Set, -Sex_Group, -`Seq Read FWD`, -`Seq Read REV`))

TraitsBW_EP = Traits %>%
  filter(Strain == "BW")

TraitsME_EP = Traits %>%
  filter(Strain == "ME")



# Read in Files + QC
EPcounts = read_xlsx("RNA_Seq_RawData/EP_Pman_ExtMMFrac_readcounts_Exon.xlsx")
EPcounts = as.data.frame(EPcounts)
EPcounts = `row.names<-`(EPcounts, EPcounts$Geneid)
EPcounts <- EPcounts[,-c(1:6)]


# Check read table vs sample info
Check = Traits$Seq_Name
colnames(EPcounts) == Check

colnames(EPcounts) = rownames(Traits)

dPman_0 <- DGEList(EPcounts)
dPman_0 <- calcNormFactors(dPman_0)
dim(dPman_0)

keep = rowSums(cpm(dPman_0) > 0.5 ) >= 60
dPman = dPman_0[keep,]
dim(dPman)

rownames(Traits) == colnames(dPman)

# WGCNA Analysis
ExprData_EP = cpm(dPman, log=TRUE, prior.count=2, normalized.lib.sizes=TRUE) #cpm normalized and log transformed expression dataset

# Filter ME and BW
ExprData_MEEP = ExprData_EP[,colnames(ExprData_EP) %in% rownames(TraitsME_EP)]
ExprData_BWEP = ExprData_EP[,colnames(ExprData_EP) %in% rownames(TraitsBW_EP)]

# Transpose expression data for further analysis
# ME
ExprData_MEEP = as.data.frame(t(ExprData_MEEP))
ExprData_MEEP$treatment = TraitsME_EP$Strain
ExprData_MEEP = ExprData_MEEP %>% select(-treatment)

# BW
ExprData_BWEP = as.data.frame(t(ExprData_BWEP))
ExprData_BWEP$treatment = TraitsBW_EP$Strain
ExprData_BWEP = ExprData_BWEP %>% select(-treatment)

match(rownames(TraitsME_EP), rownames(ExprData_MEEP))# should be in order if datasets are matched
match(rownames(TraitsBW_EP), rownames(ExprData_BWEP))

# WGCNA

# Lowlanders, i.e., BW ####
# Cluster the samples to check for outliers
sampleTreeBW = hclust(dist(ExprData_BWEP), method = "average")
# Plot the sample trees
plot(sampleTreeBW, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2) 

## Set power threshold
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sftBW = pickSoftThreshold(ExprData_BWEP, powerVector = powers, verbose = 5, 
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

# POWER = 12
NetBWEP = blockwiseModules(ExprData_BWEP, 
                           power = 12, 
                           maxBlockSize = 14000,
                           TOMType = "signed", networkType = "signed hybrid",
                           minModuleSize = 30,
                           reassignThreshold = 0, 
                           mergeCutHeight = 0.25,
                           numericLabels = TRUE, 
                           pamRespectsDendro = FALSE,
                           saveTOMs = TRUE,
                           saveTOMFileBase = "EP_WGCNA_Output_Test/EPBWExprTOM",
                           verbose = 3)

load("EP_WGCNA_Output_Test/EPBWExprTOM-block.1.RData"); load("EP_WGCNA_Output_Test/EPBWExprTOM-block.2.RData")
table(NetBWEP$colors)
moduleLabelsBWEP = NetBWEP$colors
moduleColorsBWEP = labels2colors(NetBWEP$colors)

MEsBW_EP = NetBWEP$MEs # These are module eigengenes:first principal component of the expression matrix, 
# i.e., weighted average expression profile of a sample in a module
geneTreeBWEP = NetBWEP$dendrograms[[1]]
table(moduleColorsBWEP)
dim(table(moduleColorsBWEP))

# POWER = 16
NetBWEP = blockwiseModules(ExprData_BWEP, 
                           power = 16, 
                           maxBlockSize = 14000,
                           TOMType = "signed", networkType = "signed hybrid",
                           minModuleSize = 30,
                           reassignThreshold = 0, 
                           mergeCutHeight = 0.25,
                           numericLabels = TRUE, 
                           pamRespectsDendro = FALSE,
                           saveTOMs = TRUE,
                           saveTOMFileBase = "EP_WGCNA_Output_Test/EPBWExprTOM",
                           verbose = 3)

load("EP_WGCNA_Output_Test/EPBWExprTOM-block.1.RData"); load("EP_WGCNA_Output_Test/EPBWExprTOM-block.2.RData")
table(NetBWEP$colors)
moduleLabelsBWEP = NetBWEP$colors
moduleColorsBWEP = labels2colors(NetBWEP$colors)

MEsBW_EP = NetBWEP$MEs # These are module eigengenes:first principal component of the expression matrix, 
# i.e., weighted average expression profile of a sample in a module
geneTreeBWEP = NetBWEP$dendrograms[[1]]
table(moduleColorsBWEP)
dim(table(moduleColorsBWEP))

# POWER = 18
NetBWEP = blockwiseModules(ExprData_BWEP, 
                           power = 18, 
                           maxBlockSize = 14000,
                           TOMType = "signed", networkType = "signed hybrid",
                           minModuleSize = 30,
                           reassignThreshold = 0, 
                           mergeCutHeight = 0.25,
                           numericLabels = TRUE, 
                           pamRespectsDendro = FALSE,
                           saveTOMs = TRUE,
                           saveTOMFileBase = "EP_WGCNA_Output_Test/EPBWExprTOM",
                           verbose = 3)

load("EP_WGCNA_Output_Test/EPBWExprTOM-block.1.RData"); load("EP_WGCNA_Output_Test/EPBWExprTOM-block.2.RData")
table(NetBWEP$colors)
moduleLabelsBWEP = NetBWEP$colors
moduleColorsBWEP = labels2colors(NetBWEP$colors)

MEsBW_EP = NetBWEP$MEs # These are module eigengenes:first principal component of the expression matrix, 
# i.e., weighted average expression profile of a sample in a module
geneTreeBWEP = NetBWEP$dendrograms[[1]]
table(moduleColorsBWEP)
dim(table(moduleColorsBWEP))


# POWER = 20
NetBWEP = blockwiseModules(ExprData_BWEP, 
                           power = 20, 
                           maxBlockSize = 14000,
                           TOMType = "signed", networkType = "signed hybrid",
                           minModuleSize = 30,
                           reassignThreshold = 0, 
                           mergeCutHeight = 0.25,
                           numericLabels = TRUE, 
                           pamRespectsDendro = FALSE,
                           saveTOMs = TRUE,
                           saveTOMFileBase = "EP_WGCNA_Output_Test/EPBWExprTOM",
                           verbose = 3)

load("EP_WGCNA_Output_Test/EPBWExprTOM-block.1.RData"); load("EP_WGCNA_Output_Test/EPBWExprTOM-block.2.RData")
table(NetBWEP$colors)
moduleLabelsBWEP = NetBWEP$colors
moduleColorsBWEP = labels2colors(NetBWEP$colors)

MEsBW_EP = NetBWEP$MEs # These are module eigengenes:first principal component of the expression matrix, 
# i.e., weighted average expression profile of a sample in a module
geneTreeBWEP = NetBWEP$dendrograms[[1]]
table(moduleColorsBWEP)
dim(table(moduleColorsBWEP))






# Highlanders WGCNA network build, i.e., ME ####

#cluster the samples to check for outliers
sampleTreeMEEP = hclust(dist(ExprData_MEEP), method = "average")
# Plot the sample trees
plot(sampleTreeMEEP, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2) 

## Set power threshold
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sftME = pickSoftThreshold(ExprData_MEEP, powerVector = powers, verbose = 5, 
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
NetMEEP = blockwiseModules(ExprData_MEEP, 
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
                           saveTOMFileBase = "EP_WGCNA_Output_Test/EPMEExprTOM",
                           verbose = 3)

load("EP_WGCNA_Output_Test/EPMEExprTOM-block.1.RData"); load("EP_WGCNA_Output_Test/EPMEExprTOM-block.2.RData")
table(NetMEEP$colors)
moduleLabelsMEEP = NetMEEP$colors
moduleColorsMEEP = labels2colors(NetMEEP$colors)
MEsMEEP = NetMEEP$MEs # These are module eigengenes:first principal component of the expression matrix, 
# i.e., weighted average expression profile of a sample in a module

geneTreeMEEP = NetMEEP$dendrograms[[1]]
table(moduleColorsMEEP)
dim(table(moduleColorsMEEP))

# POWER = 8
NetMEEP = blockwiseModules(ExprData_MEEP, 
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
                           saveTOMFileBase = "EP_WGCNA_Output_Test/EPMEExprTOM",
                           verbose = 3)

load("EP_WGCNA_Output_Test/EPMEExprTOM-block.1.RData"); load("EP_WGCNA_Output_Test/EPMEExprTOM-block.2.RData")
table(NetMEEP$colors)
moduleLabelsMEEP = NetMEEP$colors
moduleColorsMEEP = labels2colors(NetMEEP$colors)
MEsMEEP = NetMEEP$MEs # These are module eigengenes:first principal component of the expression matrix, 
# i.e., weighted average expression profile of a sample in a module

geneTreeMEEP = NetMEEP$dendrograms[[1]]
table(moduleColorsMEEP)
dim(table(moduleColorsMEEP))

# POWER = 12
NetMEEP = blockwiseModules(ExprData_MEEP, 
                           power = 12, 
                           maxBlockSize = 14000,
                           TOMType = "signed", 
                           networkType = "signed hybrid",
                           minModuleSize = 30,
                           reassignThreshold = 0, 
                           mergeCutHeight = 0.25,
                           numericLabels = TRUE, 
                           pamRespectsDendro = FALSE,
                           saveTOMs = TRUE,
                           saveTOMFileBase = "EP_WGCNA_Output_Test/EPMEExprTOM",
                           verbose = 3)

load("EP_WGCNA_Output_Test/EPMEExprTOM-block.1.RData"); load("EP_WGCNA_Output_Test/EPMEExprTOM-block.2.RData")
table(NetMEEP$colors)
moduleLabelsMEEP = NetMEEP$colors
moduleColorsMEEP = labels2colors(NetMEEP$colors)
MEsMEEP = NetMEEP$MEs # These are module eigengenes:first principal component of the expression matrix, 
# i.e., weighted average expression profile of a sample in a module

geneTreeMEEP = NetMEEP$dendrograms[[1]]
table(moduleColorsMEEP)
dim(table(moduleColorsMEEP))

# POWER = 16
NetMEEP = blockwiseModules(ExprData_MEEP, 
                           power = 16, 
                           maxBlockSize = 14000,
                           TOMType = "signed", 
                           networkType = "signed hybrid",
                           minModuleSize = 30,
                           reassignThreshold = 0, 
                           mergeCutHeight = 0.25,
                           numericLabels = TRUE, 
                           pamRespectsDendro = FALSE,
                           saveTOMs = TRUE,
                           saveTOMFileBase = "EP_WGCNA_Output_Test/EPMEExprTOM",
                           verbose = 3)

load("EP_WGCNA_Output_Test/EPMEExprTOM-block.1.RData"); load("EP_WGCNA_Output_Test/EPMEExprTOM-block.2.RData")
table(NetMEEP$colors)
moduleLabelsMEEP = NetMEEP$colors
moduleColorsMEEP = labels2colors(NetMEEP$colors)
MEsMEEP = NetMEEP$MEs # These are module eigengenes:first principal component of the expression matrix, 
# i.e., weighted average expression profile of a sample in a module

geneTreeMEEP = NetMEEP$dendrograms[[1]]
table(moduleColorsMEEP)
dim(table(moduleColorsMEEP))
