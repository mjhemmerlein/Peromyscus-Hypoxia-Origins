
library(WGCNA)
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
                           saveTOMFileBase = "EP_WGCNA_Output/EPBWExprTOM",
                           verbose = 3)

load("EP_WGCNA_Output/EPBWExprTOM-block.1.RData"); load("EP_WGCNA_Output/EPBWExprTOM-block.2.RData")
table(NetBWEP$colors)
moduleLabelsBWEP = NetBWEP$colors
moduleColorsBWEP = labels2colors(NetBWEP$colors)

MEsBW_EP = NetBWEP$MEs # These are module eigengenes:first principal component of the expression matrix, 
# i.e., weighted average expression profile of a sample in a module
geneTreeBWEP = NetBWEP$dendrograms[[1]]
table(moduleColorsBWEP)
dim(table(moduleColorsBWEP))

plotDendroAndColors(
  NetBWEP$dendrograms[[1]],
  moduleColorsBWEP[NetBWEP$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03, 
  addGuide = TRUE,
  guideHang = 0.05)

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(ExprData_BWEP, moduleColorsBWEP)$eigengenes
MEsBW_EP = orderMEs(MEs0) #the rownames of this dataset are equal to Expr

# write.csv(MEsBW_EP, "EP_WGCNA_Output/EP_BW_moduleEigengenes.csv")
# save(NetBWEP, MEsBW_EP, moduleLabelsBWEP, moduleColorsBWEP, geneTreeBWEP, file = "EP_WGCNA_Output/BWEP_network.RData")

load("EP_WGCNA_Output/BWEP_network.RData")


## BW Trait Module Associations and Correlations ####
ModuleME_Info_BWEP = read.csv("EP_WGCNA_Output/EP_BW_moduleEigengenes.csv")
rownames(ModuleME_Info_BWEP) = ModuleME_Info_BWEP$X
ModuleME_Info_BWEP = ModuleME_Info_BWEP %>% dplyr::select(-X)
ModuleME_Info_BWEP = ModuleME_Info_BWEP %>% dplyr::select(-MEgrey) #remove grey module

# Correlate modules to hypoxia OVERALL w/ mom as random effect
Strain = TraitsBW_EP$Strain
O2 = TraitsBW_EP$O2
FetalMass = TraitsBW_EP$Fetus
Sex = TraitsBW_EP$Sex
# summaryBWEP = data.frame()
all_modules = colnames(ModuleME_Info_BWEP)

summaryBWEP02 = data.frame()

for (p in all_modules) {
  
  test_data = ModuleME_Info_BWEP[, p]
  test_model = lmer(test_data ~ O2 + (1|(TraitsBW_EP$Mom)))
  output = anova(test_model)
  output_line <- data.frame(
    Module = p,
    F.O2 = output$`F value`,
    Pvalue.O2 = output$`Pr(>F)`
  )
  summaryBWEP02 = rbind(summaryBWEP02, output_line)
}

summaryBWEP02$Pvalue.O2corr <- p.adjust(summaryBWEP02$Pvalue.O2, method = 'BH')

# write.csv(summaryBWEP02, "EP_WGCNA_Output/EP_BW_Module_Trait_ModelSummary.csv")

# Create a dataset containing all gene-specific information
genes=names(ExprData_BWEP)
geneInfoBWEP = data.frame(Gene = genes,
                          moduleColor = moduleColorsBWEP)

# write.csv(geneInfoBWEP, "EP_WGCNA_Output/EP_BW_Modules.csv")

## Plot eigenes against hypoxia/normoxia in BW animals ####
ModuleME_Info_BWEP$hypoxia = TraitsBW_EP$O2

#brown
ggplot(data = ModuleME_Info_BWEP, aes(x = ModuleME_Info_BWEP$hypoxia, y = ModuleME_Info_BWEP$MEbrown)) + 
  geom_boxplot(aes(fill = ModuleME_Info_BWEP$hypoxia), outlier.shape = NA) + theme_classic() + 
  geom_jitter(width = 0.25) + 
  scale_fill_manual(values = c("red", "blue")) + 
  guides(fill = guide_legend(title = "Oxygen")) + 
  xlab("Oxygen") + ylab("Brown Module Eigengene")



## Highlanders WGCNA network build, i.e., ME ####

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
NetMEEP = blockwiseModules(ExprData_MEEP, 
                           power = 20, 
                           maxBlockSize = 14000,
                           TOMType = "signed", 
                           networkType = "signed hybrid",
                           minModuleSize = 30,
                           reassignThreshold = 0, 
                           mergeCutHeight = 0.25,
                           numericLabels = TRUE, 
                           pamRespectsDendro = FALSE,
                           saveTOMs = TRUE,
                           saveTOMFileBase = "EP_WGCNA_Output/EPMEExprTOM",
                           verbose = 3)

load("EP_WGCNA_Output/EPMEExprTOM-block.1.RData"); load("EP_WGCNA_Output/EPMEExprTOM-block.2.RData")
table(NetMEEP$colors)
moduleLabelsMEEP = NetMEEP$colors
moduleColorsMEEP = labels2colors(NetMEEP$colors)
MEsMEEP = NetMEEP$MEs # These are module eigengenes:first principal component of the expression matrix, 
# i.e., weighted average expression profile of a sample in a module

geneTreeMEEP = NetMEEP$dendrograms[[1]]
table(moduleColorsMEEP)
dim(table(moduleColorsMEEP))

plotDendroAndColors(
  NetMEEP$dendrograms[[1]],
  moduleColorsMEEP[NetMEEP$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03, 
  addGuide = TRUE,
  guideHang = 0.05)

# Recalculate MEs with color labels
MEs0ME = moduleEigengenes(ExprData_MEEP, moduleColorsMEEP)$eigengenes
MEsMEEP = orderMEs(MEs0ME) #the rownames of this dataset are equal to Expr

#Add metadata to this dataframe
MEsMEEP$hypoxia = TraitsME_EP$O2
MEsMEEP$fetalMass = TraitsME_EP$Fetus

# write.csv(MEsMEEP, "EP_WGCNA_Output/EP_ME_moduleEigengenes.csv")

# save(NetMEEP, MEsMEEP, moduleLabelsMEEP, moduleColorsMEEP, geneTreeMEEP, file = "EP_WGCNA_Output/MEEP_network.RData")
load("EP_WGCNA_Output/MEEP_network.RData")


## ME Trait Module Associations and Correlations ####
ModuleME_Info_MEEP = read.csv("EP_WGCNA_Output/EP_ME_moduleEigengenes.csv")
rownames(ModuleME_Info_MEEP) = ModuleME_Info_MEEP$X
ModuleME_Info_MEEP = ModuleME_Info_MEEP %>% dplyr::select(-X)
ModuleME_Info_MEEP = ModuleME_Info_MEEP %>% dplyr::select(-MEgrey) #remove grey module
ModuleME_Info_MEEP = ModuleME_Info_MEEP %>% select(-hypoxia)

# Correlate modules to hypoxia OVERALL w/ mom as random effect
Strain = TraitsME_EP$Strain
O2 = TraitsME_EP$O2
FetalMass = TraitsME_EP$Fetus
Sex = TraitsME_EP$Sex
# summaryMEEP = data.frame()
all_modules = colnames(ModuleME_Info_MEEP)

summaryMEEP02 = data.frame()

for (p in all_modules) {
  
  test_data = ModuleME_Info_MEEP[, p]
  test_model = lmer(test_data ~ O2 + (1|(TraitsME_EP$Mom)))
  output = anova(test_model)
  output_line <- data.frame(
    Module = p,
    F.O2 = output$`F value`,
    Pvalue.O2 = output$`Pr(>F)`
  )
  summaryMEEP02 = rbind(summaryMEEP02, output_line)
}

summaryMEEP02$Pvalue.O2corr <- p.adjust(summaryMEEP02$Pvalue.O2, method = 'BH')

# write.csv(summaryMEEP02, "EP_WGCNA_Output/EP_ME_Module_Trait_ModelSummary.csv")

# Create a dataset containing all gene-specific information
genes=names(ExprData_MEEP)
geneInfoMEEP = data.frame(Gene = genes,
                          moduleColor = moduleColorsMEEP)
