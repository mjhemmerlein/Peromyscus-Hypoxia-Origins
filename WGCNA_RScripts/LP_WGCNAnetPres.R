
library(WGCNA)
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

load("LP_WGCNA_Output/BWLP_network.RData")

table(NetBWLP$colors)
moduleLabelsBWLP = NetBWLP$colors
moduleColorsBWLP = labels2colors(NetBWLP$colors)

MEsBW_LP = NetBWLP$MEs # These are module eigengenes:first principal component of the expression matrix, 
# i.e., weighted average expression profile of a sample in a module
geneTreeBWLP = NetBWLP$dendrograms[[1]]
table(moduleColorsBWLP)
dim(table(moduleColorsBWLP))


## Highlanders WGCNA network build, i.e., ME ####

load("LP_WGCNA_Output/MELP_network.RData")

table(NetMELP$colors)
moduleLabelsMELP = NetMELP$colors
moduleColorsMELP = labels2colors(NetMELP$colors)
MEsMELP = NetMELP$MEs # These are module eigengenes:first principal component of the expression matrix, 
# i.e., weighted average expression profile of a sample in a module

geneTreeMELP = NetMELP$dendrograms[[1]]
table(moduleColorsMELP)
dim(table(moduleColorsMELP))

# Network Preservation analysis ####
# Number of data sets that we work with
nSets = 2;
# Object that will contain the expression data
multiExprLP = list();
multiExprLP[[1]] = list(data = ExprData_BWLP);
multiExprLP[[2]] = list(data = ExprData_MELP);

# Names for the two sets
setLabels = c("Lowlander", "Highlander");

# Important: components of multiExpr must carry identificating names
names(multiExprLP) = setLabels

# Display the dimensions of the expression data (if you are confused by this construct, ignore it):
lapply(multiExprLP, lapply, dim)

# Create an object (list) holding the module labels for each set:
colorList = list(moduleColorsBWLP, moduleColorsMELP)

# Components of the list must be named so that the names can be matched to the names of multiExpr
names(colorList) = setLabels

# Calculation of module preservation statistics
system.time( {
  mpLP = WGCNA::modulePreservation(multiExprLP, colorList,
                                   referenceNetworks = 1,
                                   loadPermutedStatistics = FALSE,
                                   nPermutations = 200,
                                   verbose = 3, 
                                   maxModuleSize = 4619) #set max size to max in lowland network
} )

save(mpLP, file = "LP_WGCNA_Output/LPHighLander-Lowlander_modPreservation.RData")
load(file = "LP_WGCNA_Output/LPHighLander-Lowlander_modPreservation.RData")

library(impute)
# Impute missing data and calculate eigengenes
impExprLP = list();
for (set in 1:nSets)
{
  impExprLP[[set]] = list(data = t(impute.knn(t(multiExprLP[[set]]$data))$data));
}
eigengenes = list();
for (set in 1:nSets)
{
  eigengenes[[set]] = multiSetMEs(impExprLP, universalColors = colorList[[set]], excludeGrey = TRUE);
  for (ss in 1:nSets)
  {
    rownames(eigengenes[[set]][[ss]]$data) = rownames(multiExprLP[[ss]]$data);
  }
}

## Analysis and display of module preservation results

#  Isolate the observed statistics and their Z scores:
# Load the module preservation statistics
ref = 1 # Select the lowland data as reference
test = 2 # Select the highland data as test
statsObsLP = cbind(mpLP$quality$observed[[ref]][[test]][, -1], mpLP$preservation$observed[[ref]][[test]][, -1])
statsZLP = cbind(mpLP$quality$Z[[ref]][[test]][, -1], mpLP$preservation$Z[[ref]][[test]][, -1])
Z.PreservationStatsLP = mpLP$preservation$Z[[ref]][[test]]
Obs.PreservationStatsLP = mpLP$preservation$observed[[ref]][[test]]
# look at the main output: the preservation Zsummary score.
print(signif(statsZLP[, "Zsummary.pres", drop = FALSE],2))
# Compare preservation to quality:
print(signif(statsZLP[, c("Zsummary.pres", "Zsummary.qual")], 2))

# look at the main output: the preservation Zsummary score.
print(signif(statsObsLP[, "Zsummary.pres", drop = FALSE],2))
# Compare preservation to quality:
print(signif(statsObsLP[, c("Zsummary.pres", "Zsummary.qual")], 2))
statsObsLP
statsZLP


## Visualize Data ####
##  Caption: Preservation of highland placenta modules in lowland placenta. 
#The right panel shows the composite statistic medianRank (Eq.9.20) versus module size. 
#he higher the medianRank the less preserved is the module relative to other modules. 
#Since medianRank is based on the observed preservation statistics (as opposed to Z statistics or p-values) 
#we find that it is much less dLPendent on module size. 
#The upper right panel shows the composite statistic Zsummary (Eq.9.1). 
#If Zsummary> 10 there is strong evidence that the module is preserved (Langfelder et al 2011). 
#If Zsummary<2, there is no evidence that the module preserved. 
#Note that Zsummary shows a strong dLPendence on module size. 
##
modColorsLP = rownames(Obs.PreservationStatsLP)
moduleSizeLP = Obs.PreservationStatsLP$moduleSize
# we will omit the grey module (background genes)
# and the gold module (random sample of genes)
selectModules = !(modColorsLP %in% c("grey", "gold"))
# Text labels for points
point.labelLP = modColorsLP[selectModules]

#Composite preservation statistics
medianRankLP = Obs.PreservationStatsLP$medianRank.pres
ZsummaryLP = Z.PreservationStatsLP$Zsummary.pres

stats = as.data.frame((signif(statsZLP[, c("Zsummary.pres", "Zsummary.qual")], 2)))

### medianRanks and Zsummary Plot ####
par(mfrow=c(1,2),mar = c(4.5,4.5,2.5,1))
# plot medianRank versus module size: The is useful for comparing relative 
# preservation among multiple modules: a module with lower median rank tends to 
# exhibit stronger observed preservation statistics than a module with a higher 
# median rank.
plot(moduleSizeLP[selectModules],medianRankLP[selectModules],col=1,
     bg=modColorsLP[selectModules],pch = 21,main="medianRank Preservation",
     cex = 1.5, ylab ="medianRank",xlab="Module size", log="x")
labelPoints(moduleSizeLP[selectModules],medianRankLP[selectModules],point.labelLP,cex=1,offs=0.03)

# plot Zsummary versus module size: If Zsummary is greater than 10 it means the module is preserved
# if 2<Zsummary<10 there is weak to moderate preservation, if less than 2
# there is no evidence
plot(moduleSizeLP[selectModules],ZsummaryLP[selectModules], col = 1,
     bg=modColorsLP[selectModules],pch = 21,main="LP Zsummary Preservation",
     cex=1.5,ylab ="Zsummary", xlab = "Module size", log = "x", ylim = c(0,100))
labelPoints(moduleSizeLP[selectModules],ZsummaryLP[selectModules],point.labelLP,cex=1,offs=0.01)
# Add threshold lines for Zsummary
abline(h=0); abline(h=2, col = "blue", lty = 2); abline(h=10, col = "red", lty = 2)


#### Add Zsummary scores to BWtables (i.e., preservation stats for highlanders)
BWpreservedStatsLP = print(signif(statsZLP[, "Zsummary.pres", drop = FALSE],2))
