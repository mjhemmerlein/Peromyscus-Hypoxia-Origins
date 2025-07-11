
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

load("EP_WGCNA_Output/BWEP_network.RData")

table(NetBWEP$colors)
moduleLabelsBWEP = NetBWEP$colors
moduleColorsBWEP = labels2colors(NetBWEP$colors)

MEsBW_EP = NetBWEP$MEs # These are module eigengenes:first principal component of the expression matrix, 
# i.e., weighted average expression profile of a sample in a module
geneTreeBWEP = NetBWEP$dendrograms[[1]]
table(moduleColorsBWEP)
dim(table(moduleColorsBWEP))

## Highlanders WGCNA network build, i.e., ME ####

load("EP_WGCNA_Output/MEEP_network.RData")

table(NetMEEP$colors)
moduleLabelsMEEP = NetMEEP$colors
moduleColorsMEEP = labels2colors(NetMEEP$colors)
MEsMEEP = NetMEEP$MEs # These are module eigengenes:first principal component of the expression matrix, 
# i.e., weighted average expression profile of a sample in a module

geneTreeMEEP = NetMEEP$dendrograms[[1]]
table(moduleColorsMEEP)
dim(table(moduleColorsMEEP))

# Network Preservation analysis ####
# Number of data sets that we work with
nSets = 2;
# Object that will contain the expression data
multiExprEP = list();
multiExprEP[[1]] = list(data = ExprData_BWEP);
multiExprEP[[2]] = list(data = ExprData_MEEP);

# Names for the two sets
setLabels = c("Lowlander", "Highlander");

# Important: components of multiExpr must carry identificating names
names(multiExprEP) = setLabels

# Display the dimensions of the expression data (if you are confused by this construct, ignore it):
lapply(multiExprEP, lapply, dim)

# Create an object (list) holding the module labels for each set:
colorList = list(moduleColorsBWEP, moduleColorsMEEP)

# Components of the list must be named so that the names can be matched to the names of multiExpr
names(colorList) = setLabels

# Calculation of module preservation statistics
system.time( {
  mpEP = WGCNA::modulePreservation(multiExprEP, colorList,
                                   referenceNetworks = 1,
                                   loadPermutedStatistics = FALSE,
                                   nPermutations = 200,
                                   verbose = 3, 
                                   maxModuleSize = 5301) #set max size to max in lowland network
} )

save(mpEP, file = "EP_WGCNA_Output/EPHighLander-Lowlander_modPreservation.RData")

load(file = "EP_WGCNA_Output/EPHighLander-Lowlander_modPreservation.RData")

library(impute)
# Impute missing data and calculate eigengenes
impExprEP = list();
for (set in 1:nSets)
{
  impExprEP[[set]] = list(data = t(impute.knn(t(multiExprEP[[set]]$data))$data));
}
eigengenes = list();
for (set in 1:nSets)
{
  eigengenes[[set]] = multiSetMEs(impExprEP, universalColors = colorList[[set]], excludeGrey = TRUE);
  for (ss in 1:nSets)
  {
    rownames(eigengenes[[set]][[ss]]$data) = rownames(multiExprEP[[ss]]$data);
  }
}

## Analysis and display of module preservation results

#  Isolate the observed statistics and their Z scores:
# Load the module preservation statistics
ref = 1 # Select the lowland data as reference
test = 2 # Select the highland data as test
statsObsEP = cbind(mpEP$quality$observed[[ref]][[test]][, -1], mpEP$preservation$observed[[ref]][[test]][, -1])
statsZEP = cbind(mpEP$quality$Z[[ref]][[test]][, -1], mpEP$preservation$Z[[ref]][[test]][, -1])
Z.PreservationStatsEP = mpEP$preservation$Z[[ref]][[test]]
Obs.PreservationStatsEP = mpEP$preservation$observed[[ref]][[test]]
# look at the main output: the preservation Zsummary score.
print(signif(statsZEP[, "Zsummary.pres", drop = FALSE],2))
# Compare preservation to quality:
print(signif(statsZEP[, c("Zsummary.pres", "Zsummary.qual")], 2))

# look at the main output: the preservation Zsummary score.
print(signif(statsObsEP[, "Zsummary.pres", drop = FALSE],2))
# Compare preservation to quality:
print(signif(statsObsEP[, c("Zsummary.pres", "Zsummary.qual")], 2))
statsObsEP
statsZEP


## Visualize Data ####
##  Caption: Preservation of highland placenta modules in lowland placenta. 
#The right panel shows the composite statistic medianRank (Eq.9.20) versus module size. 
#he higher the medianRank the less preserved is the module relative to other modules. 
#Since medianRank is based on the observed preservation statistics (as opposed to Z statistics or p-values) 
#we find that it is much less dependent on module size. 
#The upper right panel shows the composite statistic Zsummary (Eq.9.1). 
#If Zsummary> 10 there is strong evidence that the module is preserved (Langfelder et al 2011). 
#If Zsummary<2, there is no evidence that the module preserved. 
#Note that Zsummary shows a strong dependence on module size. 
##
modColorsEP = rownames(Obs.PreservationStatsEP)
moduleSizeEP = Obs.PreservationStatsEP$moduleSize
# we will omit the grey module (background genes)
# and the gold module (random sample of genes)
selectModules = !(modColorsEP %in% c("grey", "gold"))
# Text labels for points
point.labelEP = modColorsEP[selectModules]

#Composite preservation statistics
medianRankEP = Obs.PreservationStatsEP$medianRank.pres
ZsummaryEP = Z.PreservationStatsEP$Zsummary.pres

stats = as.data.frame((signif(statsZEP[, c("Zsummary.pres", "Zsummary.qual")], 2)))

### medianRanks and Zsummary Plot ####
par(mfrow=c(1,2),mar = c(4.5,4.5,2.5,1))
# plot medianRank versus module size: The is useful for comparing relative 
# preservation among multiple modules: a module with lower median rank tends to 
# exhibit stronger observed preservation statistics than a module with a higher 
# median rank.
plot(moduleSizeEP[selectModules],medianRankEP[selectModules],col=1,
     bg=modColorsEP[selectModules],pch = 21,main="medianRank Preservation",
     cex = 1.5, ylab ="medianRank",xlab="Module size", log="x")
labelPoints(moduleSizeEP[selectModules],medianRankEP[selectModules],point.labelEP,cex=1,offs=0.03)

# plot Zsummary versus module size: If Zsummary is greater than 10 it means the module is preserved
# if 2<Zsummary<10 there is weak to moderate preservation, if less than 2
# there is no evidence
plot(moduleSizeEP[selectModules],ZsummaryEP[selectModules], col = 1,
     bg=modColorsEP[selectModules],pch = 21,main="EP Zsummary Preservation",
     cex=1.5,ylab ="Zsummary", xlab = "Module size", log = "x", ylim = c(0,100))
labelPoints(moduleSizeEP[selectModules],ZsummaryEP[selectModules],point.labelEP,cex=1,offs=0.01)
# Add threshold lines for Zsummary
abline(h=0); abline(h=2, col = "blue", lty = 2); abline(h=10, col = "red", lty = 2)


#### Add Zsummary scores to BWtables (i.e., preservation stats for highlanders)
BWpreservedStatsEP = print(signif(statsZEP[, "Zsummary.pres", drop = FALSE],2))



