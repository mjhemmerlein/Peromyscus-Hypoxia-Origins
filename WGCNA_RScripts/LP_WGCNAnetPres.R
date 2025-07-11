
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
     main = paste("Scale indLPendence"));
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
                           saveTOMFileBase = "LP_WGCNA_Output/LPBWExprTOM",
                           verbose = 3)

load("LP_WGCNA_Output/LPBWExprTOM-block.1.RData"); load("LP_WGCNA_Output/LPBWExprTOM-block.2.RData")
table(NetBWLP$colors)
moduleLabelsBWLP = NetBWLP$colors
moduleColorsBWLP = labels2colors(NetBWLP$colors)

MEsBW_LP = NetBWLP$MEs # These are module eigengenes:first principal component of the expression matrix, 
# i.e., weighted average expression profile of a sample in a module
geneTreeBWLP = NetBWLP$dendrograms[[1]]
table(moduleColorsBWLP)
dim(table(moduleColorsBWLP))

plotDendroAndColors(
  NetBWLP$dendrograms[[1]],
  moduleColorsBWLP[NetBWLP$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03, 
  addGuide = TRUE,
  guideHang = 0.05)

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(ExprData_BWLP, moduleColorsBWLP)$eigengenes
MEsBW_LP = orderMEs(MEs0) #the rownames of this dataset are equal to Expr
# write.csv(MEsBW_LP, "LP_WGCNA_Ouput/LP_BW_moduleEigengenes.csv")
# save(NetBWLP, MEsBW_LP, moduleLabelsBWLP, moduleColorsBWLP, geneTreeBWLP, file = "LP_WGCNA_Ouput/BWLP_network.RData")

load("LP_WGCNA_Ouput/BWLP_network.RData")


## BW Trait Module Associations and Correlations ####
ModuleME_Info_BWLP = read.csv("LP_WGCNA_Ouput/LP_BW_moduleEigengenes.csv")
rownames(ModuleME_Info_BWLP) = ModuleME_Info_BWLP$X
ModuleME_Info_BWLP = ModuleME_Info_BWLP %>% dplyr::select(-X)
ModuleME_Info_BWLP = ModuleME_Info_BWLP %>% dplyr::select(-MEgrey) #remove grey module

# Correlate modules to hypoxia OVERALL w/ mom as random effect
Strain = TraitsBW_LP$Strain
O2 = TraitsBW_LP$O2
FetalMass = TraitsBW_LP$Fetus
Sex = TraitsBW_LP$Sex
# summaryBWLP = data.frame()
all_modules = colnames(ModuleME_Info_BWLP)

summaryBWLP02 = data.frame()

for (p in all_modules) {
  
  test_data = ModuleME_Info_BWLP[, p]
  test_model = lmer(test_data ~ O2 + (1|(TraitsBW_LP$Mom)))
  output = anova(test_model)
  output_line <- data.frame(
    Module = p,
    F.O2 = output$`F value`,
    Pvalue.O2 = output$`Pr(>F)`
  )
  summaryBWLP02 = rbind(summaryBWLP02, output_line)
}

summaryBWLP02$Pvalue.O2corr <- p.adjust(summaryBWLP02$Pvalue.O2, method = 'BH')

# write.csv(summaryBWLP02, "LP_WGCNA_Ouput/LP_BW_Module_Trait_ModelSummary.csv")

# Create a dataset containing all gene-specific information
genes=names(ExprData_BWLP)
geneInfoBWLP = data.frame(Gene = genes,
                          moduleColor = moduleColorsBWLP)

write.csv(geneInfoBWLP, "LP_WGCNA_Ouput/LP_BW_Modules.csv")

## Plot eigenes against hypoxia/normoxia in BW animals ####
ModuleME_Info_BWLP$hypoxia = TraitsBW_LP$O2

#brown
ggplot(data = ModuleME_Info_BWLP, aes(x = ModuleME_Info_BWLP$hypoxia, y = ModuleME_Info_BWLP$MEbrown)) + 
  geom_boxplot(aes(fill = ModuleME_Info_BWLP$hypoxia), outlier.shape = NA) + theme_classic() + 
  geom_jitter(width = 0.25) + 
  scale_fill_manual(values = c("red", "blue")) + 
  guides(fill = guide_legend(title = "Oxygen")) + 
  xlab("Oxygen") + ylab("Brown Module Eigengene")



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
     main = paste("Scale indLPendence"));
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
NetMELP = blockwiseModules(ExprData_MELP, 
                           power = 14, 
                           maxBlockSize = 14000,
                           TOMType = "signed", 
                           networkType = "signed hybrid",
                           minModuleSize = 30,
                           reassignThreshold = 0, 
                           mergeCutHeight = 0.25,
                           numericLabels = TRUE, 
                           pamRespectsDendro = FALSE,
                           saveTOMs = TRUE,
                           saveTOMFileBase = "LP_WGCNA_Ouput/LPMEExprTOM",
                           verbose = 3)

load("LP_WGCNA_Ouput/LPMEExprTOM-block.1.RData"); load("LP_WGCNA_Ouput/LPMEExprTOM-block.2.RData")
table(NetMELP$colors)
moduleLabelsMELP = NetMELP$colors
moduleColorsMELP = labels2colors(NetMELP$colors)
MEsMELP = NetMELP$MEs # These are module eigengenes:first principal component of the expression matrix, 
# i.e., weighted average expression profile of a sample in a module

geneTreeMELP = NetMELP$dendrograms[[1]]
table(moduleColorsMELP)
dim(table(moduleColorsMELP))

M1plotDendroAndColors(
  NetMELP$dendrograms[[1]],
  moduleColorsMELP[NetMELP$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03, 
  addGuide = TRUE,
  guideHang = 0.05)

# Recalculate MEs with color labels
MEs0ME = moduleEigengenes(ExprData_MELP, moduleColorsMELP)$eigengenes
MEsMELP = orderMEs(MEs0ME) #the rownames of this dataset are equal to Expr

#Add metadata to this dataframe
MEsMELP$hypoxia = TraitsME_LP$O2
MEsMELP$fetalMass = TraitsME_LP$Fetus

write.csv(MEsMELP, "LP_WGCNA_Ouput/LP_ME_moduleEigengenes.csv")

save(NetMELP, MEsMELP, moduleLabelsMELP, moduleColorsMELP, geneTreeMELP, file = "LP_WGCNA_Ouput/MELP_network.RData")
load("LP_WGCNA_Ouput/MELP_network.RData")


## ME Trait Module Associations and Correlations ####
ModuleME_Info_MELP = read.csv("LP_WGCNA_Ouput/LP_ME_moduleEigengenes.csv")
rownames(ModuleME_Info_MELP) = ModuleME_Info_MELP$X
ModuleME_Info_MELP = ModuleME_Info_MELP %>% dplyr::select(-X)
ModuleME_Info_MELP = ModuleME_Info_MELP %>% dplyr::select(-MEgrey) #remove grey module
ModuleME_Info_MELP = ModuleME_Info_MELP %>% select(-hypoxia)

# Correlate modules to hypoxia OVERALL w/ mom as random effect
Strain = TraitsME_LP$Strain
O2 = TraitsME_LP$O2
FetalMass = TraitsME_LP$Fetus
Sex = TraitsME_LP$Sex
# summaryMELP = data.frame()
all_modules = colnames(ModuleME_Info_MELP)

summaryMELP02 = data.frame()

for (p in all_modules) {
  
  test_data = ModuleME_Info_MELP[, p]
  test_model = lmer(test_data ~ O2 + (1|(TraitsME_LP$Mom)))
  output = anova(test_model)
  output_line <- data.frame(
    Module = p,
    F.O2 = output$`F value`,
    Pvalue.O2 = output$`Pr(>F)`
  )
  summaryMELP02 = rbind(summaryMELP02, output_line)
}

summaryMELP02$Pvalue.O2corr <- p.adjust(summaryMELP02$Pvalue.O2, method = 'BH')

# write.csv(summaryMELP02, "LP_WGCNA_Ouput/LP_ME_Module_Trait_ModelSummary.csv")

# Create a dataset containing all gene-specific information
genes=names(ExprData_MELP)
geneInfoMELP = data.frame(Gene = genes,
                          moduleColor = moduleColorsMELP)

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
  mpLP = WGCNA::modulLPreservation(multiExprLP, colorList,
                                   referenceNetworks = 1,
                                   loadPermutedStatistics = FALSE,
                                   nPermutations = 200,
                                   verbose = 3, 
                                   maxModuleSize = 5301) #set max size to max in lowland network
} )

save(mpLP, file = "LP_WGCNA_Ouput/LPHighLander-Lowlander_modPreservation.RData")
load(file = "LP_WGCNA_Ouput/LPHighLander-Lowlander_modPreservation.RData")

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


# GO Analysis for ME modules ####

musGeneinfo = read_xlsx("RNA_Seq_Output/LP_ISO_Ortho_Summary.xlsx")

## Match and add gene names from Mus from Kates musGeneInfo
## Filter genes from module colors that are in musGeneInfo
library(gprofiler2)
background = as.vector(musGeneinfo$Pman_GeneID)

# black
blackME = geneInfoBWLP %>% filter(moduleColor == "black")
as.data.frame(blackME)
blackME = blackME[, 1:2]

rownames(blackME) = blackME$Gene
intersection2 = musGeneinfo$Pman_GeneID  %in% rownames(blackME)
blackME0 = musGeneinfo[which(intersection2), ] 
blackMEGO = as.vector(blackME0$Pman_GeneID)

blackGO = gost(blackMEGO,
               organism = "mmusculus",
               user_threshold = 0.05,
               custom_bg = background,
               ordered_query = TRUE,
               correction_method = "fdr")

blackGOResults = data.frame(
  Cluster = blackGO$result$query,
  Term.ID = blackGO$result$term_id,
  Term.Name = blackGO$result$term_name,
  geneid = blackGO$result$intersection,
  P.value = blackGO$result$p_value,
  Source = blackGO$result$source,
  Term.Size = blackGO$result$term_size,
  Precision = blackGO$result$precision,
  intersection_size = blackGO$result$intersection_size,
  query_size = blackGO$result$query_size,
  Gene_ratio = as.numeric((blackGO$result$intersection_size/blackGO$result$query_size)))

# write.csv(blackGOResults, 'LP_WGCNA_Ouput/GO_Results/blackGOResults.csv')

# yellow
yellowME = geneInfoBWLP %>% filter(moduleColor == "yellow")
as.data.frame(yellowME)
yellowME = yellowME[, 1:2]

rownames(yellowME) = yellowME$Gene
intersection2 = musGeneinfo$Pman_GeneID  %in% rownames(yellowME)
yellowME0 = musGeneinfo[which(intersection2), ] 
yellowMEGO = as.vector(yellowME0$Pman_GeneID)

yellowGO = gost(yellowMEGO,
                organism = "mmusculus",
                user_threshold = 0.05,
                custom_bg = background,
                ordered_query = TRUE,
                correction_method = "fdr")

yellowGOResults = data.frame(
  Cluster = yellowGO$result$query,
  Term.ID = yellowGO$result$term_id,
  Term.Name = yellowGO$result$term_name,
  geneid = yellowGO$result$intersection,
  P.value = yellowGO$result$p_value,
  Source = yellowGO$result$source,
  Term.Size = yellowGO$result$term_size,
  Precision = yellowGO$result$precision,
  intersection_size = yellowGO$result$intersection_size,
  query_size = yellowGO$result$query_size,
  Gene_ratio = as.numeric((yellowGO$result$intersection_size/yellowGO$result$query_size)))

# write.csv(yellowGOResults, 'LP_WGCNA_Ouput/GO_Results/yellowGOResults.csv')


# magenta
magentaME = geneInfoBWLP %>% filter(moduleColor == "magenta")
as.data.frame(magentaME)
magentaME = magentaME[, 1:2]

rownames(magentaME) = magentaME$Gene
intersection2 = musGeneinfo$Pman_GeneID  %in% rownames(magentaME)
magentaME0 = musGeneinfo[which(intersection2), ] 
magentaMEGO = as.vector(magentaME0$Pman_GeneID)

magentaGO = gost(magentaMEGO,
                 organism = "mmusculus",
                 user_threshold = 0.05,
                 custom_bg = background,
                 ordered_query = TRUE,
                 correction_method = "fdr")

magentaGOResults = data.frame(
  Cluster = magentaGO$result$query,
  Term.ID = magentaGO$result$term_id,
  Term.Name = magentaGO$result$term_name,
  geneid = magentaGO$result$intersection,
  P.value = magentaGO$result$p_value,
  Source = magentaGO$result$source,
  Term.Size = magentaGO$result$term_size,
  Precision = magentaGO$result$precision,
  intersection_size = magentaGO$result$intersection_size,
  query_size = magentaGO$result$query_size,
  Gene_ratio = as.numeric((magentaGO$result$intersection_size/magentaGO$result$query_size)))

# write.csv(magentaGOResults, 'LP_WGCNA_Ouput/GO_Results/magentaGOResults.csv')

# purple
purpleME = geneInfoBWLP %>% filter(moduleColor == "purple")
as.data.frame(purpleME)
purpleME = purpleME[, 1:2]

rownames(purpleME) = purpleME$Gene
intersection2 = musGeneinfo$Pman_GeneID  %in% rownames(purpleME)
purpleME0 = musGeneinfo[which(intersection2), ] 
purpleMEGO = as.vector(purpleME0$Pman_GeneID)

purpleGO = gost(purpleMEGO,
                organism = "mmusculus",
                user_threshold = 0.05,
                custom_bg = background,
                ordered_query = TRUE,
                correction_method = "fdr")

purpleGOResults = data.frame(
  Cluster = purpleGO$result$query,
  Term.ID = purpleGO$result$term_id,
  Term.Name = purpleGO$result$term_name,
  geneid = purpleGO$result$intersection,
  P.value = purpleGO$result$p_value,
  Source = purpleGO$result$source,
  Term.Size = purpleGO$result$term_size,
  Precision = purpleGO$result$precision,
  intersection_size = purpleGO$result$intersection_size,
  query_size = purpleGO$result$query_size,
  Gene_ratio = as.numeric((purpleGO$result$intersection_size/purpleGO$result$query_size)))

# write.csv(purpleGOResults, 'LP_WGCNA_Ouput/GO_Results/purpleGOResults.csv')

# red
redME = geneInfoBWLP %>% filter(moduleColor == "red")
as.data.frame(redME)
redME = redME[, 1:2]

rownames(redME) = redME$Gene
intersection2 = musGeneinfo$Pman_GeneID  %in% rownames(redME)
redME0 = musGeneinfo[which(intersection2), ] 
redMEGO = as.vector(redME0$Pman_GeneID)

redGO = gost(redMEGO,
             organism = "mmusculus",
             user_threshold = 0.05,
             custom_bg = background,
             ordered_query = TRUE,
             correction_method = "fdr")

redGOResults = data.frame(
  Cluster = redGO$result$query,
  Term.ID = redGO$result$term_id,
  Term.Name = redGO$result$term_name,
  geneid = redGO$result$intersection,
  P.value = redGO$result$p_value,
  Source = redGO$result$source,
  Term.Size = redGO$result$term_size,
  Precision = redGO$result$precision,
  intersection_size = redGO$result$intersection_size,
  query_size = redGO$result$query_size,
  Gene_ratio = as.numeric((redGO$result$intersection_size/redGO$result$query_size)))

# write.csv(redGOResults, 'LP_WGCNA_Ouput/GO_Results/redGOResults.csv')

# pink
pinkME = geneInfoBWLP %>% filter(moduleColor == "pink")
as.data.frame(pinkME)
pinkME = pinkME[, 1:2]

rownames(pinkME) = pinkME$Gene
intersection2 = musGeneinfo$Pman_GeneID  %in% rownames(pinkME)
pinkME0 = musGeneinfo[which(intersection2), ] 
pinkMEGO = as.vector(pinkME0$Pman_GeneID)

pinkGO = gost(pinkMEGO,
              organism = "mmusculus",
              user_threshold = 0.05,
              custom_bg = background,
              ordered_query = TRUE,
              correction_method = "fdr")

pinkGOResults = data.frame(
  Cluster = pinkGO$result$query,
  Term.ID = pinkGO$result$term_id,
  Term.Name = pinkGO$result$term_name,
  geneid = pinkGO$result$intersection,
  P.value = pinkGO$result$p_value,
  Source = pinkGO$result$source,
  Term.Size = pinkGO$result$term_size,
  Precision = pinkGO$result$precision,
  intersection_size = pinkGO$result$intersection_size,
  query_size = pinkGO$result$query_size,
  Gene_ratio = as.numeric((pinkGO$result$intersection_size/pinkGO$result$query_size)))

# write.csv(pinkGOResults, 'LP_WGCNA_Ouput/GO_Results/pinkGOResults.csv')

# turquoise
turquoiseME = geneInfoBWLP %>% filter(moduleColor == "turquoise")
as.data.frame(turquoiseME)
turquoiseME = turquoiseME[, 1:2]

rownames(turquoiseME) = turquoiseME$Gene
intersection2 = musGeneinfo$Pman_GeneID  %in% rownames(turquoiseME)
turquoiseME0 = musGeneinfo[which(intersection2), ] 
turquoiseMEGO = as.vector(turquoiseME0$Pman_GeneID)

turquoiseGO = gost(turquoiseMEGO,
                   organism = "mmusculus",
                   user_threshold = 0.05,
                   custom_bg = background,
                   ordered_query = TRUE,
                   correction_method = "fdr")

turquoiseGOResults = data.frame(
  Cluster = turquoiseGO$result$query,
  Term.ID = turquoiseGO$result$term_id,
  Term.Name = turquoiseGO$result$term_name,
  geneid = turquoiseGO$result$intersection,
  P.value = turquoiseGO$result$p_value,
  Source = turquoiseGO$result$source,
  Term.Size = turquoiseGO$result$term_size,
  Precision = turquoiseGO$result$precision,
  intersection_size = turquoiseGO$result$intersection_size,
  query_size = turquoiseGO$result$query_size,
  Gene_ratio = as.numeric((turquoiseGO$result$intersection_size/turquoiseGO$result$query_size)))

# write.csv(turquoiseGOResults, 'LP_WGCNA_Ouput/GO_Results/turquoiseGOResults.csv')

# blue
blueME = geneInfoBWLP %>% filter(moduleColor == "blue")
as.data.frame(blueME)
blueME = blueME[, 1:2]

rownames(blueME) = blueME$Gene
intersection2 = musGeneinfo$Pman_GeneID  %in% rownames(blueME)
blueME0 = musGeneinfo[which(intersection2), ] 
blueMEGO = as.vector(blueME0$Pman_GeneID)

blueGO = gost(blueMEGO,
              organism = "mmusculus",
              user_threshold = 0.05,
              custom_bg = background,
              ordered_query = TRUE,
              correction_method = "fdr")

blueGOResults = data.frame(
  Cluster = blueGO$result$query,
  Term.ID = blueGO$result$term_id,
  Term.Name = blueGO$result$term_name,
  geneid = blueGO$result$intersection,
  P.value = blueGO$result$p_value,
  Source = blueGO$result$source,
  Term.Size = blueGO$result$term_size,
  Precision = blueGO$result$precision,
  intersection_size = blueGO$result$intersection_size,
  query_size = blueGO$result$query_size,
  Gene_ratio = as.numeric((blueGO$result$intersection_size/blueGO$result$query_size)))

# write.csv(blueGOResults, 'LP_WGCNA_Ouput/GO_Results/blueGOResults.csv')

# greenyellow
greenyellowME = geneInfoBWLP %>% filter(moduleColor == "greenyellow")
as.data.frame(greenyellowME)
greenyellowME = greenyellowME[, 1:2]

rownames(greenyellowME) = greenyellowME$Gene
intersection2 = musGeneinfo$Pman_GeneID  %in% rownames(greenyellowME)
greenyellowME0 = musGeneinfo[which(intersection2), ] 
greenyellowMEGO = as.vector(greenyellowME0$Pman_GeneID)

greenyellowGO = gost(greenyellowMEGO,
                     organism = "mmusculus",
                     user_threshold = 0.05,
                     custom_bg = background,
                     ordered_query = TRUE,
                     correction_method = "fdr")

greenyellowGOResults = data.frame(
  Cluster = greenyellowGO$result$query,
  Term.ID = greenyellowGO$result$term_id,
  Term.Name = greenyellowGO$result$term_name,
  geneid = greenyellowGO$result$intersection,
  P.value = greenyellowGO$result$p_value,
  Source = greenyellowGO$result$source,
  Term.Size = greenyellowGO$result$term_size,
  Precision = greenyellowGO$result$precision,
  intersection_size = greenyellowGO$result$intersection_size,
  query_size = greenyellowGO$result$query_size,
  Gene_ratio = as.numeric((greenyellowGO$result$intersection_size/greenyellowGO$result$query_size)))

# write.csv(greenyellowGOResults, 'LP_WGCNA_Ouput/GO_Results/greenyellowGOResults.csv')

# green
greenME = geneInfoBWLP %>% filter(moduleColor == "green")
as.data.frame(greenME)
greenME = greenME[, 1:2]

rownames(greenME) = greenME$Gene
intersection2 = musGeneinfo$Pman_GeneID  %in% rownames(greenME)
greenME0 = musGeneinfo[which(intersection2), ] 
greenMEGO = as.vector(greenME0$Pman_GeneID)

greenGO = gost(greenMEGO,
               organism = "mmusculus",
               user_threshold = 0.05,
               custom_bg = background,
               ordered_query = TRUE,
               correction_method = "fdr")

greenGOResults = data.frame(
  Cluster = greenGO$result$query,
  Term.ID = greenGO$result$term_id,
  Term.Name = greenGO$result$term_name,
  geneid = greenGO$result$intersection,
  P.value = greenGO$result$p_value,
  Source = greenGO$result$source,
  Term.Size = greenGO$result$term_size,
  Precision = greenGO$result$precision,
  intersection_size = greenGO$result$intersection_size,
  query_size = greenGO$result$query_size,
  Gene_ratio = as.numeric((greenGO$result$intersection_size/greenGO$result$query_size)))

# write.csv(greenGOResults, 'LP_WGCNA_Ouput/GO_Results/greenGOResults.csv')
