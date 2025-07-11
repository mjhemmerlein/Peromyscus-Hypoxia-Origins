setwd("~/Library/CloudStorage/OneDrive-Colostate/CSU/PlasticityBioinformatics_KWCD/Dream analysis/Wilsterman_resample/WGCNANetworkPreservation/")
# Load Packages ####
#BiocManager::install("WGCNA")
library(WGCNA)
library(dplyr)
library(edgeR)
library(lme4)
library(lmerTest)
###############################################################
  
# Separate WGCNA for BW and ME #####

  #Because Langfelder says so
options(stringsAsFactors = FALSE)

#Read in sample info metadata
Sample_Info_combined = read.csv("Metadata_LP_filtered.csv")

Sample_Info_combined = Sample_Info_combined %>% select(-Placenta)
Traits = Sample_Info_combined
Traits$Sample_ID_LZ = as.factor(Traits$Sample_ID_LZ)
Traits$LZ_Lane = as.factor(Traits$LZ_Lane)
Traits$Sample_ID_JZ= as.factor(Traits$Sample_ID_JZ)
Traits$JZ_Lane = as.factor(Traits$JZ_Lane)
Traits$TubeID = as.factor(Traits$TubeID)
Traits$Mom = as.factor(Traits$Mom)
Traits$Strain = as.factor(Traits$Strain)
Traits$O2 = as.factor(Traits$O2)
Traits$Fetus = as.numeric(Traits$Fetus)
Traits$Group = as.factor(Traits$Group)
Traits$Sex = as.factor(Traits$Sex)
Traits = Traits %>% select(c(-LS, -C, -Stage, -Viable))

TraitsBW_LZ = Traits %>% filter(Strain == "BW")
TraitsBW_LZ = TraitsBW_LZ[1:38, ]
rownames(TraitsBW_LZ) = TraitsBW_LZ$Sample_ID_LZ

TraitsME_LZ = Traits %>% filter(Strain == "ME")
TraitsME_LZ = TraitsME_LZ[1:40, ]
rownames(TraitsME_LZ) = TraitsME_LZ$Sample_ID_LZ

#normalize count data
## input read counts for LZ
LZcounts = read.table("../LZ_allreads.txt", header = TRUE)
rownames(LZcounts) = LZcounts$Geneid
LZcounts = LZcounts[, 7:85]
LZcounts = subset(LZcounts, select = -c(LZ089))

ExprData0 = DGEList(LZcounts)
ExprData0 = calcNormFactors(ExprData0)

keep = rowSums(cpm(ExprData0) > 0.5) >= 59
ExprData = ExprData0[keep,]
dim(ExprData)

### WGCNA Analysis
ExprData_LZ = cpm(ExprData, log=TRUE, prior.count=2, normalized.lib.sizes=TRUE) #cpm normalized and log transformed expression dataset

#Filter ME and BW
ExprData_MELZ = ExprData_LZ[, 39:78]
ExprData_BWLZ = ExprData_LZ[, 1:38]

#transpose expression data for further analysis
#ME
ExprData_MELZ = as.data.frame(t(ExprData_MELZ))
ExprData_MELZ$treatment = TraitsME_LZ$Strain
ExprData_MELZ = ExprData_MELZ %>% select(-treatment)

#BW
ExprData_BWLZ = as.data.frame(t(ExprData_BWLZ))
ExprData_BWLZ$treatment = TraitsBW_LZ$Strain
ExprData_BWLZ = ExprData_BWLZ %>% select(-treatment)

match(rownames(TraitsME_LZ), rownames(ExprData_MELZ))# should be in order if datasets are matched
match(rownames(TraitsBW_LZ), rownames(ExprData_BWLZ))

## WGCNA

### Lowlanders, i.e., BW ####
#cluster the samples to check for outliers
sampleTreeBW = hclust(dist(ExprData_BWLZ), method = "average")
# Plot the sample trees
plot(sampleTreeBW, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2) 
## Set power threshold
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sftBW = pickSoftThreshold(ExprData_BWLZ, powerVector = powers, verbose = 5, 
                          networkType = "signed hybrid") 
# Plot the results:
#pdf('wgcna/rv_beta_plot.pdf', h=4, w=7)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sftBW$fitIndices[,1], -sign(sftBW$fitIndices[,3])*sftBW$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sftBW$fitIndices[,1], -sign(sftBW$fitIndices[,3])*sftBW$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.8,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sftBW$fitIndices[,1], sftBW$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sftBW$fitIndices[,1], sftBW$fitIndices[,5], labels=powers, cex=cex1,col="red")

# Call the network topology analysis function
NetBWLZ = blockwiseModules(ExprData_BWLZ, 
                           power = 5, 
                           maxBlockSize = 14000,
                           TOMType = "signed", networkType = "signed hybrid",
                           minModuleSize = 30,
                           reassignThreshold = 0,
                           mergeCutHeight = 0.25,
                           numericLabels = TRUE, pamRespectsDendro = FALSE,
                           saveTOMs = TRUE,
                           saveTOMFileBase = "LZBWExprTOM",
                           verbose = 3)

load("LZBWExprTOM-block.1.RData"); load("LZBWExprTOM-block.2.RData")
table(NetBWLZ$colors)
moduleLabelsBWLZ = NetBWLZ$colors
moduleColorsBWLZ = labels2colors(NetBWLZ$colors)
MEsBW_LZ = NetBWLZ$MEs # These are module eigengenes:first principal component of the expression matrix, 
# i.e., weighted average expression profile of a sample in a module
geneTreeBWLZ = NetBWLZ$dendrograms[[1]]
table(moduleColorsBWLZ)
dim(table(moduleColorsBWLZ))

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(ExprData_BWLZ, moduleColorsBWLZ)$eigengenes
MEsBW_LZ = orderMEs(MEs0) #the rownames of this dataset are equal to Expr
write.csv(MEsBW_LZ, "Labyrinth Zone/BWLZ_moduleEigengenes.csv")

save(NetBWLZ, MEsBW_LZ, moduleLabelsBWLZ, moduleColorsBWLZ, geneTreeBWLZ, file = "Labyrinth Zone/BWLZ_network.RData")
load("BWLZ_network.RData")

## BW Trait Module Associations and Correlations ####
ModuleME_Info_BWLZ = read.csv("Labyrinth Zone/BWLZ_moduleEigengenes.csv")
rownames(ModuleME_Info_BWLZ) = ModuleME_Info_BWLZ$X
ModuleME_Info_BWLZ = ModuleME_Info_BWLZ %>% dplyr::select(-X)
ModuleME_Info_BWLZ = ModuleME_Info_BWLZ %>% dplyr::select(-MEgrey) #remove grey module

# Correlate modules to fetal weight OVERALL w/ mom as random effect
Strain = TraitsBW_LZ$Strain
O2 = TraitsBW_LZ$O2
FetalMass = TraitsBW_LZ$Fetus
Sex = TraitsBW_LZ$Sex
summaryBWLZ = data.frame()
all_modules = colnames(ModuleME_Info_BWLZ)

for (p in all_modules) {
  
  test_data = ModuleME_Info_BWLZ[, p]
  test_model = lmer(test_data ~ FetalMass + (1|(TraitsBW_LZ$Mom)))
  output = anova(test_model)
  output_line <- data.frame(
    Module = p,
    F.FetalMass = output$`F value`,
    Pvalue.FetalMass = output$`Pr(>F)`
  )
  summaryBWLZ = rbind(summaryBWLZ, output_line)
}

# Correlate modules to hypoxia OVERALL w/ mom as random effect
summaryBWLZ02 = data.frame()

for (p in all_modules) {
  
  test_data = ModuleME_Info_BWLZ[, p]
  test_model = lmer(test_data ~ O2 + (1|(TraitsBW_LZ$Mom)))
  output = anova(test_model)
  output_line <- data.frame(
    Module = p,
    F.O2 = output$`F value`,
    Pvalue.O2 = output$`Pr(>F)`
  )
  summaryBWLZ02 = rbind(summaryBWLZ02, output_line)
}

summaryBWLZ = cbind(summaryBWLZ, summaryBWLZ02[, 2:3])


# Correlate modules to Sex OVERALL w/ mom as random effect
summaryBWLZSex = data.frame()

for (p in all_modules) {
  
  test_data = ModuleME_Info_BWLZ[, p]
  test_model = lmer(test_data ~ Sex + (1|(TraitsBW_LZ$Mom)))
  output = anova(test_model)
  output_line <- data.frame(
    Module = p,
    F.Sex = output$`F value`,
    Pvalue.Sex = output$`Pr(>F)`
  )
  summaryBWLZSex = rbind(summaryBWLZSex, output_line)
}

summaryBWLZ = cbind(summaryBWLZ, summaryBWLZSex[, 2:3])

# Correlate modules to Mass*O2 OVERALL w/ mom as random effect
summaryBWLZmassO2 = data.frame()

for (p in all_modules) {
  
  test_data = ModuleME_Info_BWLZ[, p]
  test_model = lmer(test_data ~ FetalMass*O2 + (1|(TraitsBW_LZ$Mom)))
  output = anova(test_model)
  output_line <- data.frame(
    Module = p,
    F.massO2 = output[3,]$`F value`,
    Pvalue.massO2 = output[3,]$`Pr(>F)`
  )
  summaryBWLZmassO2 = rbind(summaryBWLZmassO2, output_line)
}

summaryBWLZ = cbind(summaryBWLZ, summaryBWLZmassO2[, 2:3])

# Correlate modules to Mass*Sex OVERALL w/ mom as random effect
summaryBWLZmassSex = data.frame()

for (p in all_modules) {
  
  test_data = ModuleME_Info_BWLZ[, p]
  test_model = lmer(test_data ~ FetalMass*Sex + (1|(TraitsBW_LZ$Mom)))
  output = anova(test_model)
  output_line <- data.frame(
    Module = p,
    F.massSex = output[3,]$`F value`,
    Pvalue.massSex = output[3,]$`Pr(>F)`
  )
  summaryBWLZmassSex = rbind(summaryBWLZmassSex, output_line)
}

summaryBWLZ = cbind(summaryBWLZ, summaryBWLZmassSex[, 2:3])

# Correlate modules to O2*Sex OVERALL w/ mom as random effect
summaryBWLZO2Sex = data.frame()

for (p in all_modules) {
  
  test_data = ModuleME_Info_BWLZ[, p]
  test_model = lmer(test_data ~ O2*Sex + (1|(TraitsBW_LZ$Mom)))
  output = anova(test_model)
  output_line <- data.frame(
    Module = p,
    F.O2Sex = output[3,]$`F value`,
    Pvalue.O2Sex = output[3,]$`Pr(>F)`
  )
  summaryBWLZO2Sex = rbind(summaryBWLZO2Sex, output_line)
}

summaryBWLZ = cbind(summaryBWLZ, summaryBWLZO2Sex[, 2:3])

# Correlate modules to O2*Sex OVERALL w/ mom as random effect
summaryBWLZMassO2Sex = data.frame()

for (p in all_modules) {
  
  test_data = ModuleME_Info_BWLZ[, p]
  test_model = lmer(test_data ~ FetalMass*O2*Sex + (1|(TraitsBW_LZ$Mom)))
  output = anova(test_model)
  output_line <- data.frame(
    Module = p,
    F.O2Sex = output[3,]$`F value`,
    Pvalue.MassO2Sex = output[3,]$`Pr(>F)`
  )
  summaryBWLZMassO2Sex = rbind(summaryBWLZMassO2Sex, output_line)
}

summaryBWLZ = cbind(summaryBWLZ, summaryBWLZMassO2Sex[, 2:3])
moduleColorsBWLZTable = as.data.frame(table(moduleColorsBWLZ))
moduleColorsBWLZTable = moduleColorsBWLZTable %>% filter(moduleColorsBWLZTable$moduleColorsBWLZ != "grey")
summaryBWLZ$ModuleSize = moduleColorsBWLZTable$Freq
summaryBWLZ = summaryBWLZ %>% relocate(ModuleSize, .after = Module)

summaryBWLZ$Pvalue.FetalMasscorr <- p.adjust(summaryBWLZ$Pvalue.FetalMass, method = 'BH')
summaryBWLZ$Pvalue.O2corr <- p.adjust(summaryBWLZ$Pvalue.O2, method = 'BH')
summaryBWLZ$Pvalue.Sexcorr <- p.adjust(summaryBWLZ$Pvalue.Sex, method = 'BH')
summaryBWLZ$Pvalue.massO2corr <- p.adjust(summaryBWLZ$Pvalue.massO2, method = 'BH')
summaryBWLZ$Pvalue.massSexcorr <- p.adjust(summaryBWLZ$Pvalue.massSex, method = 'BH')
summaryBWLZ$Pvalue.O2Sexcorr <- p.adjust(summaryBWLZ$Pvalue.O2Sex, method = 'BH')
summaryBWLZ$Pvalue.MassO2Sexcorr <- p.adjust(summaryBWLZ$Pvalue.MassO2Sex, method = 'BH')

write.csv(summaryBWLZ, "Labyrinth Zone/BWLZModule_Trait_ModelSummary.csv")

## Plot eigenes against hypoxia/normoxia in BW animals ####
ModuleME_Info_BWLZ$hypoxia = TraitsBW_LZ$O2
ModuleME_Info_BWLZ$fetalMass = TraitsBW_JZ$Fetus
ModuleME_Info_BWLZ$Sex = TraitsBW_JZ$Sex

#lightyellow
ggplot(data = ModuleME_Info_BWLZ, aes(x = ModuleME_Info_BWLZ$hypoxia, y = ModuleME_Info_BWLZ$MElightyellow)) + 
  geom_boxplot(aes(fill = ModuleME_Info_BWLZ$hypoxia), outlier.shape = NA) + theme_classic() + 
  geom_jitter(width = 0.25) + 
  scale_fill_manual(values = c("red", "blue")) + 
  guides(fill = guide_legend(title = "Oxygen")) + 
  xlab("Oxygen") + ylab("Lightyellow Module Eigengene")

ggplot(data = ModuleME_Info_BWLZ, aes(x = ModuleME_Info_BWLZ$fetalMass, y = ModuleME_Info_BWLZ$MElightyellow)) + 
  geom_point(aes(fill = ModuleME_Info_BWLZ$Sex, color = ModuleME_Info_BWLZ$Sex)) + 
  stat_smooth(method = "lm",
              formula = y ~ x, aes(group = Sex, fill = Sex),
              geom = "smooth", color = "black") + theme_classic() + 
  scale_fill_manual(values = c("red", "blue")) + 
  scale_color_manual(values = c("red", "blue")) + 
  guides(color = guide_legend(title = "Sex")) + 
  guides(fill = guide_legend(title = "Sex")) + 
  xlab("Fetal Mass") +
  ylab("Lightyellow Module Eigengene") 
  
ggplot(data = ModuleME_Info_BWLZ, aes(x = ModuleME_Info_BWLZ$hypoxia, y = ModuleME_Info_BWLZ$MElightyellow)) + 
  geom_boxplot(aes(fill = Sex), outlier.shape = NA) + 
  geom_point(aes(color = Sex), position = position_jitterdodge()) + 
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c('black', "black")) + 
  guides(fill = guide_legend(title = "Sex")) + 
  xlab("Oxygen") +
  ylab("Lightyellow Module Eigengene") + theme_classic()

#midnightblueINX
ggplot(data = ModuleME_Info_BWLZ, aes(x = ModuleME_Info_BWLZ$fetalMass, y = ModuleME_Info_BWLZ$MEmidnightblue)) + 
  geom_point(aes(fill = ModuleME_Info_BWLZ$Sex, color = ModuleME_Info_BWLZ$Sex)) + 
  stat_smooth(method = "lm",
              formula = y ~ x, aes(group = Sex, fill = Sex),
              geom = "smooth", color = "black") + theme_classic() + 
  scale_fill_manual(values = c("red", "blue")) + 
  scale_color_manual(values = c("red", "blue")) + 
  guides(color = guide_legend(title = "Sex")) + 
  guides(fill = guide_legend(title = "Sex")) + 
  xlab("Fetal Mass") +
  ylab("Midnightblue Module Eigengene") 

#redINX
ggplot(data = ModuleME_Info_BWLZ, aes(x = ModuleME_Info_BWLZ$fetalMass, y = ModuleME_Info_BWLZ$MEred)) + 
  geom_point(aes(fill = ModuleME_Info_BWLZ$Sex, color = ModuleME_Info_BWLZ$Sex)) + 
  stat_smooth(method = "lm",
              formula = y ~ x, aes(group = Sex, fill = Sex),
              geom = "smooth", color = "black") + theme_classic() + 
  scale_fill_manual(values = c("red", "blue")) + 
  scale_color_manual(values = c("red", "blue")) + 
  guides(color = guide_legend(title = "Sex")) + 
  guides(fill = guide_legend(title = "Sex")) + 
  xlab("Fetal Mass") +
  ylab("Red Module Eigengene") 

#greenyellowO2
ggplot(data = ModuleME_Info_BWLZ, aes(x = ModuleME_Info_BWLZ$hypoxia, y = ModuleME_Info_BWLZ$MEgreenyellow)) + 
  geom_boxplot(aes(fill = ModuleME_Info_BWLZ$hypoxia), outlier.shape = NA) + theme_classic() + 
  geom_jitter(width = 0.25) + 
  scale_fill_manual(values = c("red", "blue")) + 
  guides(fill = guide_legend(title = "Oxygen")) + 
  xlab("Oxygen") + ylab("Greenyellow Module Eigengene")

#darkgreen
ggplot(data = ModuleME_Info_BWLZ, aes(x = ModuleME_Info_BWLZ$hypoxia, y = ModuleME_Info_BWLZ$MEdarkgreen)) + 
  geom_boxplot(aes(fill = ModuleME_Info_BWLZ$hypoxia),outlier.shape = NA) + theme_classic() + 
  geom_jitter(width = 0.25) + 
  scale_fill_manual(values = c("red", "blue")) + 
  guides(fill = guide_legend(title = "Oxygen")) + 
  xlab("Oxygen") + ylab("Darkgreen Module Eigengene")

#green
ggplot(data = ModuleME_Info_BWLZ, aes(x = ModuleME_Info_BWLZ$hypoxia, y = ModuleME_Info_BWLZ$MEgreen)) + 
  geom_boxplot(aes(fill = ModuleME_Info_BWLZ$hypoxia), outlier.shape = NA) + theme_classic() + 
  geom_jitter(width = 0.25) + 
  scale_fill_manual(values = c("red", "blue")) + 
  guides(fill = guide_legend(title = "Oxygen")) + 
  xlab("Oxygen") + ylab("Green Module Eigengene")

# Create a dataset containing all gene-specific information
genes=names(ExprData_BWLZ)
geneInfoBWLZ = data.frame(Gene = genes,
                       moduleColor = moduleColorsBWLZ)

## Highlanders WGCNA network build, i.e., ME ####
#cluster the samples to check for outliers
sampleTreeMELZ = hclust(dist(ExprData_MELZ), method = "average")
# Plot the sample trees
plot(sampleTreeMELZ, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2) 
## Set power threshold
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sftME = pickSoftThreshold(ExprData_MELZ, powerVector = powers, verbose = 5, 
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
NetMELZ = blockwiseModules(ExprData_MELZ, power = 5, maxBlockSize = 14000,
                         TOMType = "signed", networkType = "signed hybrid",
                         minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = "ExprTOM",
                         verbose = 3)

load("ExprTOM-block.1.RData"); load("ExprTOM-block.2.RData")
table(NetMELZ$colors)
moduleLabelsMELZ = NetMELZ$colors
moduleColorsMELZ = labels2colors(NetMELZ$colors)
MEsMELZ = NetMELZ$MEs # These are module eigengenes:first principal component of the expression matrix, 
# i.e., weighted average expression profile of a sample in a module
geneTreeMELZ = NetMELZ$dendrograms[[1]]
table(moduleColorsMELZ)
dim(table(moduleColorsMELZ))

# Recalculate MEs with color labels
MEs0ME = moduleEigengenes(ExprData_MELZ, moduleColorsMELZ)$eigengenes
MEsMELZ = orderMEs(MEs0ME) #the rownames of this dataset are equal to Expr

#Add metadata to this dataframe
MEsMELZ$hypoxia = TraitsME_LZ$O2
MEsMELZ$fetalMass = TraitsME_LZ$Fetus

write.csv(MEsMELZ, "Labyrinth Zone/MELZ_ModuleEigengenes.csv")

save(NetMELZ, MEsMELZ, moduleLabelsMELZ, moduleColorsMELZ, geneTreeMELZ, file = "Labyrinth Zone/MELZ_network.RData")
load("Labyrinth Zone/MELZ_network.RData")
##

## ME Trait Module Associations and Correlations ####
ModuleME_Info_MELZ = read.csv("Labyrinth Zone/MELZ_moduleEigengenes.csv")
rownames(ModuleME_Info_MELZ) = ModuleME_Info_MELZ$X
ModuleME_Info_MELZ = ModuleME_Info_MELZ %>% dplyr::select(-X)
ModuleME_Info_MELZ = ModuleME_Info_MELZ %>% dplyr::select(-MEgrey) #remove grey module
ModuleME_Info_MELZ = ModuleME_Info_MELZ %>% select(-c(hypoxia, fetalMass))

# Correlate modules to fetal weight OVERALL w/ mom as random effect
Strain = TraitsME_LZ$Strain
O2 = TraitsME_LZ$O2
FetalMass = TraitsME_LZ$Fetus
Sex = TraitsME_LZ$Sex
summaryMELZ = data.frame()
all_modules = colnames(ModuleME_Info_MELZ)

for (p in all_modules) {
  
  test_data = ModuleME_Info_MELZ[, p]
  test_model = lmer(test_data ~ FetalMass + (1|(TraitsME_LZ$Mom)))
  output = anova(test_model)
  output_line <- data.frame(
    Module = p,
    F.FetalMass = output$`F value`,
    Pvalue.FetalMass = output$`Pr(>F)`
  )
  summaryMELZ = rbind(summaryMELZ, output_line)
}

# Correlate modules to hypoxia OVERALL w/ mom as random effect
summaryME02 = data.frame()

for (p in all_modules) {
  
  test_data = ModuleME_Info_MELZ[, p]
  test_model = lmer(test_data ~ O2 + (1|(TraitsME_LZ$Mom)))
  output = anova(test_model)
  output_line <- data.frame(
    Module = p,
    F.O2 = output$`F value`,
    Pvalue.O2 = output$`Pr(>F)`
  )
  summaryME02 = rbind(summaryME02, output_line)
}

summaryMELZ = cbind(summaryMELZ, summaryME02[, 2:3])


# Correlate modules to Sex OVERALL w/ mom as random effect
summaryMESex = data.frame()

for (p in all_modules) {
  
  test_data = ModuleME_Info_MELZ[, p]
  test_model = lmer(test_data ~ Sex + (1|(TraitsME_LZ$Mom)))
  output = anova(test_model)
  output_line <- data.frame(
    Module = p,
    F.Sex = output$`F value`,
    Pvalue.Sex = output$`Pr(>F)`
  )
  summaryMESex = rbind(summaryMESex, output_line)
}

summaryMELZ = cbind(summaryMELZ, summaryMESex[, 2:3])

# Correlate modules to FetalMass*O2 inx OVERALL w/ mom as random effect
summaryMEMassO2inx = data.frame()

for (p in all_modules) {
  
  test_data = ModuleME_Info_MELZ[, p]
  test_model = lmer(test_data ~ FetalMass*O2 + (1|(TraitsME_LZ$Mom)))
  output = anova(test_model)
  output_line <- data.frame(
    Module = p,
    F.massO2inx = output[3,]$`F value`,
    Pvalue.massO2inx = output[3,]$`Pr(>F)`
  )
  summaryMEMassO2inx = rbind(summaryMEMassO2inx, output_line)
}

summaryMELZ = cbind(summaryMELZ, summaryMEMassO2inx[, 2:3])

# Correlate modules to Sex*O2 inx OVERALL w/ mom as random effect
summaryMESexO2inx = data.frame()

for (p in all_modules) {
  
  test_data = ModuleME_Info_MELZ[, p]
  test_model = lmer(test_data ~ Sex*O2 + (1|(TraitsME_LZ$Mom)))
  output = anova(test_model)
  output_line <- data.frame(
    Module = p,
    F.SexO2inx = output[3,]$`F value`,
    Pvalue.SexO2inx = output[3,]$`Pr(>F)`
  )
  summaryMESexO2inx = rbind(summaryMESexO2inx, output_line)
}

summaryMELZ = cbind(summaryMELZ, summaryMESexO2inx[, 2:3])

# Correlate modules to FetallMass*Sex inx OVERALL w/ mom as random effect
summaryMEmassSexinx = data.frame()

for (p in all_modules) {
  
  test_data = ModuleME_Info_MELZ[, p]
  test_model = lmer(test_data ~ FetalMass*Sex + (1|(TraitsME_LZ$Mom)))
  output = anova(test_model)
  output_line <- data.frame(
    Module = p,
    F.massSexinx = output[3,]$`F value`,
    Pvalue.massSexinx = output[3,]$`Pr(>F)`
  )
  summaryMEmassSexinx = rbind(summaryMEmassSexinx, output_line)
}

summaryMELZ = cbind(summaryMELZ, summaryMEmassSexinx[, 2:3])
table(moduleColorsMELZ)

moduleColorsMELZTable = as.data.frame(table(moduleColorsMELZ))
moduleColorsMELZTable = moduleColorsMELZTable %>% filter(moduleColorsMELZTable$moduleColorsMELZ != "grey")
summaryMELZ$ModuleSize = moduleColorsMELZTable$Freq
summaryMELZ = summaryMELZ %>% relocate(ModuleSize, .after = Module)

summaryMELZ$Pvalue.FetalMasscorr <- p.adjust(summaryMELZ$Pvalue.FetalMass, method = 'BH')
summaryMELZ$Pvalue.O2corr <- p.adjust(summaryMELZ$Pvalue.O2, method = 'BH')
summaryMELZ$Pvalue.Sexcorr <- p.adjust(summaryMELZ$Pvalue.Sex, method = 'BH')
summaryMELZ$Pvalue.massO2corr <- p.adjust(summaryMELZ$Pvalue.massO2inx, method = 'BH')
summaryMELZ$Pvalue.massSexcorr <- p.adjust(summaryMELZ$Pvalue.massSexinx, method = 'BH')
summaryMELZ$Pvalue.O2Sexcorr <- p.adjust(summaryMELZ$Pvalue.SexO2inx, method = 'BH')

write.csv(summaryMELZ, "Labyrinth Zone/MELZModule_Trait_ModelSummary.csv")

## Plot eigengenes against sig module trait correlations ####
ModuleME_Info_MELZ$hypoxia = TraitsME_LZ$O2
ModuleME_Info_MELZ$fetalMass = TraitsME_LZ$Fetus
ModuleME_Info_MELZ$Sex = TraitsME_LZ$Sex


#turquoise
ggplot(data = ModuleME_Info_MELZ, aes(x = ModuleME_Info_MELZ$fetalMass, y = ModuleME_Info_MELZ$MEturquoise)) + 
  geom_point(aes(color = ModuleME_Info_MELZ$hypoxia, fill =  ModuleME_Info_MELZ$hypoxia)) + theme_classic() + 
  stat_smooth(method = "lm",
              formula = y ~ x, aes(group = hypoxia, fill = hypoxia),
              geom = "smooth", color = "black") + theme_classic() + 
  scale_fill_manual(values = c("red", "blue")) + 
  scale_color_manual(values = c("red", "blue")) + 
  guides(color = guide_legend(title = "Oxygen")) + 
  guides(fill = guide_legend(title = "Oxygen")) + 
  xlab("Fetal Mass") +
  ylab("Turquoise Module Eigengene") 

#yellow
ggplot(data = ModuleME_Info_MELZ, aes(x = ModuleME_Info_MELZ$fetalMass, y = ModuleME_Info_MELZ$MEyellow)) + 
  geom_point(aes(color = ModuleME_Info_MELZ$hypoxia, fill =  ModuleME_Info_MELZ$hypoxia)) + theme_classic() + 
  stat_smooth(method = "lm",
              formula = y ~ x, aes(group = hypoxia, fill = hypoxia),
              geom = "smooth", color = "black") + theme_classic() + 
  scale_fill_manual(values = c("red", "blue")) + 
  scale_color_manual(values = c("red", "blue")) + 
  guides(color = guide_legend(title = "Oxygen")) + 
  guides(fill = guide_legend(title = "Oxygen")) + 
  xlab("Fetal Mass") +
  ylab("Yellow Module Eigengene")

# Create a dataset containing all gene-specific information
genes=names(ExprData_MELZ)
geneInfoMELZ = data.frame(Gene = genes,
                          moduleColor = moduleColorsMELZ)

# Network Preservation analysis ####
# Number of data sets that we work with
nSets = 2;
# Object that will contain the expression data
multiExprLZ = list();
multiExprLZ[[1]] = list(data = ExprData_BWLZ);
multiExprLZ[[2]] = list(data = ExprData_MELZ);
# Names for the two sets
setLabels = c("Lowlander", "Highlander");
# Important: components of multiExpr must carry identificating names
names(multiExprLZ) = setLabels
# Display the dimensions of the expression data (if you are confused by this construct, ignore it):
lapply(multiExprLZ, lapply, dim)

# Create an object (list) holding the module labels for each set:
colorList = list(moduleColorsBWLZ, moduleColorsMELZ)
# Components of the list must be named so that the names can be matched to the names of multiExpr
names(colorList) = setLabels

# Calculation of module preservation statistics
system.time( {
  mpLZ = WGCNA::modulePreservation(multiExprLZ, colorList,
                          referenceNetworks = 1,
                          loadPermutedStatistics = FALSE,
                          nPermutations = 200,
                          verbose = 3, 
                          maxModuleSize = 2869) #set max size to max in lowland network
} )

save(mpLZ, file = "Labyrinth Zone/LZHighLander-Lowlander_modPreservation.RData")
load(file = "Labyrinth Zone/LZHighLander-Lowlander_modPreservation.RData")

library(impute)
# Impute missing data and calculate eigengenes
impExprLZ = list();
for (set in 1:nSets)
{
  impExprLZ[[set]] = list(data = t(impute.knn(t(multiExprLZ[[set]]$data))$data));
}
eigengenes = list();
for (set in 1:nSets)
{
  eigengenes[[set]] = multiSetMEs(impExprLZ, universalColors = colorList[[set]], excludeGrey = TRUE);
  for (ss in 1:nSets)
  {
    rownames(eigengenes[[set]][[ss]]$data) = rownames(multiExprLZ[[ss]]$data);
  }
}

## Analysis and display of module preservation results

#  Isolate the observed statistics and their Z scores:
# Load the module preservation statistics
ref = 1 # Select the lowland data as reference
test = 2 # Select the highland data as test
statsObsLZ = cbind(mpLZ$quality$observed[[ref]][[test]][, -1], mpLZ$preservation$observed[[ref]][[test]][, -1])
statsZLZ = cbind(mpLZ$quality$Z[[ref]][[test]][, -1], mpLZ$preservation$Z[[ref]][[test]][, -1])
Z.PreservationStatsLZ = mpLZ$preservation$Z[[ref]][[test]]
Obs.PreservationStatsLZ = mpLZ$preservation$observed[[ref]][[test]]
# look at the main output: the preservation Zsummary score.
print(signif(statsZLZ[, "Zsummary.pres", drop = FALSE],2))
# Compare preservation to quality:
print(signif(statsZLZ[, c("Zsummary.pres", "Zsummary.qual")], 2))

# look at the main output: the preservation Zsummary score.
print(signif(statsObsLZ[, "Zsummary.pres", drop = FALSE],2))
# Compare preservation to quality:
print(signif(statsObsLZ[, c("Zsummary.pres", "Zsummary.qual")], 2))
statsObsLZ
statsZLZ


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
modColorsLZ = rownames(Obs.PreservationStatsLZ)
moduleSizeLZ = Obs.PreservationStatsLZ$moduleSize
# we will omit the grey module (background genes)
# and the gold module (random sample of genes)
selectModules = !(modColorsLZ %in% c("grey", "gold"))
# Text labels for points
point.labelLZ = modColorsLZ[selectModules]

#Composite preservation statistics
medianRankLZ = Obs.PreservationStatsLZ$medianRank.pres
ZsummaryLZ = Z.PreservationStatsLZ$Zsummary.pres

### medianRanks and Zsummary Plot ####
par(mfrow=c(1,2),mar = c(4.5,4.5,2.5,1))
# plot medianRank versus module size: The is useful for comparing relative 
# preservation among multiple modules: a module with lower median rank tends to 
# exhibit stronger observed preservation statistics than a module with a higher 
# median rank.
plot(moduleSizeLZ[selectModules],medianRankLZ[selectModules],col=1,
     bg=modColorsLZ[selectModules],pch = 21,main="medianRank Preservation",
     cex = 1.5, ylab ="medianRank",xlab="Module size", log="x")
labelPoints(moduleSizeLZ[selectModules],medianRankLZ[selectModules],point.labelLZ,cex=1,offs=0.03)

# plot Zsummary versus module size: If Zsummary is greater than 10 it means the module is preserved
# if 2<Zsummary<10 there is weak to moderate preservation, if less than 2
# there is no evidence
plot(moduleSizeLZ[selectModules],ZsummaryLZ[selectModules], col = 1,
     bg=modColorsLZ[selectModules],pch = 21,main="LZ Zsummary Preservation",
     cex=1.5,ylab ="Zsummary", xlab = "Module size", log = "x", ylim = c(0,100))
labelPoints(moduleSizeLZ[selectModules],ZsummaryLZ[selectModules],point.labelLZ,cex=1,offs=0.01)
# Add threshold lines for Zsummary
abline(h=0); abline(h=2, col = "blue", lty = 2); abline(h=10, col = "red", lty = 2)


#### Add Zsummary scores to BWtables (i.e., preservation stats for highlanders)
BWpreservedStatsLZ = print(signif(statsZLZ[, "Zsummary.pres", drop = FALSE],2))

# GO Analysis for ME modules ####

musGeneinfo = read.csv("DE_LZ_allGenes_wPvals_revised.csv")

## Match and add gene names from Mus from Kates musGeneInfo
## Filter genes from module colors that are in musGeneInfo
library(gprofiler2)
background = as.vector(musGeneinfo$mus_genename)

# Turq
turqME = geneInfoMELZ %>% filter(moduleColor == "turquoise")
as.data.frame(turqME)
turqME = turqME[, 1:2]
#write.csv(turqME, "turqME.csv")

rownames(turqME) = turqME$Gene
intersection = musGeneinfo$geneID  %in% rownames(turqME)
turqME0 = musGeneinfo[which(intersection), ] 
turqMEGO = as.vector(turqME0$mus_genename)

turqGO = gost(turqMEGO,
              organism = "mmusculus",
              user_threshold = 0.05,
              custom_bg = background,
              ordered_query = TRUE,
              correction_method = "bonferroni")

turqGOResults = data.frame(Term.Name = turqGO$result$term_name, 
                           P.value = turqGO$result$p_value, 
                           Source = turqGO$result$source, 
                           Term.Size = turqGO$result$term_size, 
                           Precision = turqGO$result$precision)
write.csv(turqGOResults, 'Labyrinth Zone/GOResults/turqGOResults.csv')

# green
greenME = geneInfoMELZ %>% filter(moduleColor == "green")
as.data.frame(greenME)
greenME = greenME[, 1:2]

rownames(greenME) = greenME$Gene
intersection0 = musGeneinfo$geneID  %in% rownames(greenME)
greenME0 = musGeneinfo[which(intersection0), ] 
greenMEGO = as.vector(greenME0$mus_genename)

greenGO = gost(greenMEGO,
               organism = "mmusculus",
               user_threshold = 0.05,
               custom_bg = background,
               ordered_query = TRUE,
               correction_method = "bonferroni")

greenGOResults = data.frame(Term.Name = greenGO$result$term_name, 
                            P.value = greenGO$result$p_value, 
                            Source = greenGO$result$source, 
                            Term.Size = greenGO$result$term_size, 
                            Precision = greenGO$result$precision)
write.csv(greenGOResults, 'Labyrinth Zone/GOResults/greenGOResults.csv')

# blue
blueME = geneInfoMELZ %>% filter(moduleColor == "blue")
as.data.frame(blueME)
blueME = blueME[, 1:2]

rownames(blueME) = blueME$Gene
intersection1 = musGeneinfo$geneID  %in% rownames(blueME)
blueME0 = musGeneinfo[which(intersection1), ] 
blueMEGO = as.vector(blueME0$mus_genename)

blueGO = gost(blueMEGO,
              organism = "mmusculus",
              user_threshold = 0.05,
              custom_bg = background,
              ordered_query = TRUE,
              correction_method = "bonferroni")

blueGOResults = data.frame(Term.Name = blueGO$result$term_name, 
                           P.value = blueGO$result$p_value, 
                           Source = blueGO$result$source, 
                           Term.Size = blueGO$result$term_size, 
                           Precision = blueGO$result$precision)
write.csv(blueGOResults, 'Labyrinth Zone/GOResults/blueGOResults.csv')

# black
blackME = geneInfoMELZ %>% filter(moduleColor == "black")
as.data.frame(blackME)
blackME = blackME[, 1:2]

rownames(blackME) = blackME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(blackME)
blackME0 = musGeneinfo[which(intersection2), ] 
blackMEGO = as.vector(blackME0$mus_genename)

#black
blackGO = gost(blackMEGO,
               organism = "mmusculus",
               user_threshold = 0.05,
               custom_bg = background,
               ordered_query = TRUE,
               correction_method = "bonferroni")

blackGOResults = data.frame(Term.Name = blackGO$result$term_name, 
                            P.value = blackGO$result$p_value, 
                            Source = blackGO$result$source, 
                            Term.Size = blackGO$result$term_size, 
                            Precision = blackGO$result$precision)
write.csv(blackGOResults, 'Labyrinth Zone/GOResults/blackGOResults.csv')

#brown
brownME = geneInfoMELZ %>% filter(moduleColor == "brown")
as.data.frame(brownME)
brownME = brownME[, 1:2]
#write.csv(brownME, "GOResults/Highlanders//blackME.csv")

rownames(brownME) = brownME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(brownME)
brownME0 = musGeneinfo[which(intersection2), ] 
brownMEGO = as.vector(brownME0$mus_genename)

brownMEGO = gost(brownMEGO,
                 organism = "mmusculus",
                 user_threshold = 0.05,
                 custom_bg = background,
                 ordered_query = TRUE,
                 correction_method = "bonferroni")

brownGOResults = data.frame(Term.Name = brownMEGO$result$term_name, 
                            P.value = brownMEGO$result$p_value, 
                            Source = brownMEGO$result$source, 
                            Term.Size = brownMEGO$result$term_size, 
                            Precision = brownMEGO$result$precision)
write.csv(brownGOResults, 'Labyrinth Zone/GOResults/brownGOResults.csv')

#greenyellow
greenyellowME = geneInfoMELZ %>% filter(moduleColor == "greenyellow")
as.data.frame(greenyellowME)
greenyellowME = greenyellowME[, 1:2]
#write.csv(greenyellowME, "GOResults/Highlanders/greenyellowME.csv")

rownames(greenyellowME) = greenyellowME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(greenyellowME)
greenyellowME0 = musGeneinfo[which(intersection2), ] 
greenyellowMEGO = as.vector(greenyellowME0$mus_genename)

greenyellowMEGO = gost(greenyellowMEGO,
                       organism = "mmusculus",
                       user_threshold = 0.05,
                       custom_bg = background,
                       ordered_query = TRUE,
                       correction_method = "bonferroni")

greenyellowGOResults = data.frame(Term.Name = greenyellowMEGO$result$term_name, 
                                  P.value = greenyellowMEGO$result$p_value, 
                                  Source = greenyellowMEGO$result$source, 
                                  Term.Size = greenyellowMEGO$result$term_size, 
                                  Precision = greenyellowMEGO$result$precision)
write.csv(greenyellowGOResults, 'Labyrinth Zone/GOResults/greenyellowGOResults.csv')

#magenta
magentaME = geneInfoMELZ %>% filter(moduleColor == "magenta")
as.data.frame(magentaME)
magentaME = magentaME[, 1:2]
#write.csv(magentaME, "GOResults/Highlanders/magentaME.csv")

rownames(magentaME) = magentaME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(magentaME)
magentaME0 = musGeneinfo[which(intersection2), ] 
magentaMEGO = as.vector(magentaME0$mus_genename)

magentaMEGO = gost(magentaMEGO,
                   organism = "mmusculus",
                   user_threshold = 0.05,
                   custom_bg = background,
                   ordered_query = TRUE,
                   correction_method = "bonferroni")

magentaGOResults = data.frame(Term.Name = magentaMEGO$result$term_name, 
                              P.value = magentaMEGO$result$p_value, 
                              Source = magentaMEGO$result$source, 
                              Term.Size = magentaMEGO$result$term_size, 
                              Precision = magentaMEGO$result$precision)
write.csv(magentaGOResults, 'Labyrinth Zone/GOResults/magentaGOResults.csv')

#pink
pinkME = geneInfoMELZ %>% filter(moduleColor == "pink")
as.data.frame(pinkME)
pinkME = pinkME[, 1:2]
#write.csv(pinkME, "GOResults/Highlanders/pinkME.csv")

rownames(pinkME) = pinkME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(pinkME)
pinkME0 = musGeneinfo[which(intersection2), ] 
pinkMEGO = as.vector(pinkME0$mus_genename)

pinkMEGO = gost(pinkMEGO,
                organism = "mmusculus",
                user_threshold = 0.05,
                custom_bg = background,
                ordered_query = TRUE,
                correction_method = "bonferroni")

pinkGOResults = data.frame(Term.Name = pinkMEGO$result$term_name, 
                           P.value = pinkMEGO$result$p_value, 
                           Source = pinkMEGO$result$source, 
                           Term.Size = pinkMEGO$result$term_size, 
                           Precision = pinkMEGO$result$precision)
write.csv(pinkGOResults, 'Labyrinth Zone/GOResults/pinkGOResults.csv')

#red
redME = geneInfoMELZ %>% filter(moduleColor == "red")
as.data.frame(redME)
redME = redME[, 1:2]
#write.csv(redME, "GOResults/Highlanders/redME.csv")

rownames(redME) = redME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(redME)
redME0 = musGeneinfo[which(intersection2), ] 
redMEGO = as.vector(redME0$mus_genename)

redMEGO = gost(redMEGO,
               organism = "mmusculus",
               user_threshold = 0.05,
               custom_bg = background,
               ordered_query = TRUE,
               correction_method = "bonferroni")

redGOResults = data.frame(Term.Name = redMEGO$result$term_name, 
                          P.value = redMEGO$result$p_value, 
                          Source = redMEGO$result$source, 
                          Term.Size = redMEGO$result$term_size, 
                          Precision = redMEGO$result$precision)
write.csv(redGOResults, 'Labyrinth Zone/GOResults/redGOResults.csv')

#yellow
yellowME = geneInfoMELZ %>% filter(moduleColor == "yellow")
as.data.frame(yellowME)
yellowME = yellowME[, 1:2]
#write.csv(yellowME, "GOResults/Highlanders/yellowME.csv")

rownames(yellowME) = yellowME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(yellowME)
yellowME0 = musGeneinfo[which(intersection2), ] 
yellowMEGO = as.vector(yellowME0$mus_genename)

yellowMEGO = gost(yellowMEGO,
               organism = "mmusculus",
               user_threshold = 0.05,
               custom_bg = background,
               ordered_query = TRUE,
               correction_method = "bonferroni")

yellowGOResults = data.frame(Term.Name = yellowMEGO$result$term_name, 
                          P.value = yellowMEGO$result$p_value, 
                          Source = yellowMEGO$result$source, 
                          Term.Size = yellowMEGO$result$term_size, 
                          Precision = yellowMEGO$result$precision)
write.csv(yellowGOResults, 'Labyrinth Zone/GOResults/yellowGOResults.csv')

# salmon
salmonME = geneInfoMELZ %>% filter(moduleColor == "salmon")
as.data.frame(salmonME)
salmonME = salmonME[, 1:2]
#write.csv(yellowME, "GOResults/Highlanders/yellowME.csv")

rownames(salmonME) = salmonME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(salmonME)
salmonME0 = musGeneinfo[which(intersection2), ] 
salmonMEGO = as.vector(salmonME0$mus_genename)

salmonMEGO = gost(salmonMEGO,
                  organism = "mmusculus",
                  user_threshold = 0.05,
                  custom_bg = background,
                  ordered_query = TRUE,
                  correction_method = "bonferroni")

salmonGOResults = data.frame(Term.Name = salmonMEGO$result$term_name, 
                             P.value = salmonMEGO$result$p_value, 
                             Source = salmonMEGO$result$source, 
                             Term.Size = salmonMEGO$result$term_size, 
                             Precision = salmonMEGO$result$precision)
write.csv(salmonGOResults, 'Labyrinth Zone/GOResults/salmonGOResults.csv')

#purple 
purpleME = geneInfoMELZ %>% filter(moduleColor == "purple")
as.data.frame(purpleME)
purpleME = purpleME[, 1:2]
#write.csv(yellowME, "GOResults/Highlanders/yellowME.csv")

rownames(purpleME) = purpleME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(purpleME)
purpleME0 = musGeneinfo[which(intersection2), ] 
purpleMEGO = as.vector(purpleME0$mus_genename)

purpleMEGO = gost(purpleMEGO,
                  organism = "mmusculus",
                  user_threshold = 0.05,
                  custom_bg = background,
                  ordered_query = TRUE,
                  correction_method = "bonferroni")

purpleGOResults = data.frame(Term.Name = purpleMEGO$result$term_name, 
                             P.value = purpleMEGO$result$p_value, 
                             Source = purpleMEGO$result$source, 
                             Term.Size = purpleMEGO$result$term_size, 
                             Precision = purpleMEGO$result$precision)
write.csv(purpleGOResults, 'Labyrinth Zone/GOResults/purpleGOResults.csv')

#lightcyan
lightcyanME = geneInfoMELZ %>% filter(moduleColor == "lightcyan")
as.data.frame(lightcyanME)
lightcyanME = lightcyanME[, 1:2]
#write.csv(yellowME, "GOResults/Highlanders/yellowME.csv")

rownames(lightcyanME) = lightcyanME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(lightcyanME)
lightcyanME0 = musGeneinfo[which(intersection2), ] 
lightcyanMEGO = as.vector(lightcyanME0$mus_genename)

lightcyanMEGO = gost(lightcyanMEGO,
                  organism = "mmusculus",
                  user_threshold = 0.05,
                  custom_bg = background,
                  ordered_query = TRUE,
                  correction_method = "bonferroni")

lightcyanGOResults = data.frame(Term.Name = lightcyanMEGO$result$term_name, 
                             P.value = lightcyanMEGO$result$p_value, 
                             Source = lightcyanMEGO$result$source, 
                             Term.Size = lightcyanMEGO$result$term_size, 
                             Precision = lightcyanMEGO$result$precision)
write.csv(lightcyanGOResults, 'Labyrinth Zone/GOResults/lightcyanGOResults.csv')

#lightyellow
lightyellowME = geneInfoMELZ %>% filter(moduleColor == "lightyellow")
as.data.frame(lightyellowME)
lightyellowME = lightyellowME[, 1:2]
#write.csv(yellowME, "GOResults/Highlanders/yellowME.csv")

rownames(lightyellowME) = lightyellowME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(lightyellowME)
lightyellowME0 = musGeneinfo[which(intersection2), ] 
lightyellowMEGO = as.vector(lightyellowME0$mus_genename)

lightyellowMEGO = gost(lightyellowMEGO,
                     organism = "mmusculus",
                     user_threshold = 0.05,
                     custom_bg = background,
                     ordered_query = TRUE,
                     correction_method = "bonferroni")

lightyellowGOResults = data.frame(Term.Name = lightyellowMEGO$result$term_name, 
                                P.value = lightyellowMEGO$result$p_value, 
                                Source = lightyellowMEGO$result$source, 
                                Term.Size = lightyellowMEGO$result$term_size, 
                                Precision = lightyellowMEGO$result$precision)
write.csv(lightyellowGOResults, 'Labyrinth Zone/GOResults/lightyellowGOResults.csv')

#royalblue
royalblueME = geneInfoMELZ %>% filter(moduleColor == "royalblue")
as.data.frame(royalblueME)
royalblueME = royalblueME[, 1:2]
#write.csv(yellowME, "GOResults/Highlanders/yellowME.csv")

rownames(royalblueME) = royalblueME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(royalblueME)
royalblueME0 = musGeneinfo[which(intersection2), ] 
royalblueMEGO = as.vector(royalblueME0$mus_genename)

royalblueMEGO = gost(royalblueMEGO,
                       organism = "mmusculus",
                       user_threshold = 0.05,
                       custom_bg = background,
                       ordered_query = TRUE,
                       correction_method = "bonferroni")

royalblueGOResults = data.frame(Term.Name = royalblueMEGO$result$term_name, 
                                  P.value = royalblueMEGO$result$p_value, 
                                  Source = royalblueMEGO$result$source, 
                                  Term.Size = royalblueMEGO$result$term_size, 
                                  Precision = royalblueMEGO$result$precision)
write.csv(royalblueGOResults, 'Labyrinth Zone/GOResults/royalblueGOResults.csv')

#darkgreen
darkgreenME = geneInfoMELZ %>% filter(moduleColor == "darkgreen")
as.data.frame(darkgreenME)
darkgreenME = darkgreenME[, 1:2]
#write.csv(yellowME, "GOResults/Highlanders/yellowME.csv")

rownames(darkgreenME) = darkgreenME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(darkgreenME)
darkgreenME0 = musGeneinfo[which(intersection2), ] 
darkgreenMEGO = as.vector(darkgreenME0$mus_genename)

darkgreenMEGO = gost(darkgreenMEGO,
                     organism = "mmusculus",
                     user_threshold = 0.05,
                     custom_bg = background,
                     ordered_query = TRUE,
                     correction_method = "bonferroni")

darkgreenGOResults = data.frame(Term.Name = darkgreenMEGO$result$term_name, 
                                P.value = darkgreenMEGO$result$p_value, 
                                Source = darkgreenMEGO$result$source, 
                                Term.Size = darkgreenMEGO$result$term_size, 
                                Precision = darkgreenMEGO$result$precision)
write.csv(darkgreenGOResults, 'Labyrinth Zone/GOResults/darkgreenGOResults.csv')

#grey60
grey60ME = geneInfoMELZ %>% filter(moduleColor == "grey60")
as.data.frame(grey60ME)
grey60ME = grey60ME[, 1:2]
#write.csv(yellowME, "GOResults/Highlanders/yellowME.csv")

rownames(grey60ME) = grey60ME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(grey60ME)
grey60ME0 = musGeneinfo[which(intersection2), ] 
grey60MEGO = as.vector(grey60ME0$mus_genename)

grey60MEGO = gost(grey60MEGO,
                     organism = "mmusculus",
                     user_threshold = 0.05,
                     custom_bg = background,
                     ordered_query = TRUE,
                     correction_method = "bonferroni")

grey60GOResults = data.frame(Term.Name = grey60MEGO$result$term_name, 
                                P.value = grey60MEGO$result$p_value, 
                                Source = grey60MEGO$result$source, 
                                Term.Size = grey60MEGO$result$term_size, 
                                Precision = grey60MEGO$result$precision)
write.csv(grey60GOResults, 'Labyrinth Zone/GOResults/grey60GOResults.csv')

#darkred
darkredME = geneInfoMELZ %>% filter(moduleColor == "darkred")
as.data.frame(darkredME)
darkredME = darkredME[, 1:2]
#write.csv(yellowME, "GOResults/Highlanders/yellowME.csv")

rownames(darkredME) = darkredME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(darkredME)
darkredME0 = musGeneinfo[which(intersection2), ] 
darkredMEGO = as.vector(darkredME0$mus_genename)

darkredMEGO = gost(darkredMEGO,
                  organism = "mmusculus",
                  user_threshold = 0.05,
                  custom_bg = background,
                  ordered_query = TRUE,
                  correction_method = "bonferroni")

darkredGOResults = data.frame(Term.Name = darkredMEGO$result$term_name, 
                             P.value = darkredMEGO$result$p_value, 
                             Source = darkredMEGO$result$source, 
                             Term.Size = darkredMEGO$result$term_size, 
                             Precision = darkredMEGO$result$precision)
write.csv(darkredGOResults, 'Labyrinth Zone/GOResults/darkredGOResults.csv')

#lightgreen
lightgreenME = geneInfoMELZ %>% filter(moduleColor == "lightgreen")
as.data.frame(lightgreenME)
lightgreenME = lightgreenME[, 1:2]
#write.csv(yellowME, "GOResults/Highlanders/yellowME.csv")

rownames(lightgreenME) = lightgreenME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(lightgreenME)
lightgreenME0 = musGeneinfo[which(intersection2), ] 
lightgreenMEGO = as.vector(lightgreenME0$mus_genename)

lightgreenMEGO = gost(lightgreenMEGO,
                   organism = "mmusculus",
                   user_threshold = 0.05,
                   custom_bg = background,
                   ordered_query = TRUE,
                   correction_method = "bonferroni")

lightgreenGOResults = data.frame(Term.Name = lightgreenMEGO$result$term_name, 
                              P.value = lightgreenMEGO$result$p_value, 
                              Source = lightgreenMEGO$result$source, 
                              Term.Size = lightgreenMEGO$result$term_size, 
                              Precision = lightgreenMEGO$result$precision)
write.csv(lightgreenGOResults, 'Labyrinth Zone/GOResults/lightgreenGOResults.csv')

#cyan
cyanME = geneInfoMELZ %>% filter(moduleColor == "cyan")
as.data.frame(cyanME)
cyanME = cyanME[, 1:2]
#write.csv(yellowME, "GOResults/Highlanders/yellowME.csv")

rownames(cyanME) = cyanME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(cyanME)
cyanME0 = musGeneinfo[which(intersection2), ] 
cyanMEGO = as.vector(cyanME0$mus_genename)

cyanMEGO = gost(cyanMEGO,
                      organism = "mmusculus",
                      user_threshold = 0.05,
                      custom_bg = background,
                      ordered_query = TRUE,
                      correction_method = "bonferroni")

cyanGOResults = data.frame(Term.Name = cyanMEGO$result$term_name, 
                                 P.value = cyanMEGO$result$p_value, 
                                 Source = cyanMEGO$result$source, 
                                 Term.Size = cyanMEGO$result$term_size, 
                                 Precision = cyanMEGO$result$precision)
write.csv(cyanGOResults, 'Labyrinth Zone/GOResults/cyanGOResults.csv')

#midnightblue
midnightblueME = geneInfoMELZ %>% filter(moduleColor == "midnightblue")
as.data.frame(midnightblueME)
midnightblueME = midnightblueME[, 1:2]
#write.csv(yellowME, "GOResults/Highlanders/yellowME.csv")

rownames(midnightblueME) = midnightblueME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(midnightblueME)
midnightblueME0 = musGeneinfo[which(intersection2), ] 
midnightblueMEGO = as.vector(midnightblueME0$mus_genename)

midnightblueMEGO = gost(midnightblueMEGO,
                organism = "mmusculus",
                user_threshold = 0.05,
                custom_bg = background,
                ordered_query = TRUE,
                correction_method = "bonferroni")

midnightblueGOResults = data.frame(Term.Name = midnightblueMEGO$result$term_name, 
                           P.value = midnightblueMEGO$result$p_value, 
                           Source = midnightblueMEGO$result$source, 
                           Term.Size = midnightblueMEGO$result$term_size, 
                           Precision = midnightblueMEGO$result$precision)
write.csv(midnightblueGOResults, 'Labyrinth Zone/GOResults/midnightblueGOResults.csv')

#tan
tanME = geneInfoMELZ %>% filter(moduleColor == "tan")
as.data.frame(tanME)
tanME = tanME[, 1:2]
#write.csv(yellowME, "GOResults/Highlanders/yellowME.csv")

rownames(tanME) = tanME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(tanME)
tanME0 = musGeneinfo[which(intersection2), ] 
tanMEGO = as.vector(tanME0$mus_genename)

tanMEGO = gost(tanMEGO,
                        organism = "mmusculus",
                        user_threshold = 0.05,
                        custom_bg = background,
                        ordered_query = TRUE,
                        correction_method = "bonferroni")

tanGOResults = data.frame(Term.Name = tanMEGO$result$term_name, 
                                   P.value = tanMEGO$result$p_value, 
                                   Source = tanMEGO$result$source, 
                                   Term.Size = tanMEGO$result$term_size, 
                                   Precision = tanMEGO$result$precision)
write.csv(tanGOResults, 'Labyrinth Zone/GOResults/tanGOResults.csv')

#darkgrey
darkgreyME = geneInfoMELZ %>% filter(moduleColor == "darkgrey")
as.data.frame(darkgreyME)
darkgreyME = darkgreyME[, 1:2]
#write.csv(yellowME, "GOResults/Highlanders/yellowME.csv")

rownames(darkgreyME) = darkgreyME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(darkgreyME)
darkgreyME0 = musGeneinfo[which(intersection2), ] 
darkgreyMEGO = as.vector(darkgreyME0$mus_genename)

darkgreyMEGO = gost(darkgreyMEGO,
               organism = "mmusculus",
               user_threshold = 0.05,
               custom_bg = background,
               ordered_query = TRUE,
               correction_method = "bonferroni")

darkgreyGOResults = data.frame(Term.Name = darkgreyMEGO$result$term_name, 
                          P.value = darkgreyMEGO$result$p_value, 
                          Source = darkgreyMEGO$result$source, 
                          Term.Size = darkgreyMEGO$result$term_size, 
                          Precision = darkgreyMEGO$result$precision)
write.csv(darkgreyGOResults, 'Labyrinth Zone/GOResults/Highlanders/darkgreyGOResults.csv')

#darkolivegreen
darkolivegreenME = geneInfoMELZ %>% filter(moduleColor == "darkolivegreen")
as.data.frame(darkolivegreenME)
darkolivegreenME = darkolivegreenME[, 1:2]
#write.csv(yellowME, "GOResults/Highlanders/yellowME.csv")

rownames(darkolivegreenME) = darkolivegreenME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(darkolivegreenME)
darkolivegreenME0 = musGeneinfo[which(intersection2), ] 
darkolivegreenMEGO = as.vector(darkolivegreenME0$mus_genename)

darkolivegreenMEGO = gost(darkolivegreenMEGO,
                    organism = "mmusculus",
                    user_threshold = 0.05,
                    custom_bg = background,
                    ordered_query = TRUE,
                    correction_method = "bonferroni")

darkolivegreenGOResults = data.frame(Term.Name = darkolivegreenMEGO$result$term_name, 
                               P.value = darkolivegreenMEGO$result$p_value, 
                               Source = darkolivegreenMEGO$result$source, 
                               Term.Size = darkolivegreenMEGO$result$term_size, 
                               Precision = darkolivegreenMEGO$result$precision)
write.csv(darkolivegreenGOResults, 'Labyrinth Zone/GOResults/Highlanders/darkolivegreenGOResults.csv')

#darkorange
darkorangeME = geneInfoMELZ %>% filter(moduleColor == "darkorange")
as.data.frame(darkorangeME)
darkorangeME = darkorangeME[, 1:2]
#write.csv(yellowME, "GOResults/Highlanders/yellowME.csv")

rownames(darkorangeME) = darkorangeME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(darkorangeME)
darkorangeME0 = musGeneinfo[which(intersection2), ] 
darkorangeMEGO = as.vector(darkorangeME0$mus_genename)

darkorangeMEGO = gost(darkorangeMEGO,
                          organism = "mmusculus",
                          user_threshold = 0.05,
                          custom_bg = background,
                          ordered_query = TRUE,
                          correction_method = "bonferroni")

darkorangeGOResults = data.frame(Term.Name = darkorangeMEGO$result$term_name, 
                                     P.value = darkorangeMEGO$result$p_value, 
                                     Source = darkorangeMEGO$result$source, 
                                     Term.Size = darkorangeMEGO$result$term_size, 
                                     Precision = darkorangeMEGO$result$precision)
write.csv(darkorangeGOResults, 'Labyrinth Zone/GOResults/Highlanders/darkorangeGOResults.csv')

#darkturquoise
darkturquoiseME = geneInfoMELZ %>% filter(moduleColor == "darkturquoise")
as.data.frame(darkturquoiseME)
darkturquoiseME = darkturquoiseME[, 1:2]
#write.csv(yellowME, "GOResults/Highlanders/yellowME.csv")

rownames(darkturquoiseME) = darkturquoiseME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(darkturquoiseME)
darkturquoiseME0 = musGeneinfo[which(intersection2), ] 
darkturquoiseMEGO = as.vector(darkturquoiseME0$mus_genename)

darkturquoiseMEGO = gost(darkturquoiseMEGO,
                      organism = "mmusculus",
                      user_threshold = 0.05,
                      custom_bg = background,
                      ordered_query = TRUE,
                      correction_method = "bonferroni")

darkturquoiseGOResults = data.frame(Term.Name = darkturquoiseMEGO$result$term_name, 
                                 P.value = darkturquoiseMEGO$result$p_value, 
                                 Source = darkturquoiseMEGO$result$source, 
                                 Term.Size = darkturquoiseMEGO$result$term_size, 
                                 Precision = darkturquoiseMEGO$result$precision)
write.csv(darkturquoiseGOResults, 'Labyrinth Zone/GOResults/Highlanders/darkturquoiseGOResults.csv')


#paleturquoise
paleturquoiseME = geneInfoMELZ %>% filter(moduleColor == "paleturquoise")
as.data.frame(paleturquoiseME)
paleturquoiseME = paleturquoiseME[, 1:2]
#write.csv(yellowME, "GOResults/Highlanders/yellowME.csv")

rownames(paleturquoiseME) = paleturquoiseME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(paleturquoiseME)
paleturquoiseME0 = musGeneinfo[which(intersection2), ] 
paleturquoiseMEGO = as.vector(paleturquoiseME0$mus_genename)

paleturquoiseMEGO = gost(paleturquoiseMEGO,
                         organism = "mmusculus",
                         user_threshold = 0.05,
                         custom_bg = background,
                         ordered_query = TRUE,
                         correction_method = "bonferroni")

paleturquoiseGOResults = data.frame(Term.Name = paleturquoiseMEGO$result$term_name, 
                                    P.value = paleturquoiseMEGO$result$p_value, 
                                    Source = paleturquoiseMEGO$result$source, 
                                    Term.Size = paleturquoiseMEGO$result$term_size, 
                                    Precision = paleturquoiseMEGO$result$precision)
write.csv(paleturquoiseGOResults, 'Labyrinth Zone/GOResults/Highlanders/paleturquoiseGOResults.csv')

#saddlebrown
saddlebrownME = geneInfoMELZ %>% filter(moduleColor == "saddlebrown")
as.data.frame(saddlebrownME)
saddlebrownME = saddlebrownME[, 1:2]
#write.csv(yellowME, "GOResults/Highlanders/yellowME.csv")

rownames(saddlebrownME) = saddlebrownME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(saddlebrownME)
saddlebrownME0 = musGeneinfo[which(intersection2), ] 
saddlebrownMEGO = as.vector(saddlebrownME0$mus_genename)

saddlebrownMEGO = gost(saddlebrownMEGO,
                         organism = "mmusculus",
                         user_threshold = 0.05,
                         custom_bg = background,
                         ordered_query = TRUE,
                         correction_method = "bonferroni")

saddlebrownGOResults = data.frame(Term.Name = saddlebrownMEGO$result$term_name, 
                                    P.value = saddlebrownMEGO$result$p_value, 
                                    Source = saddlebrownMEGO$result$source, 
                                    Term.Size = saddlebrownMEGO$result$term_size, 
                                    Precision = saddlebrownMEGO$result$precision)
write.csv(saddlebrownGOResults, 'Labyrinth Zone/GOResults/Highlanders/saddlebrownGOResults.csv')

#skyblue
skyblueME = geneInfoMELZ %>% filter(moduleColor == "skyblue")
as.data.frame(skyblueME)
skyblueME = skyblueME[, 1:2]
#write.csv(yellowME, "GOResults/Highlanders/yellowME.csv")

rownames(skyblueME) = skyblueME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(skyblueME)
skyblueME0 = musGeneinfo[which(intersection2), ] 
skyblueMEGO = as.vector(skyblueME0$mus_genename)

skyblueMEGO = gost(skyblueMEGO,
                       organism = "mmusculus",
                       user_threshold = 0.05,
                       custom_bg = background,
                       ordered_query = TRUE,
                       correction_method = "bonferroni")

skyblueGOResults = data.frame(Term.Name = skyblueMEGO$result$term_name, 
                                  P.value = skyblueMEGO$result$p_value, 
                                  Source = skyblueMEGO$result$source, 
                                  Term.Size = skyblueMEGO$result$term_size, 
                                  Precision = skyblueMEGO$result$precision)
write.csv(skyblueGOResults, 'Labyrinth Zone/GOResults/Highlanders/skyblueGOResults.csv')

#steelblue
steelblueME = geneInfoMELZ %>% filter(moduleColor == "steelblue")
as.data.frame(steelblueME)
steelblueME = steelblueME[, 1:2]
#write.csv(yellowME, "GOResults/Highlanders/yellowME.csv")

rownames(steelblueME) = steelblueME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(steelblueME)
steelblueME0 = musGeneinfo[which(intersection2), ] 
steelblueMEGO = as.vector(steelblueME0$mus_genename)

steelblueMEGO = gost(steelblueMEGO,
                   organism = "mmusculus",
                   user_threshold = 0.05,
                   custom_bg = background,
                   ordered_query = TRUE,
                   correction_method = "bonferroni")

steelblueGOResults = data.frame(Term.Name = steelblueMEGO$result$term_name, 
                              P.value = steelblueMEGO$result$p_value, 
                              Source = steelblueMEGO$result$source, 
                              Term.Size = steelblueMEGO$result$term_size, 
                              Precision = steelblueMEGO$result$precision)
write.csv(steelblueGOResults, 'Labyrinth Zone/GOResults/Highlanders/steelblueGOResults.csv')

#violet
violetME = geneInfoMELZ %>% filter(moduleColor == "violet")
as.data.frame(violetME)
violetME = violetME[, 1:2]
#write.csv(yellowME, "GOResults/Highlanders/yellowME.csv")

rownames(violetME) = violetME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(violetME)
violetME0 = musGeneinfo[which(intersection2), ] 
violetMEGO = as.vector(violetME0$mus_genename)

violetMEGO = gost(violetMEGO,
                     organism = "mmusculus",
                     user_threshold = 0.05,
                     custom_bg = background,
                     ordered_query = TRUE,
                     correction_method = "bonferroni")

violetGOResults = data.frame(Term.Name = violetMEGO$result$term_name, 
                                P.value = violetMEGO$result$p_value, 
                                Source = violetMEGO$result$source, 
                                Term.Size = violetMEGO$result$term_size, 
                                Precision = violetMEGO$result$precision)
write.csv(violetGOResults, 'Labyrinth Zone/GOResults/Highlanders/violetGOResults.csv')

#white
whiteME = geneInfoMELZ %>% filter(moduleColor == "white")
as.data.frame(whiteME)
whiteME = whiteME[, 1:2]
#write.csv(yellowME, "GOResults/Highlanders/yellowME.csv")

rownames(whiteME) = whiteME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(whiteME)
whiteME0 = musGeneinfo[which(intersection2), ] 
whiteMEGO = as.vector(whiteME0$mus_genename)

whiteMEGO = gost(whiteMEGO,
                  organism = "mmusculus",
                  user_threshold = 0.05,
                  custom_bg = background,
                  ordered_query = TRUE,
                  correction_method = "bonferroni")

whiteGOResults = data.frame(Term.Name = whiteMEGO$result$term_name, 
                             P.value = whiteMEGO$result$p_value, 
                             Source = whiteMEGO$result$source, 
                             Term.Size = whiteMEGO$result$term_size, 
                             Precision = whiteMEGO$result$precision)
write.csv(whiteGOResults, 'Labyrinth Zone/GOResults/Highlanders/whiteGOResults.csv')



## Placental cell enrichment analysis ####

# GO Analysis for BW LZ animals ####
musGeneinfo = read.csv("DE_LZ_allGenes_wPvals_revised.csv")

## Match and add gene names from Mus from Kates musGeneInfo
## Filter genes from module colors that are in musGeneInfo
library(gprofiler2)
background = as.vector(musGeneinfo$mus_genename)

# Turq
turqME = geneInfoBWLZ %>% filter(moduleColor == "turquoise")
as.data.frame(turqME)
turqME = turqME[, 1:2]
#write.csv(turqME, "turqME.csv")

rownames(turqME) = turqME$Gene
intersection = musGeneinfo$geneID  %in% rownames(turqME)
turqME0 = musGeneinfo[which(intersection), ] 
turqMEGO = as.vector(turqME0$mus_genename)

turqGO = gost(turqMEGO,
              organism = "mmusculus",
              user_threshold = 0.05,
              custom_bg = background,
              ordered_query = TRUE,
              correction_method = "bonferroni")

turqGOResults = data.frame(Term.Name = turqGO$result$term_name, 
                           P.value = turqGO$result$p_value, 
                           Source = turqGO$result$source, 
                           Term.Size = turqGO$result$term_size, 
                           Precision = turqGO$result$precision)
write.csv(turqGOResults, 'Labyrinth Zone/GOResults/Lowlanders/turqGOResults.csv')

# green
greenME = geneInfoBWLZ %>% filter(moduleColor == "green")
as.data.frame(greenME)
greenME = greenME[, 1:2]

rownames(greenME) = greenME$Gene
intersection0 = musGeneinfo$geneID  %in% rownames(greenME)
greenME0 = musGeneinfo[which(intersection0), ] 
greenMEGO = as.vector(greenME0$mus_genename)

greenGO = gost(greenMEGO,
               organism = "mmusculus",
               user_threshold = 0.05,
               custom_bg = background,
               ordered_query = TRUE,
               correction_method = "bonferroni")

greenGOResults = data.frame(Term.Name = greenGO$result$term_name, 
                            P.value = greenGO$result$p_value, 
                            Source = greenGO$result$source, 
                            Term.Size = greenGO$result$term_size, 
                            Precision = greenGO$result$precision)
write.csv(greenGOResults, 'Labyrinth Zone/GOResults/Lowlanders/greenGOResults.csv')

# blue
blueME = geneInfoBWLZ %>% filter(moduleColor == "blue")
as.data.frame(blueME)
blueME = blueME[, 1:2]

rownames(blueME) = blueME$Gene
intersection1 = musGeneinfo$geneID  %in% rownames(blueME)
blueME0 = musGeneinfo[which(intersection1), ] 
blueMEGO = as.vector(blueME0$mus_genename)

blueGO = gost(blueMEGO,
              organism = "mmusculus",
              user_threshold = 0.05,
              custom_bg = background,
              ordered_query = TRUE,
              correction_method = "bonferroni")

blueGOResults = data.frame(Term.Name = blueGO$result$term_name, 
                           P.value = blueGO$result$p_value, 
                           Source = blueGO$result$source, 
                           Term.Size = blueGO$result$term_size, 
                           Precision = blueGO$result$precision)
write.csv(blueGOResults, 'Labyrinth Zone/GOResults/Lowlanders/blueGOResults.csv')

# black
blackME = geneInfoBWLZ %>% filter(moduleColor == "black")
as.data.frame(blackME)
blackME = blackME[, 1:2]

rownames(blackME) = blackME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(blackME)
blackME0 = musGeneinfo[which(intersection2), ] 
blackMEGO = as.vector(blackME0$mus_genename)

#black
blackGO = gost(blackMEGO,
               organism = "mmusculus",
               user_threshold = 0.05,
               custom_bg = background,
               ordered_query = TRUE,
               correction_method = "bonferroni")

blackGOResults = data.frame(Term.Name = blackGO$result$term_name, 
                            P.value = blackGO$result$p_value, 
                            Source = blackGO$result$source, 
                            Term.Size = blackGO$result$term_size, 
                            Precision = blackGO$result$precision)
write.csv(blackGOResults, 'Labyrinth Zone/GOResults/Lowlanders/blackGOResults.csv')

#brown
brownME = geneInfoBWLZ %>% filter(moduleColor == "brown")
as.data.frame(brownME)
brownME = brownME[, 1:2]
#write.csv(brownME, "GOResults/Highlanders//blackME.csv")

rownames(brownME) = brownME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(brownME)
brownME0 = musGeneinfo[which(intersection2), ] 
brownMEGO = as.vector(brownME0$mus_genename)

brownMEGO = gost(brownMEGO,
                 organism = "mmusculus",
                 user_threshold = 0.05,
                 custom_bg = background,
                 ordered_query = TRUE,
                 correction_method = "bonferroni")

brownGOResults = data.frame(Term.Name = brownMEGO$result$term_name, 
                            P.value = brownMEGO$result$p_value, 
                            Source = brownMEGO$result$source, 
                            Term.Size = brownMEGO$result$term_size, 
                            Precision = brownMEGO$result$precision)
write.csv(brownGOResults, 'Labyrinth Zone/GOResults/Lowlanders/brownGOResults.csv')

#greenyellow
greenyellowME = geneInfoBWLZ %>% filter(moduleColor == "greenyellow")
as.data.frame(greenyellowME)
greenyellowME = greenyellowME[, 1:2]
#write.csv(greenyellowME, "GOResults/Highlanders/greenyellowME.csv")

rownames(greenyellowME) = greenyellowME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(greenyellowME)
greenyellowME0 = musGeneinfo[which(intersection2), ] 
greenyellowMEGO = as.vector(greenyellowME0$mus_genename)

greenyellowMEGO = gost(greenyellowMEGO,
                       organism = "mmusculus",
                       user_threshold = 0.05,
                       custom_bg = background,
                       ordered_query = TRUE,
                       correction_method = "bonferroni")

greenyellowGOResults = data.frame(Term.Name = greenyellowMEGO$result$term_name, 
                                  P.value = greenyellowMEGO$result$p_value, 
                                  Source = greenyellowMEGO$result$source, 
                                  Term.Size = greenyellowMEGO$result$term_size, 
                                  Precision = greenyellowMEGO$result$precision)
write.csv(greenyellowGOResults, 'Labyrinth Zone/GOResults/Lowlanders/greenyellowGOResults.csv')

#magenta
magentaME = geneInfoBWLZ %>% filter(moduleColor == "magenta")
as.data.frame(magentaME)
magentaME = magentaME[, 1:2]
#write.csv(magentaME, "GOResults/Highlanders/magentaME.csv")

rownames(magentaME) = magentaME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(magentaME)
magentaME0 = musGeneinfo[which(intersection2), ] 
magentaMEGO = as.vector(magentaME0$mus_genename)

magentaMEGO = gost(magentaMEGO,
                   organism = "mmusculus",
                   user_threshold = 0.05,
                   custom_bg = background,
                   ordered_query = TRUE,
                   correction_method = "bonferroni")

magentaGOResults = data.frame(Term.Name = magentaMEGO$result$term_name, 
                              P.value = magentaMEGO$result$p_value, 
                              Source = magentaMEGO$result$source, 
                              Term.Size = magentaMEGO$result$term_size, 
                              Precision = magentaMEGO$result$precision)
write.csv(magentaGOResults, 'Labyrinth Zone/GOResults/Lowlanders/magentaGOResults.csv')

#pink
pinkME = geneInfoBWLZ %>% filter(moduleColor == "pink")
as.data.frame(pinkME)
pinkME = pinkME[, 1:2]
#write.csv(pinkME, "GOResults/Highlanders/pinkME.csv")

rownames(pinkME) = pinkME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(pinkME)
pinkME0 = musGeneinfo[which(intersection2), ] 
pinkMEGO = as.vector(pinkME0$mus_genename)

pinkMEGO = gost(pinkMEGO,
                organism = "mmusculus",
                user_threshold = 0.05,
                custom_bg = background,
                ordered_query = TRUE,
                correction_method = "bonferroni")

pinkGOResults = data.frame(Term.Name = pinkMEGO$result$term_name, 
                           P.value = pinkMEGO$result$p_value, 
                           Source = pinkMEGO$result$source, 
                           Term.Size = pinkMEGO$result$term_size, 
                           Precision = pinkMEGO$result$precision)
write.csv(pinkGOResults, 'Labyrinth Zone/GOResults/Lowlanders/pinkGOResults.csv')

#red
redME = geneInfoBWLZ %>% filter(moduleColor == "red")
as.data.frame(redME)
redME = redME[, 1:2]
#write.csv(redME, "GOResults/Highlanders/redME.csv")

rownames(redME) = redME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(redME)
redME0 = musGeneinfo[which(intersection2), ] 
redMEGO = as.vector(redME0$mus_genename)

redMEGO = gost(redMEGO,
               organism = "mmusculus",
               user_threshold = 0.05,
               custom_bg = background,
               ordered_query = TRUE,
               correction_method = "bonferroni")

redGOResults = data.frame(Term.Name = redMEGO$result$term_name, 
                          P.value = redMEGO$result$p_value, 
                          Source = redMEGO$result$source, 
                          Term.Size = redMEGO$result$term_size, 
                          Precision = redMEGO$result$precision)
write.csv(redGOResults, 'Labyrinth Zone/GOResults/Lowlanders/redGOResults.csv')

#yellow
yellowME = geneInfoBWLZ %>% filter(moduleColor == "yellow")
as.data.frame(yellowME)
yellowME = yellowME[, 1:2]
#write.csv(yellowME, "GOResults/Highlanders/yellowME.csv")

rownames(yellowME) = yellowME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(yellowME)
yellowME0 = musGeneinfo[which(intersection2), ] 
yellowMEGO = as.vector(yellowME0$mus_genename)

yellowMEGO = gost(yellowMEGO,
                  organism = "mmusculus",
                  user_threshold = 0.05,
                  custom_bg = background,
                  ordered_query = TRUE,
                  correction_method = "bonferroni")

yellowGOResults = data.frame(Term.Name = yellowMEGO$result$term_name, 
                             P.value = yellowMEGO$result$p_value, 
                             Source = yellowMEGO$result$source, 
                             Term.Size = yellowMEGO$result$term_size, 
                             Precision = yellowMEGO$result$precision)
write.csv(yellowGOResults, 'Labyrinth Zone/GOResults/Lowlanders/yellowGOResults.csv')

# salmon
salmonME = geneInfoBWLZ %>% filter(moduleColor == "salmon")
as.data.frame(salmonME)
salmonME = salmonME[, 1:2]
#write.csv(yellowME, "GOResults/Highlanders/yellowME.csv")

rownames(salmonME) = salmonME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(salmonME)
salmonME0 = musGeneinfo[which(intersection2), ] 
salmonMEGO = as.vector(salmonME0$mus_genename)

salmonMEGO = gost(salmonMEGO,
                  organism = "mmusculus",
                  user_threshold = 0.05,
                  custom_bg = background,
                  ordered_query = TRUE,
                  correction_method = "bonferroni")

salmonGOResults = data.frame(Term.Name = salmonMEGO$result$term_name, 
                             P.value = salmonMEGO$result$p_value, 
                             Source = salmonMEGO$result$source, 
                             Term.Size = salmonMEGO$result$term_size, 
                             Precision = salmonMEGO$result$precision)
write.csv(salmonGOResults, 'Labyrinth Zone/GOResults/Lowlanders/salmonGOResults.csv')

#purple 
purpleME = geneInfoBWLZ %>% filter(moduleColor == "purple")
as.data.frame(purpleME)
purpleME = purpleME[, 1:2]
#write.csv(yellowME, "GOResults/Highlanders/yellowME.csv")

rownames(purpleME) = purpleME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(purpleME)
purpleME0 = musGeneinfo[which(intersection2), ] 
purpleMEGO = as.vector(purpleME0$mus_genename)

purpleMEGO = gost(purpleMEGO,
                  organism = "mmusculus",
                  user_threshold = 0.05,
                  custom_bg = background,
                  ordered_query = TRUE,
                  correction_method = "bonferroni")

purpleGOResults = data.frame(Term.Name = purpleMEGO$result$term_name, 
                             P.value = purpleMEGO$result$p_value, 
                             Source = purpleMEGO$result$source, 
                             Term.Size = purpleMEGO$result$term_size, 
                             Precision = purpleMEGO$result$precision)
write.csv(purpleGOResults, 'Labyrinth Zone/GOResults/Lowlanders/purpleGOResults.csv')

#lightcyan
lightcyanME = geneInfoBWLZ %>% filter(moduleColor == "lightcyan")
as.data.frame(lightcyanME)
lightcyanME = lightcyanME[, 1:2]
#write.csv(yellowME, "GOResults/Highlanders/yellowME.csv")

rownames(lightcyanME) = lightcyanME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(lightcyanME)
lightcyanME0 = musGeneinfo[which(intersection2), ] 
lightcyanMEGO = as.vector(lightcyanME0$mus_genename)

lightcyanMEGO = gost(lightcyanMEGO,
                     organism = "mmusculus",
                     user_threshold = 0.05,
                     custom_bg = background,
                     ordered_query = TRUE,
                     correction_method = "bonferroni")

lightcyanGOResults = data.frame(Term.Name = lightcyanMEGO$result$term_name, 
                                P.value = lightcyanMEGO$result$p_value, 
                                Source = lightcyanMEGO$result$source, 
                                Term.Size = lightcyanMEGO$result$term_size, 
                                Precision = lightcyanMEGO$result$precision)
write.csv(lightcyanGOResults, 'Labyrinth Zone/GOResults/Lowlanders/lightcyanGOResults.csv')

#lightyellow
lightyellowME = geneInfoBWLZ %>% filter(moduleColor == "lightyellow")
as.data.frame(lightyellowME)
lightyellowME = lightyellowME[, 1:2]
#write.csv(yellowME, "GOResults/Highlanders/yellowME.csv")

rownames(lightyellowME) = lightyellowME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(lightyellowME)
lightyellowME0 = musGeneinfo[which(intersection2), ] 
lightyellowMEGO = as.vector(lightyellowME0$mus_genename)

lightyellowMEGO = gost(lightyellowMEGO,
                       organism = "mmusculus",
                       user_threshold = 0.05,
                       custom_bg = background,
                       ordered_query = TRUE,
                       correction_method = "bonferroni")

lightyellowGOResults = data.frame(Term.Name = lightyellowMEGO$result$term_name, 
                                  P.value = lightyellowMEGO$result$p_value, 
                                  Source = lightyellowMEGO$result$source, 
                                  Term.Size = lightyellowMEGO$result$term_size, 
                                  Precision = lightyellowMEGO$result$precision)
write.csv(lightyellowGOResults, 'Labyrinth Zone/GOResults/Lowlanders/lightyellowGOResults.csv')

#royalblue
royalblueME = geneInfoBWLZ %>% filter(moduleColor == "royalblue")
as.data.frame(royalblueME)
royalblueME = royalblueME[, 1:2]
#write.csv(yellowME, "GOResults/Highlanders/yellowME.csv")

rownames(royalblueME) = royalblueME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(royalblueME)
royalblueME0 = musGeneinfo[which(intersection2), ] 
royalblueMEGO = as.vector(royalblueME0$mus_genename)

royalblueMEGO = gost(royalblueMEGO,
                     organism = "mmusculus",
                     user_threshold = 0.05,
                     custom_bg = background,
                     ordered_query = TRUE,
                     correction_method = "bonferroni")

royalblueGOResults = data.frame(Term.Name = royalblueMEGO$result$term_name, 
                                P.value = royalblueMEGO$result$p_value, 
                                Source = royalblueMEGO$result$source, 
                                Term.Size = royalblueMEGO$result$term_size, 
                                Precision = royalblueMEGO$result$precision)
write.csv(royalblueGOResults, 'Labyrinth Zone/GOResults/Lowlanders/royalblueGOResults.csv')

#darkgreen
darkgreenME = geneInfoBWLZ %>% filter(moduleColor == "darkgreen")
as.data.frame(darkgreenME)
darkgreenME = darkgreenME[, 1:2]
#write.csv(yellowME, "GOResults/Highlanders/yellowME.csv")

rownames(darkgreenME) = darkgreenME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(darkgreenME)
darkgreenME0 = musGeneinfo[which(intersection2), ] 
darkgreenMEGO = as.vector(darkgreenME0$mus_genename)

darkgreenMEGO = gost(darkgreenMEGO,
                     organism = "mmusculus",
                     user_threshold = 0.05,
                     custom_bg = background,
                     ordered_query = TRUE,
                     correction_method = "bonferroni")

darkgreenGOResults = data.frame(Term.Name = darkgreenMEGO$result$term_name, 
                                P.value = darkgreenMEGO$result$p_value, 
                                Source = darkgreenMEGO$result$source, 
                                Term.Size = darkgreenMEGO$result$term_size, 
                                Precision = darkgreenMEGO$result$precision)
write.csv(darkgreenGOResults, 'Labyrinth Zone/GOResults/Lowlanders/darkgreenGOResults.csv')

#grey60
grey60ME = geneInfoBWLZ %>% filter(moduleColor == "grey60")
as.data.frame(grey60ME)
grey60ME = grey60ME[, 1:2]
#write.csv(yellowME, "GOResults/Highlanders/yellowME.csv")

rownames(grey60ME) = grey60ME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(grey60ME)
grey60ME0 = musGeneinfo[which(intersection2), ] 
grey60MEGO = as.vector(grey60ME0$mus_genename)

grey60MEGO = gost(grey60MEGO,
                  organism = "mmusculus",
                  user_threshold = 0.05,
                  custom_bg = background,
                  ordered_query = TRUE,
                  correction_method = "bonferroni")

grey60GOResults = data.frame(Term.Name = grey60MEGO$result$term_name, 
                             P.value = grey60MEGO$result$p_value, 
                             Source = grey60MEGO$result$source, 
                             Term.Size = grey60MEGO$result$term_size, 
                             Precision = grey60MEGO$result$precision)
write.csv(grey60GOResults, 'Labyrinth Zone/GOResults/Lowlanders/grey60GOResults.csv')

#darkred
darkredME = geneInfoBWLZ %>% filter(moduleColor == "darkred")
as.data.frame(darkredME)
darkredME = darkredME[, 1:2]
#write.csv(yellowME, "GOResults/Highlanders/yellowME.csv")

rownames(darkredME) = darkredME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(darkredME)
darkredME0 = musGeneinfo[which(intersection2), ] 
darkredMEGO = as.vector(darkredME0$mus_genename)

darkredMEGO = gost(darkredMEGO,
                   organism = "mmusculus",
                   user_threshold = 0.05,
                   custom_bg = background,
                   ordered_query = TRUE,
                   correction_method = "bonferroni")

darkredGOResults = data.frame(Term.Name = darkredMEGO$result$term_name, 
                              P.value = darkredMEGO$result$p_value, 
                              Source = darkredMEGO$result$source, 
                              Term.Size = darkredMEGO$result$term_size, 
                              Precision = darkredMEGO$result$precision)
write.csv(darkredGOResults, 'Labyrinth Zone/GOResults/Lowlanders/darkredGOResults.csv')

#lightgreen
lightgreenME = geneInfoBWLZ %>% filter(moduleColor == "lightgreen")
as.data.frame(lightgreenME)
lightgreenME = lightgreenME[, 1:2]
#write.csv(yellowME, "GOResults/Highlanders/yellowME.csv")

rownames(lightgreenME) = lightgreenME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(lightgreenME)
lightgreenME0 = musGeneinfo[which(intersection2), ] 
lightgreenMEGO = as.vector(lightgreenME0$mus_genename)

lightgreenMEGO = gost(lightgreenMEGO,
                      organism = "mmusculus",
                      user_threshold = 0.05,
                      custom_bg = background,
                      ordered_query = TRUE,
                      correction_method = "bonferroni")

lightgreenGOResults = data.frame(Term.Name = lightgreenMEGO$result$term_name, 
                                 P.value = lightgreenMEGO$result$p_value, 
                                 Source = lightgreenMEGO$result$source, 
                                 Term.Size = lightgreenMEGO$result$term_size, 
                                 Precision = lightgreenMEGO$result$precision)
write.csv(lightgreenGOResults, 'Labyrinth Zone/GOResults/Lowlanders/lightgreenGOResults.csv')

#cyan
cyanME = geneInfoBWLZ %>% filter(moduleColor == "cyan")
as.data.frame(cyanME)
cyanME = cyanME[, 1:2]
#write.csv(yellowME, "GOResults/Highlanders/yellowME.csv")

rownames(cyanME) = cyanME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(cyanME)
cyanME0 = musGeneinfo[which(intersection2), ] 
cyanMEGO = as.vector(cyanME0$mus_genename)

cyanMEGO = gost(cyanMEGO,
                organism = "mmusculus",
                user_threshold = 0.05,
                custom_bg = background,
                ordered_query = TRUE,
                correction_method = "bonferroni")

cyanGOResults = data.frame(Term.Name = cyanMEGO$result$term_name, 
                           P.value = cyanMEGO$result$p_value, 
                           Source = cyanMEGO$result$source, 
                           Term.Size = cyanMEGO$result$term_size, 
                           Precision = cyanMEGO$result$precision)
write.csv(cyanGOResults, 'Labyrinth Zone/GOResults/Lowlanders/cyanGOResults.csv')

#midnightblue
midnightblueME = geneInfoBWLZ %>% filter(moduleColor == "midnightblue")
as.data.frame(midnightblueME)
midnightblueME = midnightblueME[, 1:2]
#write.csv(yellowME, "GOResults/Highlanders/yellowME.csv")

rownames(midnightblueME) = midnightblueME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(midnightblueME)
midnightblueME0 = musGeneinfo[which(intersection2), ] 
midnightblueMEGO = as.vector(midnightblueME0$mus_genename)

midnightblueMEGO = gost(midnightblueMEGO,
                        organism = "mmusculus",
                        user_threshold = 0.05,
                        custom_bg = background,
                        ordered_query = TRUE,
                        correction_method = "bonferroni")

midnightblueGOResults = data.frame(Term.Name = midnightblueMEGO$result$term_name, 
                                   P.value = midnightblueMEGO$result$p_value, 
                                   Source = midnightblueMEGO$result$source, 
                                   Term.Size = midnightblueMEGO$result$term_size, 
                                   Precision = midnightblueMEGO$result$precision)
write.csv(midnightblueGOResults, 'Labyrinth Zone/GOResults/Lowlanders/midnightblueGOResults.csv')

#tan
tanME = geneInfoBWLZ %>% filter(moduleColor == "tan")
as.data.frame(tanME)
tanME = tanME[, 1:2]
#write.csv(yellowME, "GOResults/Highlanders/yellowME.csv")

rownames(tanME) = tanME$Gene
intersection2 = musGeneinfo$geneID  %in% rownames(tanME)
tanME0 = musGeneinfo[which(intersection2), ] 
tanMEGO = as.vector(tanME0$mus_genename)

tanMEGO = gost(tanMEGO,
               organism = "mmusculus",
               user_threshold = 0.05,
               custom_bg = background,
               ordered_query = TRUE,
               correction_method = "bonferroni")

tanGOResults = data.frame(Term.Name = tanMEGO$result$term_name, 
                          P.value = tanMEGO$result$p_value, 
                          Source = tanMEGO$result$source, 
                          Term.Size = tanMEGO$result$term_size, 
                          Precision = tanMEGO$result$precision)
write.csv(tanGOResults, 'Labyrinth Zone/GOResults/Lowlanders/tanGOResults.csv')


# Do preserved or non-preserved lowlander modules have targets of selection in them?? ####
targets = read.csv("DE_LZ_allGenes_wPvals_revised.csv")
targets = data.frame(Genes = targets$geneID, musGeneName = targets$mus_genename, targets =  targets$PBSRDA_outlier)

# Preserved modules
#darkgreen
darkgreentargets = geneInfoBWLZ %>% filter(geneInfoBWLZ$moduleColor == 'darkgreen')
intersection = targets$Genes %in% darkgreentargets$Gene
darkgreentargets = targets[which(intersection), ]
darkgreentargets = darkgreentargets %>% filter(darkgreentargets$targets == "YES")
length(darkgreentargets$targets) # 1 target gene
darkgreentargets$musGeneName

#green
greentargets = geneInfoBWLZ %>% filter(geneInfoBWLZ$moduleColor == 'green')
intersection = targets$Genes %in% greentargets$Gene
greentargets = targets[which(intersection), ]
greentargets = greentargets %>% filter(greentargets$targets == "YES")
length(greentargets$targets) # 19 target gene
greentargets$musGeneName

#greenyellow
greenyellowtargets = geneInfoBWLZ %>% filter(geneInfoBWLZ$moduleColor == 'greenyellow')
intersection = targets$Genes %in% greenyellowtargets$Gene
greenyellowtargets = targets[which(intersection), ]
greenyellowtargets = greenyellowtargets %>% filter(greenyellowtargets$targets == "YES")
length(greenyellowtargets$targets) # 9 target gene
greenyellowtargets$musGeneName

#midnightblue
midnightbluetargets = geneInfoBWLZ %>% filter(geneInfoBWLZ$moduleColor == 'midnightblue')
intersection = targets$Genes %in% midnightbluetargets$Gene
midnightbluetargets = targets[which(intersection), ]
midnightbluetargets = midnightbluetargets %>% filter(midnightbluetargets$targets == "YES")
length(midnightbluetargets$targets) # 12 target gene
midnightbluetargets$musGeneName

#red
redtargets = geneInfoBWLZ %>% filter(geneInfoBWLZ$moduleColor == 'red')
intersection = targets$Genes %in% redtargets$Gene
redtargets = targets[which(intersection), ]
redtargets = redtargets %>% filter(redtargets$targets == "YES")
length(redtargets$targets) # 34 target gene
redtargets$musGeneName

#Weakly Preserved
#lightyellow
lightyellowtargets = geneInfoBWLZ %>% filter(geneInfoBWLZ$moduleColor == 'lightyellow')
intersection = targets$Genes %in% lightyellowtargets$Gene
lightyellowtargets = targets[which(intersection), ]
lightyellowtargets = lightyellowtargets %>% filter(lightyellowtargets$targets == "YES")
length(lightyellowtargets$targets) # 34 target gene
lightyellowtargets$musGeneName


# ARE GENES IN 2x2 table present in any significant module-trait correlations ####
# , or preserved/non-preserved networks? 

#Junctional Zone
JZsummary = readxl::read_xlsx("summaryJZ.xlsx")
JZsummary = as.data.frame(JZsummary)
colnames(JZsummary)[colnames(JZsummary) == ("...1")] = "Gene"
JZsummary = JZsummary %>% select(-ID)

#Labyrinth Zone
LZsummary = readxl::read_xlsx("summaryLZ.xlsx")
LZsummary = as.data.frame(LZsummary)
colnames(LZsummary)[colnames(LZsummary) == ("...1")] = "Gene"
  #Add gene names to data frames

# Filter ff, fp, pf, pp genes in JZ
ffJZ = JZsummary %>% filter(lowland == 'FIXED' & highland == 'FIXED')
  intersection = ffJZ$Gene %in% musGeneinfo$geneID
  ffJZ = ffJZ[which(intersection), ] 

fpJZ = JZsummary %>% filter(lowland == "FIXED" & highland == "PLASTIC")
  intersection = fpJZ$Gene %in% musGeneinfo$geneID
  fpJZ = fpJZ[which(intersection), ] 

pfJZ = JZsummary %>% filter(lowland == "PLASTIC" & highland == "FIXED")
  intersection = pfJZ$Gene %in% musGeneinfo$geneID
  pfJZ = pfJZ[which(intersection), ] 

ppJZ = JZsummary %>% filter(lowland == "PLASTIC" & highland == "PLASTIC")
  intersection = ppJZ$Gene %in% musGeneinfo$geneID
  ppJZ = ppJZ[which(intersection), ] 
  
# Filter ff, fp, pf, pp genes in LZ
ffLZ = LZsummary %>% filter(lowland == 'FIXED' & highland == 'FIXED')
  intersection = ffLZ$Gene %in% musGeneinfo$geneID
  ffLZ = ffLZ[which(intersection), ] 

fpLZ = LZsummary %>% filter(lowland == "FIXED" & highland == "PLASTIC")
  intersection = fpLZ$Gene %in% musGeneinfo$geneID
  fpLZ = fpLZ[which(intersection), ] 

pfLZ = LZsummary %>% filter(lowland == "PLASTIC" & highland == "FIXED")
  intersection = pfLZ$Gene %in% musGeneinfo$geneID
  pfLZ = pfLZ[which(intersection), ] 

ppLZ = LZsummary %>% filter(lowland == "PLASTIC" & highland == "PLASTIC")
  intersection = ppLZ$Gene %in% musGeneinfo$geneID
  ppLZ = ppLZ[which(intersection), ] 


  #LZ Highland and Lowland
  #ff
  intersection = geneInfoMELZ$Gene %in% ffLZ$Gene 
  ffGeneModMELZ = geneInfoMELZ[which(intersection), ]
  ffGeneModMELZ = ffGeneModMELZ %>% filter(!moduleColor %in% c("grey"))
  
  intersection = geneInfoBWLZ$Gene %in% ffLZ$Gene 
  ffGeneModBWLZ = geneInfoBWLZ[which(intersection), ]
  ffGeneModBWLZ = ffGeneModBWLZ %>% filter(!moduleColor %in% c("grey"))
  
  #fp
  intersection = geneInfoMELZ$Gene %in% fpLZ$Gene 
  fpGeneModMELZ = geneInfoMELZ[which(intersection), ]
  fpGeneModMELZ = fpGeneModMELZ %>% filter(!moduleColor %in% c("grey"))
  
  intersection = geneInfoBWLZ$Gene %in% fpLZ$Gene 
  fpGeneModBWLZ = geneInfoBWLZ[which(intersection), ]
  fpGeneModBWLZ = fpGeneModBWLZ %>% filter(!moduleColor %in% c("grey"))
  
  #pf
  intersection = geneInfoMELZ$Gene %in% pfLZ$Gene 
  pfGeneModMELZ = geneInfoMELZ[which(intersection), ]
  pfGeneModMELZ = pfGeneModMELZ %>% filter(!moduleColor %in% c("grey"))
  
  intersection = geneInfoBWJZ$Gene %in% pfLZ$Gene 
  pfGeneModBWLZ = geneInfoBWJZ[which(intersection), ]
  pfGeneModBWLZ = pfGeneModBWLZ %>% filter(!moduleColor %in% c("grey"))
  
  #pp
  intersection = geneInfoMELZ$Gene %in% ppLZ$Gene 
  ppGeneModMELZ = geneInfoMELZ[which(intersection), ]
  ppGeneModMELZ = ppGeneModMELZ %>% filter(!moduleColor %in% c("grey"))
  
  intersection = geneInfoBWLZ$Gene %in% ppLZ$Gene 
  ppGeneModBWLZ = geneInfoBWLZ[which(intersection), ]
  ppGeneModBWLZ = ppGeneModBWLZ %>% filter(!moduleColor %in% c("grey"))
  
  ## Look at 2x2 genes in lowlanders and what modules they show up in  ####
  fpcharacter_counts = table(fpGeneModBWLZ$moduleColor)
  print(fpcharacter_counts)
  
  pfcharacter_counts = table(pfGeneModBWLZ$moduleColor)
  print(pfcharacter_counts)
  
  ppcharacter_counts = table(ppGeneModBWLZ$moduleColor)
  print(ppcharacter_counts)
  
  ffcharacter_counts = table(ffGeneModBWLZ$moduleColor)
  print(ffcharacter_counts)


## GO ANALYSIS OF THESE GENES in LZ and JZ ####
## choose the module as the background gene set. 
musGeneinfo = read.csv("DE_LZ_allGenes_wPvals_revised.csv")

library(gprofiler2)
background = as.vector(musGeneinfo$mus_genename)

# Plastic/Plastic Turquoise genes
ppTurqLZ = ppGeneModMELZ %>% filter(moduleColor == "turquoise")

rownames(ppTurqLZ) = ppTurqLZ$Gene
intersection = musGeneinfo$geneID  %in% rownames(ppTurqLZ)
ppTurqLZ0 = musGeneinfo[which(intersection), ] 
ppTurqLZGO = as.vector(ppTurqLZ0$mus_genename)

ppLZturqGO = gost(ppTurqLZGO,
              organism = "mmusculus",
              user_threshold = 0.05,
              custom_bg = background,
              ordered_query = TRUE,
              correction_method = "bonferroni")

ppLZturqGOResults = data.frame(Term.Name = ppLZturqGO$result$term_name, 
                           P.value = ppLZturqGO$result$p_value, 
                           Source = ppLZturqGO$result$source, 
                           Term.Size = ppLZturqGO$result$term_size, 
                           Precision = ppLZturqGO$result$precision)
write.csv(ppLZturqGOResults, 'GOResults/Highlanders/ppLZturqGOResults.csv')

# Plastic/Plastic Blue genes
ppBlueLZ = ppGeneModMELZ %>% filter(moduleColor == "blue")

rownames(ppBlueLZ) = ppBlueLZ$Gene
intersection = musGeneinfo$geneID  %in% rownames(ppBlueLZ)
ppBlueLZ0 = musGeneinfo[which(intersection), ] 
ppBlueLZGO = as.vector(ppBlueLZ0$mus_genename)

ppLZblueGO = gost(ppBlueLZGO,
                  organism = "mmusculus",
                  user_threshold = 0.05,
                  custom_bg = background,
                  ordered_query = TRUE,
                  correction_method = "bonferroni")

ppLZblueGOResults = data.frame(Term.Name = ppLZblueGO$result$term_name, 
                               P.value = ppLZblueGO$result$p_value, 
                               Source = ppLZblueGO$result$source, 
                               Term.Size = ppLZblueGO$result$term_size, 
                               Precision = ppLZblueGO$result$precision)
write.csv(ppLZblueGOResults, 'GOResults/Highlanders/ppLZblueGOResults.csv')

# Plastic/Plastic brown genes
ppBrownLZ = ppGeneModMELZ %>% filter(moduleColor == "brown")

rownames(ppBrownLZ) = ppBrownLZ$Gene
intersection = musGeneinfo$geneID  %in% rownames(ppBrownLZ)
ppBrownLZ0 = musGeneinfo[which(intersection), ] 
ppBrownLZGO = as.vector(ppBrownLZ0$mus_genename)

ppLZbrownGO = gost(ppBrownLZGO,
                  organism = "mmusculus",
                  user_threshold = 0.05,
                  custom_bg = background,
                  ordered_query = TRUE,
                  correction_method = "bonferroni")

ppLZbrownGOResults = data.frame(Term.Name = ppLZbrownGO$result$term_name, 
                               P.value = ppLZbrownGO$result$p_value, 
                               Source = ppLZbrownGO$result$source, 
                               Term.Size = ppLZbrownGO$result$term_size, 
                               Precision = ppLZbrownGO$result$precision)
write.csv(ppLZbrownGOResults, 'GOResults/Highlanders/ppLZbrownGOResults.csv')

# Plastic/Plastic yellow genes
ppYellowLZ = ppGeneModMELZ %>% filter(moduleColor == "yellow")

rownames(ppYellowLZ) = ppYellowLZ$Gene
intersection = musGeneinfo$geneID  %in% rownames(ppYellowLZ)
ppYellowLZ0 = musGeneinfo[which(intersection), ] 
ppYellowLZGO = as.vector(ppYellowLZ0$mus_genename)

ppLZyellowGO = gost(ppYellowLZGO,
                   organism = "mmusculus",
                   user_threshold = 0.05,
                   custom_bg = background,
                   ordered_query = TRUE,
                   correction_method = "bonferroni")

ppLZyellowGOResults = data.frame(Term.Name = ppLZyellowGO$result$term_name, 
                                P.value = ppLZyellowGO$result$p_value, 
                                Source = ppLZyellowGO$result$source, 
                                Term.Size = ppLZyellowGO$result$term_size, 
                                Precision = ppLZyellowGO$result$precision)
write.csv(ppLZyellowGOResults, 'GOResults/Highlanders/ppLZyellowGOResults.csv')

# Plastic/Plastic red genes
ppRedLZ = ppGeneModMELZ %>% filter(moduleColor == "red")

rownames(ppRedLZ) = ppRedLZ$Gene
intersection = musGeneinfo$geneID  %in% rownames(ppRedLZ)
ppRedLZ0 = musGeneinfo[which(intersection), ] 
ppRedLZGO = as.vector(ppRedLZ0$mus_genename)

ppLZredGO = gost(ppRedLZGO,
                    organism = "mmusculus",
                    user_threshold = 0.05,
                    custom_bg = background,
                    ordered_query = TRUE,
                    correction_method = "bonferroni")

ppLZredGOResults = data.frame(Term.Name = ppLZredGO$result$term_name, 
                                 P.value = ppLZredGO$result$p_value, 
                                 Source = ppLZredGO$result$source, 
                                 Term.Size = ppLZredGO$result$term_size, 
                                 Precision = ppLZredGO$result$precision)
write.csv(ppLZredGOResults, 'GOResults/Highlanders/ppLZredGOResults.csv')
##
