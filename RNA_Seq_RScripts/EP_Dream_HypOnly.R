
# Packages
library(BiocParallel)
library(dplyr)
library(ggplot2)
library(lme4)
library(readxl)
library(variancePartition)
library(edgeR)
library(Matrix)

# Hypoxia Sensitive genes only

EP_Summary = read_xlsx("RNA_Seq_Output/EP_ISO_Ortho_Summary.xlsx")

gene_class <- EP_Summary %>%
  filter(O2_SIG == "SIG" | IXN_SIG == "SIG") %>%
  mutate(Category = case_when(
    IXN_SIG == "SIG" ~ "IXN_SIG",
    O2_SIG == "SIG" ~ "O2_SIG")) %>%
  distinct(Pman_GeneID, .keep_all = TRUE)

Hypoxia_genes= gene_class %>%
  filter(Category == "O2_SIG") %>%
  pull(Pman_GeneID)

# Combined Analysis  ------------------
Sample_Info <- read_xlsx("RNA_Seq_RawData/MetaData_EP.xlsx")
Sample_Info = as.data.frame(Sample_Info)
rownames(Sample_Info) = Sample_Info$Sample_ID
Sample_Info = Sample_Info %>% filter(Strain == "BW" | Strain == "ME")

# Read in Files + QC
Pman_rawreads = read_xlsx("RNA_Seq_RawData/EP_Pman_ExtMMFrac_readcounts_Exon.xlsx")
Pman_rawreads = as.data.frame(Pman_rawreads)
Pman_rawreads = `row.names<-`(Pman_rawreads, Pman_rawreads$Geneid)
Pman_rawreads <- Pman_rawreads[,-c(1:6)]


# Check read table vs sample info
Check = Sample_Info$Seq_Name
colnames(Pman_rawreads) == Check
colnames(Pman_rawreads) = rownames(Sample_Info)

# Hypoxia and interaction genes only!
Pman_rawreads = Pman_rawreads[rownames(Pman_rawreads) %in% Hypoxia_genes, ]

Pman_readcounts <- as.matrix(Pman_rawreads)
dPman_0 <- DGEList(Pman_readcounts)
dPman_0 <- calcNormFactors(dPman_0)
dim(dPman_0)
keep <- rowSums(cpm(dPman_0) > 0.5 ) >= 60
dPman <- dPman_0[keep,]
dim(dPman)
plotMDS(dPman, col = as.numeric(Sample_Info$Strain), labels = Sample_Info$Strain)

rownames(Sample_Info) == colnames(dPman)

# # Interaction
# param = SnowParam(8, "SOCK", progressbar=TRUE)
# form <- ~ Strain*O2 + (1|Mom)
# vobjDream = voomWithDreamWeights(dPman, form, Sample_Info, BPPARAM=param, plot = T)
# fitmm = dream( vobjDream, form, Sample_Info, ddf = "Kenward-Roger")
# fitmm = eBayes(fitmm)
# 
# DE_strain <- topTable( fitmm, coef='StrainME', sort.by = "P", n = Inf, )
# DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)
# DE_ixn <- topTable( fitmm, coef='StrainME:O22H', sort.by = "P", n = Inf)
# 
# # Write Out
# write.csv(DE_strain, "RNA_Seq_Output/Dream_RawFiles/EP_HypOnly_dreamISO_strainDE.csv")
# write.csv(DE_o2, "RNA_Seq_Output/Dream_RawFiles/EP_HypOnly_dreamISO_o2DE.csv")
# write.csv(DE_ixn, "RNA_Seq_Output/Dream_RawFiles/EP_HypOnly_dreamISO_ixnDE.csv")


# BW Only ------------------
Sample_Info <- read_xlsx("RNA_Seq_RawData/MetaData_EP.xlsx")
Sample_Info = as.data.frame(Sample_Info)
rownames(Sample_Info) = Sample_Info$Sample_ID
Sample_Info = Sample_Info %>% filter(Strain == "BW")
BW = Sample_Info$Seq_Name

# Read in Files + QC
Pman_rawreads = read_xlsx("RNA_Seq_RawData/EP_Pman_ExtMMFrac_readcounts_Exon.xlsx")
Pman_rawreads = as.data.frame(Pman_rawreads)
Pman_rawreads = `row.names<-`(Pman_rawreads, Pman_rawreads$Geneid)
Pman_rawreads <- Pman_rawreads[,-c(1:6)]
Pman_rawreads <- Pman_rawreads[,names(Pman_rawreads) %in% BW]

# Check read table vs sample info
Check = Sample_Info$Seq_Name
colnames(Pman_rawreads) == Check

colnames(Pman_rawreads) = rownames(Sample_Info)

# Hypoxia and interaction genes only!
Pman_rawreads = Pman_rawreads[rownames(Pman_rawreads) %in% Hypoxia_genes, ]

Pman_readcounts <- as.matrix(Pman_rawreads)
dPman_0 <- DGEList(Pman_readcounts)

dPman_0 <- calcNormFactors(dPman_0)
dim(dPman_0)

dPman_BW <- dPman_0[keep,]
dim(dPman_BW)
plotMDS(dPman_BW, col = as.numeric(Sample_Info$Strain), labels = Sample_Info$Strain)

rownames(Sample_Info) == colnames(dPman_BW)

# Interaction  
param = SnowParam(8, "SOCK", progressbar=TRUE)
form <- ~ O2 + (1|Mom)
vobjDream = voomWithDreamWeights(dPman_BW, form, Sample_Info, BPPARAM=param, plot = T)
fitmm = dream(vobjDream, form, Sample_Info, ddf = "Kenward-Roger")
fitmm = eBayes(fitmm)

DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)

length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])

# Write Out
write.csv(DE_o2, "RNA_Seq_Output/Dream_RawFiles/EP_BW_HypOnly_dreamISO_o2DE.csv")


# ME Only ------------------
Sample_Info <- read_xlsx("RNA_Seq_RawData/MetaData_EP.xlsx")
Sample_Info = as.data.frame(Sample_Info)
rownames(Sample_Info) = Sample_Info$Sample_ID
Sample_Info = Sample_Info %>% filter(Strain == "ME")
ME = Sample_Info$Seq_Name

# Read in Files + QC
Pman_rawreads = read_xlsx("RNA_Seq_RawData/EP_Pman_ExtMMFrac_readcounts_Exon.xlsx")
Pman_rawreads = as.data.frame(Pman_rawreads)
Pman_rawreads = `row.names<-`(Pman_rawreads, Pman_rawreads$Geneid)
Pman_rawreads <- Pman_rawreads[,-c(1:6)]
Pman_rawreads <- Pman_rawreads[,names(Pman_rawreads) %in% ME]

# Check read table vs sample info
Check = Sample_Info$Seq_Name
colnames(Pman_rawreads) == Check

colnames(Pman_rawreads) = rownames(Sample_Info)

# Hypoxia and interaction genes only!
Pman_rawreads = Pman_rawreads[rownames(Pman_rawreads) %in% Hypoxia_genes, ]

Pman_readcounts <- as.matrix(Pman_rawreads)
dPman_0 <- DGEList(Pman_readcounts)

dPman_0 <- calcNormFactors(dPman_0)
dim(dPman_0)

dPman_ME <- dPman_0[keep,]
dim(dPman_ME)
plotMDS(dPman_ME, col = as.numeric(Sample_Info$Strain), labels = Sample_Info$Sample_ID)

rownames(Sample_Info) == colnames(dPman_ME)

# Interaction  
param = SnowParam(8, "SOCK", progressbar=TRUE)
form <- ~ O2 + (1|Mom)
vobjDream = voomWithDreamWeights(dPman_ME, form, Sample_Info, BPPARAM=param, plot = T)
fitmm = dream(vobjDream, form, Sample_Info, ddf = "Kenward-Roger")
fitmm = eBayes(fitmm)

DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)

length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])

# Write Out
write.csv(DE_o2, "RNA_Seq_Output/Dream_RawFiles/EP_ME_HypOnly_dreamISO_o2DE.csv")
