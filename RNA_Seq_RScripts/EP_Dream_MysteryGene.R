
# Packages
library('BiocParallel')
library('dplyr')
library('ggplot2')
library('lme4')
library('readxl')
library('variancePartition')
library('edgeR')
library('Matrix')

# Combined Analysis  ------------------
Sample_Info <- read_xlsx("RNA_Seq_RawData/MetaData_EP.xlsx")
Sample_Info = as.data.frame(Sample_Info)
rownames(Sample_Info) = Sample_Info$Sample_ID
Sample_Info = Sample_Info %>% filter(Strain == "BW" | Strain == "ME")

# Read in Files + QC
Pman_rawreads = read_xlsx("RNA_Seq_RawData/EP_Pman_ExtMMFrac_readcounts.xlsx")
Pman_rawreads = as.data.frame(Pman_rawreads)
Pman_rawreads = Pman_rawreads %>%
  filter(!is.na(Geneid))
Pman_rawreads = `row.names<-`(Pman_rawreads, Pman_rawreads$Geneid)
Pman_rawreads <- Pman_rawreads[,-c(1:6)]

# Check read table vs sample info
Check = Sample_Info$Seq_Name
colnames(Pman_rawreads) == Check

colnames(Pman_rawreads) = rownames(Sample_Info)

Pman_readcounts <- as.matrix(Pman_rawreads)
dPman_0 <- DGEList(Pman_readcounts)

dPman_0 <- calcNormFactors(dPman_0)
dim(dPman_0)
keep <- rowSums(cpm(dPman_0) > 0.5 ) >= 60
dPman <- dPman_0[keep,]
dim(dPman)
plotMDS(dPman, col = as.numeric(Sample_Info$Strain), labels = Sample_Info$Strain)

rownames(Sample_Info) == colnames(dPman)

# Interaction
param = SnowParam(8, "SOCK", progressbar=TRUE)
form <- ~ Strain*O2 + (1|Mom)
vobjDream = voomWithDreamWeights(dPman, form, Sample_Info, BPPARAM=param, plot = T)
fitmm = dream( vobjDream, form, Sample_Info, ddf = "Kenward-Roger")
fitmm = eBayes(fitmm)

DE_strain <- topTable( fitmm, coef='StrainME', sort.by = "P", n = Inf, )
DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)
DE_ixn <- topTable( fitmm, coef='StrainME:O22H', sort.by = "P", n = Inf)

length(DE_ixn$logFC[which(DE_ixn$adj.P.Val < 0.05)])
length(DE_strain$logFC[which(DE_strain$adj.P.Val < 0.05)])
length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])

# Summarize Reads
library(emmeans)

Strain <- Sample_Info$Strain
O2 <- Sample_Info$O2
summary <- data.frame()
all_genes <- row.names(DE_ixn)

for (p in all_genes) {
  
  gene_id <- p
  test_data <- vobjDream$E[gene_id,]
  test_model <- lmer(test_data ~ Strain*O2 + (1|(Sample_Info$Mom)))
  anova(test_model)
  output <- summary(pairs(emmeans(test_model, ~ Strain*O2), adjust = "BH"))
  output_line <- output$p.value
  output_line[7] <- gene_id
  
  summary <- rbind(summary, output_line)
}

output_colnames <- output$contrast
colnames(summary) <- output_colnames

corrCounts <- t(vobjDream$E)
corrCounts <- as.data.frame(corrCounts)
corrCounts$Treatment <- Sample_Info$Group
corrCounts$ID <- row.names(corrCounts)
corrCounts$Check <- Sample_Info$Sample_ID
corrCounts$Check == row.names(corrCounts)

MeanCounts <- corrCounts %>%
  group_by(corrCounts$Treatment) %>%
  summarise_at(colnames(corrCounts), funs(mean(., na.rm=TRUE)))
meanCounts <- t(MeanCounts)
meanCounts <- as.data.frame(meanCounts)
colnames(meanCounts) <- meanCounts[1,]
meanCounts <- meanCounts[-1,]
meanCounts$ID <- row.names(meanCounts)

final <- merge(meanCounts, summary, by.x = "ID", by.y = c(7))

# Write Out
write.csv(final, "RNA_Seq_Output/Test_Dream_RawFiles/EP_dreamISO_Contrasts.csv")
write.csv(DE_strain, "RNA_Seq_Output/Test_Dream_RawFiles/EP_dreamISO_strainDE.csv")
write.csv(DE_o2, "RNA_Seq_Output/Test_Dream_RawFiles/EP_dreamISO_o2DE.csv")
write.csv(DE_ixn, "RNA_Seq_Output/Test_Dream_RawFiles/EP_dreamISO_ixnDE.csv")


# BW Only ------------------
Sample_Info <- read_xlsx("RNA_Seq_RawData/MetaData_EP.xlsx")
Sample_Info = as.data.frame(Sample_Info)
rownames(Sample_Info) = Sample_Info$Sample_ID
Sample_Info = Sample_Info %>% filter(Strain == "BW")
BW = Sample_Info$Seq_Name

# Read in Files + QC
Pman_rawreads = read_xlsx("RNA_Seq_RawData/EP_Pman_ExtMMFrac_readcounts.xlsx")
Pman_rawreads = as.data.frame(Pman_rawreads)
Pman_rawreads = Pman_rawreads %>%
  filter(!is.na(Geneid))
Pman_rawreads = `row.names<-`(Pman_rawreads, Pman_rawreads$Geneid)
Pman_rawreads <- Pman_rawreads[,-c(1:6)]
Pman_rawreads <- Pman_rawreads[,names(Pman_rawreads) %in% BW]

# Check read table vs sample info
Check = Sample_Info$Seq_Name
colnames(Pman_rawreads) == Check

colnames(Pman_rawreads) = rownames(Sample_Info)

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
write.csv(DE_o2, "RNA_Seq_Output/Test_Dream_RawFiles/EP_BW_dreamISO_o2DE.csv")


# ME Only ------------------
Sample_Info <- read_xlsx("RNA_Seq_RawData/MetaData_EP.xlsx")
Sample_Info = as.data.frame(Sample_Info)
rownames(Sample_Info) = Sample_Info$Sample_ID
Sample_Info = Sample_Info %>% filter(Strain == "ME")
ME = Sample_Info$Seq_Name

# Read in Files + QC
Pman_rawreads = read_xlsx("RNA_Seq_RawData/EP_Pman_ExtMMFrac_readcounts.xlsx")
Pman_rawreads = as.data.frame(Pman_rawreads)
Pman_rawreads = Pman_rawreads %>%
  filter(!is.na(Geneid))
Pman_rawreads = `row.names<-`(Pman_rawreads, Pman_rawreads$Geneid)
Pman_rawreads <- Pman_rawreads[,-c(1:6)]
Pman_rawreads <- Pman_rawreads[,names(Pman_rawreads) %in% ME]

# Check read table vs sample info
Check = Sample_Info$Seq_Name
colnames(Pman_rawreads) == Check

colnames(Pman_rawreads) = rownames(Sample_Info)

Pman_readcounts <- as.matrix(Pman_rawreads)
dPman_0 <- DGEList(Pman_readcounts)

dPman_0 <- calcNormFactors(dPman_0)
dim(dPman_0)

dPman_ME <- dPman_0[keep,]
dim(dPman_ME)
plotMDS(dPman_ME, col = as.numeric(Sample_Info$Strain), labels = Sample_Info$Strain)

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
write.csv(DE_o2, "RNA_Seq_Output/Test_Dream_RawFiles/EP_ME_dreamISO_o2DE.csv")
