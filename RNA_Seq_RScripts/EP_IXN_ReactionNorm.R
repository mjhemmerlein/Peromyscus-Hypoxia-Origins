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
library(reshape2)
library(viridis)

# Early pregnancy raw data
EP_Sample_Info <- read_xlsx("RNA_Seq_RawData/MetaData_EP.xlsx")
EP_Sample_Info <- as.data.frame(EP_Sample_Info)
rownames(EP_Sample_Info) <- EP_Sample_Info$Sample_ID
EP_Sample_Info <- EP_Sample_Info %>% filter(Strain == "BW" | Strain == "ME")
EP_Sample_Info$Group <- as.factor(EP_Sample_Info$Group)
EP_Sample_Info$Group_Timepoint = paste0("EP_", EP_Sample_Info$Group)

Pman_rawreads_EP <- read_xlsx("RNA_Seq_RawData/EP_Pman_ExtMMFrac_readcounts_Exon.xlsx")
Pman_rawreads_EP <- as.data.frame(Pman_rawreads_EP)
Pman_rawreads_EP <- Pman_rawreads_EP %>% filter(!is.na(Geneid))
row.names(Pman_rawreads_EP) <- Pman_rawreads_EP$Geneid
Pman_rawreads_EP <- Pman_rawreads_EP[,-c(1:6)]

# Check and match columns
Check_EP <- EP_Sample_Info$Seq_Name
colnames(Pman_rawreads_EP) == Check_EP
colnames(Pman_rawreads_EP) <- rownames(EP_Sample_Info)

# Filter data
Pman_readcounts_EP <- as.matrix(Pman_rawreads_EP)
dPman_0_EP <- DGEList(Pman_readcounts_EP)
dPman_0_EP <- calcNormFactors(dPman_0_EP)
keep_EP <- rowSums(cpm(dPman_0_EP) > 0.5) >= 60
dPman_EP <- dPman_0_EP[keep_EP,]
dim(dPman_EP)

# Log-CPM values
logcpm <- cpm(dPman_EP, log = TRUE)

EP_Summary = read_xlsx("RNA_Seq_Output/EP_ISO_Ortho_Summary.xlsx")

IXN = EP_Summary %>%
  filter(IXN_SIG == "SIG")

ixn_genes <- IXN$Pman_GeneID

logcpm_subset = logcpm[rownames(logcpm) %in% ixn_genes, , drop = FALSE]

logcpm_df = as.data.frame(logcpm_subset)

logcpm_df_long = logcpm_df %>%
  rownames_to_column("Gene_ID") %>%
  pivot_longer(-Gene_ID, names_to = "Sample_ID", values_to = "Expression")

logcpm_df_long = logcpm_df_long %>%
  left_join(EP_Sample_Info, by = "Sample_ID")

logcpm_df_long$Strain = as.factor(logcpm_df_long$Strain)
logcpm_df_long$O2 = as.factor(logcpm_df_long$O2)

logcpm_df_long = logcpm_df_long %>%
  mutate(Strain = recode(Strain,
                         "BW" = "Lowland",
                         "ME" = "Highland"),
         O2 = recode(O2,
                     "1N" = "Normoxia",
                     "2H" = "Hypoxia"))


# calculate mean expression per gene × group
mean_expr_df <- logcpm_df_long %>%
  group_by(Gene_ID, Group) %>%
  summarise(mean_expr = mean(Expression, na.rm = TRUE), .groups = "drop")

logcpm_df_long <- logcpm_df_long %>%
  left_join(mean_expr_df, by = c("Gene_ID", "Group"))

ggplot(logcpm_df_long, aes(x = O2, y = mean_expr, group = Strain, color = Strain)) +
  scale_color_manual(values = c("#F4C552", "#247BB1")) +
  geom_line(alpha = 0.5) +
  geom_point(size = 2) +
  theme_minimal() +
  facet_wrap(~ Gene_ID, scales = "free_y", ncol = 4) +
  labs(
    x = "",
    y = "Mean Expression",
    title = "Reaction Norm Plot"
  )

#ggsave("Plots/Faceted_Gene_Plots/IXN_Reaction_Norm.pdf", width = 11, height = 10, units = "in", dpi = 300)

