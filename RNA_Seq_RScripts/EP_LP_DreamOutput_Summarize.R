# Summarize raw Dream Output file

# Libraries -----------
library(tidyverse)
library(readr)
library(readxl)

# PBS
PBS = read_xlsx("RNA_Seq_RawData/PBS_RDA_outliers.xlsx")
PBS = PBS %>%
  rename(Gene_ID = gene_pman) %>%
  select(Gene_ID) %>%
  mutate(PBS = "TRUE")

Apriori = read_xlsx("RNA_Seq_RawData/Pman_Apriori.xlsx")
Apriori = Apriori %>%
  rename(Gene_ID = geneName) %>% 
  select(Gene_ID) %>%
  mutate(APRI = "TRUE")

# IsoQuant Early Pregnancy -----------
# Filtered at 0.5 CPM

EP_ISO_Strain = read_csv("RNA_Seq_Output/Dream_RawFiles/EP_dreamISO_strainDE.csv")
colnames(EP_ISO_Strain)[1] = "Gene_ID"

EP_ISO_Strain = EP_ISO_Strain %>%
  select(Gene_ID, logFC, P.Value, adj.P.Val) %>%
  mutate(strain_SIG = ifelse(adj.P.Val < 0.05, "SIG", "NA")) %>%
  rename(strain_logFC = logFC) %>% 
  rename(strain_P.Val = P.Value) %>% 
  rename(strain_adj_P.Val = adj.P.Val)

EP_ISO_O2 = read_csv("RNA_Seq_Output/Dream_RawFiles/EP_dreamISO_o2DE.csv")
colnames(EP_ISO_O2)[1] = "Gene_ID"

EP_ISO_O2 = EP_ISO_O2 %>%
  select(Gene_ID, logFC, P.Value, adj.P.Val) %>%
  mutate(O2_SIG = ifelse(adj.P.Val < 0.05, "SIG", "NA")) %>% 
  rename(O2_logFC = logFC) %>% 
  rename(O2_P.Val = P.Value) %>% 
  rename(O2_adj_P.Val = adj.P.Val)

EP_ISO_IXN = read_csv("RNA_Seq_Output/Dream_RawFiles/EP_dreamISO_ixnDE.csv")
colnames(EP_ISO_IXN)[1] = "Gene_ID"

EP_ISO_IXN = EP_ISO_IXN %>%
  select(Gene_ID, logFC, P.Value, adj.P.Val) %>%
  mutate(IXN_SIG = ifelse(adj.P.Val < 0.05, "SIG", "NA")) %>% 
  rename(IXN_logFC = logFC) %>% 
  rename(IXN_P.Val = P.Value) %>% 
  rename(IXN_adj_P.Val = adj.P.Val)

EP_ISO_BW = read_csv("RNA_Seq_Output/Dream_RawFiles/EP_BW_dreamISO_o2DE.csv")
colnames(EP_ISO_BW)[1] = "Gene_ID"

EP_ISO_BW = EP_ISO_BW %>%
  select(Gene_ID, logFC, P.Value, adj.P.Val) %>%
  mutate(BW_O2_SIG = ifelse(adj.P.Val < 0.05, "SIG", "NA")) %>% 
  rename(BW_O2_logFC = logFC) %>% 
  rename(BW_O2_P.Val = P.Value) %>% 
  rename(BW_O2_adj_P.Val = adj.P.Val)

EP_ISO_ME = read_csv("RNA_Seq_Output/Dream_RawFiles/EP_ME_dreamISO_o2DE.csv")
colnames(EP_ISO_ME)[1] = "Gene_ID"

EP_ISO_ME = EP_ISO_ME %>%
  select(Gene_ID, logFC, P.Value, adj.P.Val) %>%
  mutate(ME_O2_SIG = ifelse(adj.P.Val < 0.05, "SIG", "NA")) %>% 
  rename(ME_O2_logFC = logFC) %>% 
  rename(ME_O2_P.Val = P.Value) %>% 
  rename(ME_O2_adj_P.Val = adj.P.Val)

EP_ISO_1 = left_join(EP_ISO_Strain, EP_ISO_O2, by = "Gene_ID")
EP_ISO_2 = left_join(EP_ISO_1, EP_ISO_IXN, by = "Gene_ID" )
EP_ISO_3 = left_join(EP_ISO_2, EP_ISO_BW, by = "Gene_ID")
EP_ISO_4 = left_join(EP_ISO_3, EP_ISO_ME, by = "Gene_ID")

EP_ISO_5 = left_join(EP_ISO_4, PBS, by = "Gene_ID") %>% 
  mutate(PBS = ifelse(is.na(PBS), "NA", PBS))

EP_ISO_Summary = left_join(EP_ISO_5, Apriori, by = "Gene_ID") %>% 
  mutate(APRI = ifelse(is.na(APRI), "NA", PBS))

# write.csv(EP_ISO_Summary, "RNA_Seq_Output/Dream_RawFiles/EP_ISO_Summary.csv")

Orthologs = read_csv("Dream_Raw_Data/EP_Pman_orthologs.csv")
Orthologs = Orthologs %>%
  select(Gene_ID, P_maniculatus_ID, P_maniculatus_attr, M_musculus_ID, M_musculus_attr)

EP_ISO_Summary <- EP_ISO_Summary %>%
  left_join(Orthologs, by = "Gene_ID") %>%
  relocate(25:27, .after = 1) %>% # Move columns with "Ortholog" to positions 2-5
  select(-P_maniculatus_ID)

EP_ISO_Summary <- EP_ISO_Summary %>%
  mutate(
    Pman_GeneID = ifelse(grepl("^LOC", Gene_ID), M_musculus_ID, 
                         ifelse(!is.na(Gene_ID), Gene_ID, M_musculus_ID)),
    Pman_GeneID = ifelse(is.na(Pman_GeneID), Gene_ID, Pman_GeneID)
  ) %>%
  relocate(Pman_GeneID, .before = everything())

# Summary
sum(EP_ISO_Summary$strain_SIG == "SIG")
sum(EP_ISO_Summary$O2_SIG == "SIG")
sum(EP_ISO_Summary$IXN_SIG == "SIG")
sum(EP_ISO_Summary$BW_O2_SIG == "SIG")
sum(EP_ISO_Summary$ME_O2_SIG == "SIG")
sum(EP_ISO_Summary$PBS == "TRUE")
sum(EP_ISO_Summary$APRI == "TRUE")

#write.csv(EP_ISO_Summary, "Dream_Output/EP_ISO_Summary.csv")


# IsoQuant Late Pregnancy -----------
# Filtered at 0.5 CPM

LP_ISO_Strain = read_csv("Dream_Output/Dream_RawFiles/LP_dreamISO_strainDE.csv")
colnames(LP_ISO_Strain)[1] = "Gene_ID"

LP_ISO_Strain = LP_ISO_Strain %>%
  select(Gene_ID, logFC, P.Value, adj.P.Val) %>%
  mutate(strain_SIG = ifelse(adj.P.Val < 0.05, "SIG", "NA")) %>%
  rename(strain_logFC = logFC) %>% 
  rename(strain_P.Val = P.Value) %>% 
  rename(strain_adj_P.Val = adj.P.Val)

LP_ISO_O2 = read_csv("Dream_Output/Dream_RawFiles/LP_dreamISO_o2DE.csv")
colnames(LP_ISO_O2)[1] = "Gene_ID"

LP_ISO_O2 = LP_ISO_O2 %>%
  select(Gene_ID, logFC, P.Value, adj.P.Val) %>%
  mutate(O2_SIG = ifelse(adj.P.Val < 0.05, "SIG", "NA")) %>% 
  rename(O2_logFC = logFC) %>% 
  rename(O2_P.Val = P.Value) %>% 
  rename(O2_adj_P.Val = adj.P.Val)

LP_ISO_IXN = read_csv("Dream_Output/Dream_RawFiles/LP_dreamISO_ixnDE.csv")
colnames(LP_ISO_IXN)[1] = "Gene_ID"

LP_ISO_IXN = LP_ISO_IXN %>%
  select(Gene_ID, logFC, P.Value, adj.P.Val) %>%
  mutate(IXN_SIG = ifelse(adj.P.Val < 0.05, "SIG", "NA")) %>% 
  rename(IXN_logFC = logFC) %>% 
  rename(IXN_P.Val = P.Value) %>% 
  rename(IXN_adj_P.Val = adj.P.Val)

LP_ISO_BW = read_csv("Dream_Output/Dream_RawFiles/LP_BW_dreamISO_o2DE.csv")
colnames(LP_ISO_BW)[1] = "Gene_ID"

LP_ISO_BW = LP_ISO_BW %>%
  select(Gene_ID, logFC, P.Value, adj.P.Val) %>%
  mutate(BW_O2_SIG = ifelse(adj.P.Val < 0.05, "SIG", "NA")) %>% 
  rename(BW_O2_logFC = logFC) %>% 
  rename(BW_O2_P.Val = P.Value) %>% 
  rename(BW_O2_adj_P.Val = adj.P.Val)

LP_ISO_ME = read_csv("Dream_Output/Dream_RawFiles/LP_ME_dreamISO_o2DE.csv")
colnames(LP_ISO_ME)[1] = "Gene_ID"

LP_ISO_ME = LP_ISO_ME %>%
  select(Gene_ID, logFC, P.Value, adj.P.Val) %>%
  mutate(ME_O2_SIG = ifelse(adj.P.Val < 0.05, "SIG", "NA")) %>% 
  rename(ME_O2_logFC = logFC) %>% 
  rename(ME_O2_P.Val = P.Value) %>% 
  rename(ME_O2_adj_P.Val = adj.P.Val)

LP_ISO_1 = left_join(LP_ISO_Strain, LP_ISO_O2, by = "Gene_ID")
LP_ISO_2 = left_join(LP_ISO_1, LP_ISO_IXN, by = "Gene_ID" )
LP_ISO_3 = left_join(LP_ISO_2, LP_ISO_BW, by = "Gene_ID")
LP_ISO_4 = left_join(LP_ISO_3, LP_ISO_ME, by = "Gene_ID")

LP_ISO_5 = left_join(LP_ISO_4, PBS, by = "Gene_ID") %>% 
  mutate(PBS = ifelse(is.na(PBS), "NA", PBS))

LP_ISO_Summary = left_join(LP_ISO_5, Apriori, by = "Gene_ID") %>% 
  mutate(APRI = ifelse(is.na(APRI), "NA", PBS))

Orthologs = read_csv("Dream_Raw_Data/LP_Pman_orthologs.csv")
Orthologs = Orthologs %>%
  select(Gene_ID, P_maniculatus_ID, P_maniculatus_attr, M_musculus_ID, M_musculus_attr)

LP_ISO_Summary <- LP_ISO_Summary %>%
  left_join(Orthologs, by = "Gene_ID") %>%
  relocate(25:27, .after = 1) %>% # Move columns with "Ortholog" to positions 2-5
  select(-P_maniculatus_ID)

LP_ISO_Summary <- LP_ISO_Summary %>%
  mutate(
    Pman_GeneID = ifelse(grepl("^LOC", Gene_ID), M_musculus_ID, 
                         ifelse(!is.na(Gene_ID), Gene_ID, M_musculus_ID)),
    Pman_GeneID = ifelse(is.na(Pman_GeneID), Gene_ID, Pman_GeneID)
  ) %>%
  relocate(Pman_GeneID, .before = everything())

# Summary
sum(LP_ISO_Summary$strain_SIG == "SIG")
sum(LP_ISO_Summary$O2_SIG == "SIG")
sum(LP_ISO_Summary$IXN_SIG == "SIG")
sum(LP_ISO_Summary$BW_O2_SIG == "SIG")
sum(LP_ISO_Summary$ME_O2_SIG == "SIG")
sum(LP_ISO_Summary$PBS == "TRUE")
sum(LP_ISO_Summary$APRI == "TRUE")

#write.csv(LP_ISO_Summary, "Dream_Output/LP_ISO_Summary.csv")


# Overlap between EP and LP in ISOQUANT


count(intersect(EP_ISO_Strain
                %>% filter(strain_SIG == "SIG") 
                %>% select(Gene_ID), 
                LP_ISO_Strain %>% 
                  filter(strain_SIG == "SIG") 
                %>% select(Gene_ID)))

count(intersect(EP_ISO_O2
                %>% filter(O2_SIG == "SIG") 
                %>% select(Gene_ID), 
                LP_ISO_O2 %>% 
                  filter(O2_SIG == "SIG") 
                %>% select(Gene_ID)))

O2_overlap = intersect(EP_ISO_O2
                             %>% filter(O2_SIG == "SIG") 
                             %>% select(Gene_ID), 
                             LP_ISO_O2 %>% 
                               filter(O2_SIG == "SIG") 
                             %>% select(Gene_ID))

count(intersect(EP_ISO_IXN
                %>% filter(IXN_SIG == "SIG") 
                %>% select(Gene_ID), 
                LP_ISO_IXN %>% 
                  filter(IXN_SIG == "SIG") 
                %>% select(Gene_ID)))

ixn_overlap = intersect(EP_ISO_IXN
                %>% filter(IXN_SIG == "SIG") 
                %>% select(Gene_ID), 
                LP_ISO_IXN %>% 
                  filter(IXN_SIG == "SIG") 
                %>% select(Gene_ID))

count(intersect(EP_ISO_BW
                %>% filter(BW_O2_SIG == "SIG") 
                %>% select(Gene_ID), 
                LP_ISO_BW %>% 
                  filter(BW_O2_SIG == "SIG") 
                %>% select(Gene_ID)))

BW_overlap = intersect(EP_ISO_BW
                %>% filter(BW_O2_SIG == "SIG") 
                %>% select(Gene_ID), 
                LP_ISO_BW %>% 
                  filter(BW_O2_SIG == "SIG") 
                %>% select(Gene_ID))

count(intersect(EP_ISO_ME
                %>% filter(ME_O2_SIG == "SIG") 
                %>% select(Gene_ID), 
                LP_ISO_ME %>% 
                  filter(ME_O2_SIG == "SIG") 
                %>% select(Gene_ID)))




# GCF Early Pregnancy -----------
# Filtered at 0.5 CPM

EP_GCF_Strain = read_csv("Dream_Output/Dream_RawFiles/EP_dreamGCF_strainDE.csv")
colnames(EP_GCF_Strain)[1] = "Gene_ID"

EP_GCF_Strain = EP_GCF_Strain %>%
  select(Gene_ID, logFC, P.Value, adj.P.Val) %>%
  mutate(strain_SIG = ifelse(adj.P.Val < 0.05, "SIG", "NA")) %>%
  rename(strain_logFC = logFC) %>% 
  rename(strain_P.Val = P.Value) %>% 
  rename(strain_adj_P.Val = adj.P.Val)

EP_GCF_O2 = read_csv("Dream_Output/Dream_RawFiles/EP_dreamGCF_o2DE.csv")
colnames(EP_GCF_O2)[1] = "Gene_ID"

EP_GCF_O2 = EP_GCF_O2 %>%
  select(Gene_ID, logFC, P.Value, adj.P.Val) %>%
  mutate(O2_SIG = ifelse(adj.P.Val < 0.05, "SIG", "NA")) %>% 
  rename(O2_logFC = logFC) %>% 
  rename(O2_P.Val = P.Value) %>% 
  rename(O2_adj_P.Val = adj.P.Val)

EP_GCF_IXN = read_csv("Dream_Output/Dream_RawFiles/EP_dreamGCF_ixnDE.csv")
colnames(EP_GCF_IXN)[1] = "Gene_ID"

EP_GCF_IXN = EP_GCF_IXN %>%
  select(Gene_ID, logFC, P.Value, adj.P.Val) %>%
  mutate(IXN_SIG = ifelse(adj.P.Val < 0.05, "SIG", "NA")) %>% 
  rename(IXN_logFC = logFC) %>% 
  rename(IXN_P.Val = P.Value) %>% 
  rename(IXN_adj_P.Val = adj.P.Val)

EP_GCF_BW = read_csv("Dream_Output/Dream_RawFiles/EP_BW_dreamGCF_o2DE.csv")
colnames(EP_GCF_BW)[1] = "Gene_ID"

EP_GCF_BW = EP_GCF_BW %>%
  select(Gene_ID, logFC, P.Value, adj.P.Val) %>%
  mutate(BW_O2_SIG = ifelse(adj.P.Val < 0.05, "SIG", "NA")) %>% 
  rename(BW_O2_logFC = logFC) %>% 
  rename(BW_O2_P.Val = P.Value) %>% 
  rename(BW_O2_adj_P.Val = adj.P.Val)

EP_GCF_ME = read_csv("Dream_Output/Dream_RawFiles/EP_ME_dreamGCF_o2DE.csv")
colnames(EP_GCF_ME)[1] = "Gene_ID"

EP_GCF_ME = EP_GCF_ME %>%
  select(Gene_ID, logFC, P.Value, adj.P.Val) %>%
  mutate(ME_O2_SIG = ifelse(adj.P.Val < 0.05, "SIG", "NA")) %>% 
  rename(ME_O2_logFC = logFC) %>% 
  rename(ME_O2_P.Val = P.Value) %>% 
  rename(ME_O2_adj_P.Val = adj.P.Val)

EP_GCF_1 = left_join(EP_GCF_Strain, EP_GCF_O2, by = "Gene_ID")
EP_GCF_2 = left_join(EP_GCF_1, EP_GCF_IXN, by = "Gene_ID" )
EP_GCF_3 = left_join(EP_GCF_2, EP_GCF_BW, by = "Gene_ID")
EP_GCF_4 = left_join(EP_GCF_3, EP_GCF_ME, by = "Gene_ID")

EP_GCF_5 = left_join(EP_GCF_4, PBS, by = "Gene_ID") %>% 
  mutate(PBS = ifelse(is.na(PBS), "NA", PBS))

EP_GCF_Summary = left_join(EP_GCF_5, Apriori, by = "Gene_ID") %>% 
  mutate(APRI = ifelse(is.na(PBS), "NA", PBS))

# Summary
sum(EP_GCF_Summary$strain_SIG == "SIG")
sum(EP_GCF_Summary$O2_SIG == "SIG")
sum(EP_GCF_Summary$IXN_SIG == "SIG")
sum(EP_GCF_Summary$BW_O2_SIG == "SIG")
sum(EP_GCF_Summary$ME_O2_SIG == "SIG")


# GCF Late Pregnancy -----------
# Filtered at 0.5 CPM

LP_GCF_Strain = read_csv("Dream_Output/Dream_RawFiles/LP_dreamGCF_strainDE.csv")
colnames(LP_GCF_Strain)[1] = "Gene_ID"

LP_GCF_Strain = LP_GCF_Strain %>%
  select(Gene_ID, logFC, P.Value, adj.P.Val) %>%
  mutate(strain_SIG = ifelse(adj.P.Val < 0.05, "SIG", "NA")) %>%
  rename(strain_logFC = logFC) %>% 
  rename(strain_P.Val = P.Value) %>% 
  rename(strain_adj_P.Val = adj.P.Val)

LP_GCF_O2 = read_csv("Dream_Output/Dream_RawFiles/LP_dreamGCF_o2DE.csv")
colnames(LP_GCF_O2)[1] = "Gene_ID"

LP_GCF_O2 = LP_GCF_O2 %>%
  select(Gene_ID, logFC, P.Value, adj.P.Val) %>%
  mutate(O2_SIG = ifelse(adj.P.Val < 0.05, "SIG", "NA")) %>% 
  rename(O2_logFC = logFC) %>% 
  rename(O2_P.Val = P.Value) %>% 
  rename(O2_adj_P.Val = adj.P.Val)

LP_GCF_IXN = read_csv("Dream_Output/Dream_RawFiles/LP_dreamGCF_ixnDE.csv")
colnames(LP_GCF_IXN)[1] = "Gene_ID"

LP_GCF_IXN = LP_GCF_IXN %>%
  select(Gene_ID, logFC, P.Value, adj.P.Val) %>%
  mutate(IXN_SIG = ifelse(adj.P.Val < 0.05, "SIG", "NA")) %>% 
  rename(IXN_logFC = logFC) %>% 
  rename(IXN_P.Val = P.Value) %>% 
  rename(IXN_adj_P.Val = adj.P.Val)

LP_GCF_BW = read_csv("Dream_Output/Dream_RawFiles/LP_BW_dreamGCF_o2DE.csv")
colnames(LP_GCF_BW)[1] = "Gene_ID"

LP_GCF_BW = LP_GCF_BW %>%
  select(Gene_ID, logFC, P.Value, adj.P.Val) %>%
  mutate(BW_O2_SIG = ifelse(adj.P.Val < 0.05, "SIG", "NA")) %>% 
  rename(BW_O2_logFC = logFC) %>% 
  rename(BW_O2_P.Val = P.Value) %>% 
  rename(BW_O2_adj_P.Val = adj.P.Val)

LP_GCF_ME = read_csv("Dream_Output/Dream_RawFiles/LP_ME_dreamGCF_o2DE.csv")
colnames(LP_GCF_ME)[1] = "Gene_ID"

LP_GCF_ME = LP_GCF_ME %>%
  select(Gene_ID, logFC, P.Value, adj.P.Val) %>%
  mutate(ME_O2_SIG = ifelse(adj.P.Val < 0.05, "SIG", "NA")) %>% 
  rename(ME_O2_logFC = logFC) %>% 
  rename(ME_O2_P.Val = P.Value) %>% 
  rename(ME_O2_adj_P.Val = adj.P.Val)

LP_GCF_1 = left_join(LP_GCF_Strain, LP_GCF_O2, by = "Gene_ID")
LP_GCF_2 = left_join(LP_GCF_1, LP_GCF_IXN, by = "Gene_ID" )
LP_GCF_3 = left_join(LP_GCF_2, LP_GCF_BW, by = "Gene_ID")
LP_GCF_4 = left_join(LP_GCF_3, LP_GCF_ME, by = "Gene_ID")

LP_GCF_5 = left_join(LP_GCF_4, PBS, by = "Gene_ID") %>% 
  mutate(PBS = ifelse(is.na(PBS), "NA", PBS))

LP_GCF_Summary = left_join(LP_GCF_5, Apriori, by = "Gene_ID") %>% 
  mutate(APRI = ifelse(is.na(PBS), "NA", PBS))

# Summary
sum(LP_GCF_Summary$strain_SIG == "SIG")
sum(LP_GCF_Summary$O2_SIG == "SIG")
sum(LP_GCF_Summary$IXN_SIG == "SIG")
sum(LP_GCF_Summary$BW_O2_SIG == "SIG")
sum(LP_GCF_Summary$ME_O2_SIG == "SIG")


# Overlap between GCF and Isoquant-----------
# Both filtered at 0.05

# Strain
GCF_strain = LP_GCF_Strain %>% 
  filter(strain_SIG == "SIG") %>% 
  select(Gene_ID)

ISO_strain = LP_ISO_Strain %>% 
  filter(strain_SIG == "SIG") %>% 
  select(Gene_ID)

overlap_strain = intersect(GCF_strain, ISO_strain)
count(overlap_strain)

# Hypoxia
GCF_O2 = LP_GCF_O2 %>% 
  filter(O2_SIG == "SIG") %>% 
  select(Gene_ID)

ISO_O2 = LP_ISO_O2 %>% 
  filter(O2_SIG == "SIG") %>% 
  select(Gene_ID)

overlap_O2 = intersect(GCF_O2, ISO_O2)
count(overlap_O2)

# IXN
GCF_IXN = LP_GCF_IXN %>% 
  filter(IXN_SIG == "SIG") %>% 
  select(Gene_ID)

ISO_IXN = LP_ISO_IXN %>% 
  filter(IXN_SIG == "SIG") %>% 
  select(Gene_ID)

overlap_IXN = intersect(GCF_IXN, ISO_IXN)
count(overlap_IXN)

# BW
GCF_BW = LP_GCF_BW %>% 
  filter(BW_O2_SIG == "SIG") %>% 
  select(Gene_ID)

ISO_BW = LP_ISO_BW %>% 
  filter(BW_O2_SIG == "SIG") %>% 
  select(Gene_ID)

overlap_BW = intersect(GCF_BW, ISO_BW)
count(overlap_BW)

# ME
GCF_ME = LP_GCF_ME %>% 
  filter(ME_O2_SIG == "SIG") %>% 
  select(Gene_ID)

ISO_ME = LP_ISO_ME %>% 
  filter(ME_O2_SIG == "SIG") %>% 
  select(Gene_ID)

overlap_ME = intersect(GCF_ME, ISO_ME)
count(overlap_ME)

# Within GCF -- Hypoxia vs BW/ME Only
overlap_GCF_O2 = intersect(GCF_O2, GCF_BW)
count(overlap_GCF_O2)

overlap_GCF_O2 = intersect(GCF_O2, GCF_ME)
count(overlap_GCF_O2)

# Within ISO -- Hypoxia vs BW/ME Only
overlap_ISO_O2 = intersect(ISO_O2, ISO_BW)
count(overlap_ISO_O2)

overlap_ISO_O2 = intersect(ISO_O2, ISO_ME)
count(overlap_ISO_O2)


