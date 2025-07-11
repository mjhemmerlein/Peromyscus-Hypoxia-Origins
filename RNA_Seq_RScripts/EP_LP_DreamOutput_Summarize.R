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

Orthologs = read_csv("RNA_Seq_RawData/EP_orthologs.csv")
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

# write.csv(EP_ISO_Summary, "RNA_Seq_Output/EP_ISO_Ortho_Summary.csv")


# IsoQuant Late Pregnancy -----------
# Filtered at 0.5 CPM

LP_ISO_Strain = read_csv("RNA_Seq_Output/Dream_RawFiles/LP_dreamISO_strainDE.csv")
colnames(LP_ISO_Strain)[1] = "Gene_ID"

LP_ISO_Strain = LP_ISO_Strain %>%
  select(Gene_ID, logFC, P.Value, adj.P.Val) %>%
  mutate(strain_SIG = ifelse(adj.P.Val < 0.05, "SIG", "NA")) %>%
  rename(strain_logFC = logFC) %>% 
  rename(strain_P.Val = P.Value) %>% 
  rename(strain_adj_P.Val = adj.P.Val)

LP_ISO_O2 = read_csv("RNA_Seq_Output/Dream_RawFiles/LP_dreamISO_o2DE.csv")
colnames(LP_ISO_O2)[1] = "Gene_ID"

LP_ISO_O2 = LP_ISO_O2 %>%
  select(Gene_ID, logFC, P.Value, adj.P.Val) %>%
  mutate(O2_SIG = ifelse(adj.P.Val < 0.05, "SIG", "NA")) %>% 
  rename(O2_logFC = logFC) %>% 
  rename(O2_P.Val = P.Value) %>% 
  rename(O2_adj_P.Val = adj.P.Val)

LP_ISO_IXN = read_csv("RNA_Seq_Output/Dream_RawFiles/LP_dreamISO_ixnDE.csv")
colnames(LP_ISO_IXN)[1] = "Gene_ID"

LP_ISO_IXN = LP_ISO_IXN %>%
  select(Gene_ID, logFC, P.Value, adj.P.Val) %>%
  mutate(IXN_SIG = ifelse(adj.P.Val < 0.05, "SIG", "NA")) %>% 
  rename(IXN_logFC = logFC) %>% 
  rename(IXN_P.Val = P.Value) %>% 
  rename(IXN_adj_P.Val = adj.P.Val)

LP_ISO_BW = read_csv("RNA_Seq_Output/Dream_RawFiles/LP_BW_dreamISO_o2DE.csv")
colnames(LP_ISO_BW)[1] = "Gene_ID"

LP_ISO_BW = LP_ISO_BW %>%
  select(Gene_ID, logFC, P.Value, adj.P.Val) %>%
  mutate(BW_O2_SIG = ifelse(adj.P.Val < 0.05, "SIG", "NA")) %>% 
  rename(BW_O2_logFC = logFC) %>% 
  rename(BW_O2_P.Val = P.Value) %>% 
  rename(BW_O2_adj_P.Val = adj.P.Val)

LP_ISO_ME = read_csv("RNA_Seq_Output/Dream_RawFiles/LP_ME_dreamISO_o2DE.csv")
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

# write.csv(LP_ISO_Summary, "RNA_Seq_Output/Dream_RawFiles/LP_LZ_ISO_Summary.csv")

Orthologs = read_csv("RNA_Seq_RawData/LP_orthologs.csv")
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

# write.csv(LP_ISO_Summary, "RNA_Seq_Output//LP_ISO_Ortho_Summary.csv")


# IsoQuant Late Pregnancy JZ -----------
# Filtered at 0.5 CPM

LP_JZ_Strain = read_csv("RNA_Seq_Output/Dream_RawFiles/LP_JZ_dreamISO_strainDE.csv")
colnames(LP_JZ_Strain)[1] = "Gene_ID"

LP_JZ_Strain = LP_JZ_Strain %>%
  select(Gene_ID, logFC, P.Value, adj.P.Val) %>%
  mutate(strain_SIG = ifelse(adj.P.Val < 0.05, "SIG", "NA")) %>%
  rename(strain_logFC = logFC) %>% 
  rename(strain_P.Val = P.Value) %>% 
  rename(strain_adj_P.Val = adj.P.Val)

LP_JZ_O2 = read_csv("RNA_Seq_Output/Dream_RawFiles/LP_JZ_dreamISO_o2DE.csv")
colnames(LP_JZ_O2)[1] = "Gene_ID"

LP_JZ_O2 = LP_JZ_O2 %>%
  select(Gene_ID, logFC, P.Value, adj.P.Val) %>%
  mutate(O2_SIG = ifelse(adj.P.Val < 0.05, "SIG", "NA")) %>% 
  rename(O2_logFC = logFC) %>% 
  rename(O2_P.Val = P.Value) %>% 
  rename(O2_adj_P.Val = adj.P.Val)

LP_JZ_IXN = read_csv("RNA_Seq_Output/Dream_RawFiles/LP_JZ_dreamISO_ixnDE.csv")
colnames(LP_JZ_IXN)[1] = "Gene_ID"

LP_JZ_IXN = LP_JZ_IXN %>%
  select(Gene_ID, logFC, P.Value, adj.P.Val) %>%
  mutate(IXN_SIG = ifelse(adj.P.Val < 0.05, "SIG", "NA")) %>% 
  rename(IXN_logFC = logFC) %>% 
  rename(IXN_P.Val = P.Value) %>% 
  rename(IXN_adj_P.Val = adj.P.Val)

LP_JZ_BW = read_csv("RNA_Seq_Output/Dream_RawFiles/LP_JZ_BW_dreamISO_o2DE.csv")
colnames(LP_JZ_BW)[1] = "Gene_ID"

LP_JZ_BW = LP_JZ_BW %>%
  select(Gene_ID, logFC, P.Value, adj.P.Val) %>%
  mutate(BW_O2_SIG = ifelse(adj.P.Val < 0.05, "SIG", "NA")) %>% 
  rename(BW_O2_logFC = logFC) %>% 
  rename(BW_O2_P.Val = P.Value) %>% 
  rename(BW_O2_adj_P.Val = adj.P.Val)

LP_JZ_ME = read_csv("RNA_Seq_Output/Dream_RawFiles/LP_JZ_ME_dreamISO_o2DE.csv")
colnames(LP_JZ_ME)[1] = "Gene_ID"

LP_JZ_ME = LP_JZ_ME %>%
  select(Gene_ID, logFC, P.Value, adj.P.Val) %>%
  mutate(ME_O2_SIG = ifelse(adj.P.Val < 0.05, "SIG", "NA")) %>% 
  rename(ME_O2_logFC = logFC) %>% 
  rename(ME_O2_P.Val = P.Value) %>% 
  rename(ME_O2_adj_P.Val = adj.P.Val)

LP_JZ_1 = left_join(LP_JZ_Strain, LP_JZ_O2, by = "Gene_ID")
LP_JZ_2 = left_join(LP_JZ_1, LP_JZ_IXN, by = "Gene_ID" )
LP_JZ_3 = left_join(LP_JZ_2, LP_JZ_BW, by = "Gene_ID")
LP_JZ_4 = left_join(LP_JZ_3, LP_JZ_ME, by = "Gene_ID")

LP_JZ_5 = left_join(LP_JZ_4, PBS, by = "Gene_ID") %>% 
  mutate(PBS = ifelse(is.na(PBS), "NA", PBS))

LP_JZ_Summary = left_join(LP_JZ_5, Apriori, by = "Gene_ID") %>% 
  mutate(APRI = ifelse(is.na(APRI), "NA", PBS))

# write.csv(LP_JZ_Summary, "RNA_Seq_Output/Dream_RawFiles/LP_JZ_ISO_Summary.csv")

Orthologs = read_csv("RNA_Seq_Output/LP_orthologs.csv")
Orthologs = Orthologs %>%
  select(Gene_ID, P_maniculatus_ID, P_maniculatus_attr, M_musculus_ID, M_musculus_attr)

LP_JZ_Summary <- LP_JZ_Summary %>%
  left_join(Orthologs, by = "Gene_ID") %>%
  relocate(25:27, .after = 1) %>% # Move columns with "Ortholog" to positions 2-5
  select(-P_maniculatus_ID)

LP_JZ_Summary <- LP_JZ_Summary %>%
  mutate(
    Pman_GeneID = ifelse(grepl("^LOC", Gene_ID), M_musculus_ID, 
                         ifelse(!is.na(Gene_ID), Gene_ID, M_musculus_ID)),
    Pman_GeneID = ifelse(is.na(Pman_GeneID), Gene_ID, Pman_GeneID)
  ) %>%
  relocate(Pman_GeneID, .before = everything())

# Summary
sum(LP_JZ_Summary$strain_SIG == "SIG")
sum(LP_JZ_Summary$O2_SIG == "SIG")
sum(LP_JZ_Summary$IXN_SIG == "SIG")
sum(LP_JZ_Summary$BW_O2_SIG == "SIG")
sum(LP_JZ_Summary$ME_O2_SIG == "SIG")
sum(LP_JZ_Summary$PBS == "TRUE")
sum(LP_JZ_Summary$APRI == "TRUE")

# write.csv(LP_ISO_Summary, "RNA_Seq_Output//LP_ISO_Ortho_Summary.csv")




# Overlap between EP and LP in ISOQUANT

# EP
O2_overlap <- EP_ISO_Summary %>%
  filter(O2_SIG == "SIG", BW_O2_SIG == "SIG") %>%
  select(Gene_ID)

O2_overlap <- EP_ISO_Summary %>%
  filter(O2_SIG == "SIG", ME_O2_SIG == "SIG") %>%
  select(Gene_ID)

# LP LZ
O2_overlap <- LP_ISO_Summary %>%
  filter(O2_SIG == "SIG", BW_O2_SIG == "SIG") %>%
  select(Gene_ID)

O2_overlap <- LP_ISO_Summary %>%
  filter(O2_SIG == "SIG", ME_O2_SIG == "SIG") %>%
  select(Gene_ID)

# LP JZ
O2_overlap <- LP_JZ_Summary %>%
  filter(O2_SIG == "SIG", BW_O2_SIG == "SIG") %>%
  select(Gene_ID)

O2_overlap <- LP_JZ_Summary %>%
  filter(O2_SIG == "SIG", ME_O2_SIG == "SIG") %>%
  select(Gene_ID)






count_sig_overlap <- function(df1, df2, sig_col) {
  genes1 <- df1 %>% filter(.data[[sig_col]] == "SIG") %>% pull(Gene_ID)
  genes2 <- df2 %>% filter(.data[[sig_col]] == "SIG") %>% pull(Gene_ID)
  length(intersect(genes1, genes2))
}

count_sig_overlap(EP_ISO_Strain, LP_ISO_Strain, "strain_SIG")
count_sig_overlap(EP_ISO_O2, LP_ISO_O2, "O2_SIG")
count_sig_overlap(EP_ISO_IXN, LP_ISO_IXN, "IXN_SIG")
count_sig_overlap(EP_ISO_BW, LP_ISO_BW, "BW_O2_SIG")
count_sig_overlap(EP_ISO_ME, LP_ISO_ME, "ME_O2_SIG")

count_sig_overlap(EP_ISO_Strain, LP_JZ_Strain, "strain_SIG")
count_sig_overlap(EP_ISO_O2, LP_JZ_O2, "O2_SIG")
count_sig_overlap(EP_ISO_IXN, LP_JZ_IXN, "IXN_SIG")
count_sig_overlap(EP_ISO_BW, LP_JZ_BW, "BW_O2_SIG")
count_sig_overlap(EP_ISO_ME, LP_JZ_ME, "ME_O2_SIG")

count_sig_overlap(LP_ISO_Strain, LP_JZ_Strain, "strain_SIG")
count_sig_overlap(LP_ISO_O2, LP_JZ_O2, "O2_SIG")
count_sig_overlap(LP_ISO_IXN, LP_JZ_IXN, "IXN_SIG")
count_sig_overlap(LP_ISO_BW, LP_JZ_BW, "BW_O2_SIG")
count_sig_overlap(LP_ISO_ME, LP_JZ_ME, "ME_O2_SIG")





library(UpSetR)

# Make sure the data frame is clean and logical
upset_data <- data.frame(
  Gene_ID = unique(c(
    EP_ISO_Strain$Gene_ID,
    LP_ISO_Strain$Gene_ID
  ))
) %>%
  mutate(
    EP = Gene_ID %in% (EP_ISO_Strain %>% filter(strain_SIG == "SIG") %>% pull(Gene_ID)),
    LP = Gene_ID %in% (LP_ISO_Strain %>% filter(strain_SIG == "SIG") %>% pull(Gene_ID)),
  )

upset_data <- upset_data %>%
  mutate(across(c(EP, LP), as.integer))

upset(upset_data[, c("EP", "LP")], sets = c("EP", "LP"))


# Make sure the data frame is clean and logical
upset_data <- data.frame(
  Gene_ID = unique(c(
    EP_ISO_O2$Gene_ID,
    LP_ISO_O2$Gene_ID
  ))
) %>%
  mutate(
    EP = Gene_ID %in% (EP_ISO_O2 %>% filter(O2_SIG == "SIG") %>% pull(Gene_ID)),
    LP = Gene_ID %in% (LP_ISO_O2 %>% filter(O2_SIG == "SIG") %>% pull(Gene_ID)),
  )

upset_data <- upset_data %>%
  mutate(across(c(EP, LP), as.integer))

upset(upset_data[, c("EP", "LP")], sets = c("EP", "LP"))


# Make sure the data frame is clean and logical
upset_data <- data.frame(
  Gene_ID = unique(c(
    EP_ISO_O2$Gene_ID,
    EP_ISO_BW$Gene_ID,
    EP_ISO_ME$Gene_ID
  ))
) %>%
  mutate(
    EP = Gene_ID %in% (EP_ISO_O2 %>% filter(O2_SIG == "SIG") %>% pull(Gene_ID)),
    EP_BW = Gene_ID %in% (EP_ISO_BW %>% filter(BW_O2_SIG == "SIG") %>% pull(Gene_ID)),
    EP_ME = Gene_ID %in% (EP_ISO_ME %>% filter(ME_O2_SIG == "SIG") %>% pull(Gene_ID))
  )

upset_data <- upset_data %>%
  mutate(across(c(EP, EP_BW, EP_ME), as.integer))

upset(upset_data[, c("EP", "EP_BW", "EP_ME")], sets = c("EP", "EP_BW", "EP_ME"))

# Make sure the data frame is clean and logical
upset_data <- data.frame(
  Gene_ID = unique(c(
    LP_ISO_O2$Gene_ID,
    LP_ISO_BW$Gene_ID,
    LP_ISO_ME$Gene_ID
  ))
) %>%
  mutate(
    LP = Gene_ID %in% (LP_ISO_O2 %>% filter(O2_SIG == "SIG") %>% pull(Gene_ID)),
    LP_BW = Gene_ID %in% (LP_ISO_BW %>% filter(BW_O2_SIG == "SIG") %>% pull(Gene_ID)),
    LP_ME = Gene_ID %in% (LP_ISO_ME %>% filter(ME_O2_SIG == "SIG") %>% pull(Gene_ID))
  )

upset_data <- upset_data %>%
  mutate(across(c(LP, LP_BW, LP_ME), as.integer))

upset(upset_data[, c("LP", "LP_BW", "LP_ME")], sets = c("LP", "LP_BW", "LP_ME"))
