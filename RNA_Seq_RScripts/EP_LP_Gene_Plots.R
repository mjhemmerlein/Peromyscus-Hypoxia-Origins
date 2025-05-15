
# Packages
library('BiocParallel')
library('dplyr')
library('ggplot2')
library('lme4')
library('readxl')
library('variancePartition')
library('edgeR')
library(forcats)

# Early Pregnancy -----
EP_Sample_Info <- read_xlsx("RNA_Seq_RawData/MetaData_EP.xlsx")
EP_Sample_Info = as.data.frame(EP_Sample_Info)
rownames(EP_Sample_Info) = EP_Sample_Info$Sample_ID
EP_Sample_Info = EP_Sample_Info %>% filter(Strain == "BW" | Strain == "ME")

# Read in Files + QC
Pman_rawreads <- read_xlsx("RNA_Seq_RawData/EP_Pman_ExtMMFrac_readcounts.xlsx")
Pman_rawreads = as.data.frame(Pman_rawreads)
Pman_rawreads = na.omit(Pman_rawreads, col = "Geneid") # 1 NA
Pman_rawreads <- `row.names<-`(Pman_rawreads, Pman_rawreads$Geneid)
Pman_rawreads <- Pman_rawreads[,-c(1:6)]

# Check read table vs sample info
Check = EP_Sample_Info$Seq_Name
colnames(Pman_rawreads) == Check

colnames(Pman_rawreads) = rownames(EP_Sample_Info)

Pman_readcounts <- as.matrix(Pman_rawreads)
dPman_0 <- DGEList(Pman_readcounts)

dPman_0 <- calcNormFactors(dPman_0)
dim(dPman_0)
keep <- rowSums(cpm(dPman_0) > 0.5 ) >= 60
dPman <- dPman_0[keep,]
dim(dPman)
plotMDS(dPman, col = as.numeric(EP_Sample_Info$Strain), labels = EP_Sample_Info$Strain)

rownames(EP_Sample_Info) == colnames(dPman)

# Interaction
param = SnowParam(8, "SOCK", progressbar=TRUE)
form <- ~ Strain*O2 + (1|Mom)
EP_vobjdream = voomWithDreamWeights(dPman, form, EP_Sample_Info, BPPARAM=param, plot = T)
fitmm = dream( EP_vobjdream, form, EP_Sample_Info, ddf = "Kenward-Roger")
fitmm = eBayes(fitmm)

DE_strain <- topTable( fitmm, coef='StrainME', sort.by = "P", n = Inf, )
DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)
DE_ixn <- topTable( fitmm, coef='StrainME:O22H', sort.by = "P", n = Inf)

length(DE_strain$logFC[which(DE_strain$adj.P.Val < 0.05)])
length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])
length(DE_ixn$logFC[which(DE_ixn$adj.P.Val < 0.05)])


# Late Pregnancy ------------
Sample_Info <- read_xlsx("RNA_Seq_RawData/MetaData_LP.xlsx")
Sample_Info = as.data.frame(Sample_Info)
Sample_Info = Sample_Info[-c(80,81),]
rownames(Sample_Info) = Sample_Info$Sample_ID_LZ
Sample_Info = Sample_Info %>% filter(Strain == "BW" | Strain == "ME")
Sample_Info = Sample_Info[-which(rownames(Sample_Info) == "LZ089"),]

# Read in Files + QC
Pman_rawreads = read_xlsx("RNA_Seq_RawData/LP_Pman_ExtMMFrac_readcounts.xlsx")
Pman_rawreads = as.data.frame(Pman_rawreads)
Pman_rawreads = Pman_rawreads %>%
  filter(!is.na(Geneid))
Pman_rawreads = `row.names<-`(Pman_rawreads, Pman_rawreads$Geneid)
Pman_rawreads <- Pman_rawreads[,-c(1:6)]
Pman_rawreads <- subset(Pman_rawreads, select = -c(RNA201216ZC_LZ089_S16_L001_fastp_pman_Halign_liberal.bam))

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
LP_vobjDream = voomWithDreamWeights(dPman, form, LP_Sample_Info, BPPARAM=param, plot = T)
fitmm = dream( LP_vobjDream, form, LP_Sample_Info, ddf = "Kenward-Roger")
fitmm = eBayes(fitmm)

DE_strain <- topTable( fitmm, coef='StrainME', sort.by = "P", n = Inf, )
DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)
DE_ixn <- topTable( fitmm, coef='StrainME:O22H', sort.by = "P", n = Inf)

length(DE_strain$logFC[which(DE_strain$adj.P.Val < 0.05)])
length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])
length(DE_ixn$logFC[which(DE_ixn$adj.P.Val < 0.05)])







# Plotting functions -------

# Function to combine and prepare data for faceted plotting
prepare_combined_data <- function(gene_id, EP_vobjdream, LP_vobjDream, 
                                  EP_Sample_Info, LP_Sample_Info) {
  # Extract early pregnancy data
  ep_data <- data.frame(
    Expression = EP_vobjdream$E[gene_id, ],
    Group = EP_Sample_Info$Group) %>%
    mutate(
      Group = fct_recode(Group, 
                         "Lowland\nNormoxia" = "BW1N",
                         "Lowland\nHypoxia" = "BW2H", 
                         "Highland\nNormoxia" = "ME1N",
                         "Highland\nHypoxia" = "ME2H"),
      Strain = EP_Sample_Info$Strain,
      Stage = "Early Pregnancy"
    )
  
  # Extract late pregnancy data
  lp_data <- data.frame(
    Expression = LP_vobjDream$E[gene_id, ],
    Group = LP_Sample_Info$Group) %>%
    mutate(
      Group = fct_recode(Group, 
                         "Lowland\nNormoxia" = "BW1N",
                         "Lowland\nHypoxia" = "BW2H", 
                         "Highland\nNormoxia" = "ME1N",
                         "Highland\nHypoxia" = "ME2H"),
      Strain = LP_Sample_Info$Strain,
      Stage = "Late Pregnancy"
    )
  
  # Combine the datasets
  rbind(ep_data, lp_data)
}

# Modified plotting function using facets
plot_gene_expression_faceted <- function(gene_id, gene_name = NULL,
                                         EP_vobjdream, LP_vobjDream,
                                         EP_Sample_Info, LP_Sample_Info) {
  # Prepare combined data
  combined_data <- prepare_combined_data(gene_id, EP_vobjdream, LP_vobjDream,
                                         EP_Sample_Info, LP_Sample_Info)
  
  # Set the gene name for the y-axis
  if (is.null(gene_name)) {
    gene_name <- gene_id
  }
  
  # Create the faceted plot
  ggplot(combined_data, aes(x = Group, y = Expression, fill = Strain)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, color = "black", size = 3.5, alpha = 0.5) +
    facet_grid(. ~ Stage) +
    theme_bw() +
    scale_fill_manual(values = c("goldenrod1", "skyblue")) +
    labs(
      x = NULL,
      y = paste(gene_name, "\nNormalized Transcript\nAbundance")) +
    theme(
      legend.position = "none",
      legend.title = element_blank(),
      strip.background = element_rect(fill = "snow2"),
      strip.text = element_text(size = 20, face = "bold"),
      axis.text.y = element_text(face = "bold", size = 12, color = "black"),
      axis.title.y = element_text(face = "bold", size = 15, color = "black"),
      axis.text.x = element_text(face = "bold", size = 12, color = "black"),
      panel.spacing = unit(1, "lines"))
}

# Function to save faceted plots
save_faceted_plot <- function(gene_id, gene_name = NULL,
                              EP_vobjdream, LP_vobjDream,
                              EP_Sample_Info, LP_Sample_Info,
                              output_dir = "Dream_Output/Faceted_Gene_Plots/") {
  
  plot <- plot_gene_expression_faceted(gene_id, gene_name,
                                       EP_vobjdream, LP_vobjDream,
                                       EP_Sample_Info, LP_Sample_Info)
  
  # Generate filename
  filename <- if (is.null(gene_name)) {
    paste0(output_dir, gene_id, "_faceted.png")
  } else {
    paste0(output_dir, gene_name, "_faceted.png")
  }
  
  # Save the plot
  ggsave(filename, plot, width = 13, height = 5, units = "in", dpi = 300)
}


# Viewing:

# Glucose transporters
plot_gene_expression_faceted("Slc2a1", "Slc2a1", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("Slc1a2", "Slc1a2", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("Slc2a8", "Slc2a8", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("Slc2a12", "Slc2a12", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("Slc22a3", "Slc22a3", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("Slc6a13", "Slc6a13", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("Slc25a23", "Slc25a23", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)

# Transporters
plot_gene_expression_faceted("P2rx1", "P2rx1", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("Clcn2", "Clcn2", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("Trpm6", "Trpm6", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("Tmem150c", "Tmem150c", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("Daam2", "Daam2", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)

# ER transport
plot_gene_expression_faceted("Sec22a", "Sec22a", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)



# Glycolysis

plot_gene_expression_faceted("Gpi", "Gpi", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)

plot_gene_expression_faceted("Pfkm", "Pfkm", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("Pfkl", "Pfkl", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("Pfkfb3", "Pfkfb3", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("Pfkfb1", "Pfkfb1",EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("Pfkfb1", "Pfkfb1", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)

plot_gene_expression_faceted("LOC102925554", "Aldoa-A-like", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("Aldoa", "Aldoa", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)

plot_gene_expression_faceted("LOC102904208", "GAPDH-like", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("LOC102916082", "GAPDH-like", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("Gapdh", "Gapdh", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)

plot_gene_expression_faceted("Gapdhs", "GAPDHS", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)

plot_gene_expression_faceted("Pgam1", "Pgam1", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)

plot_gene_expression_faceted("Eno1", "Eno1", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("Eno4", "Eno4", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)

plot_gene_expression_faceted("LOC102923285", "PKM-like", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("Pklr", "Pklr", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)

# Lactate dehydrogenase
plot_gene_expression_faceted("LOC102928417", "LDHA", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("Ldha", "Ldha", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("Pdk1", "Pdk1", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)

# Classic Genes
plot_gene_expression_faceted("Hif1a", "Hif1a", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("Vegfa", "Vegfa", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("Epas1", "Epas1", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("Flt1", "Flt1", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)

# Testing
plot_gene_expression_faceted("Ubxn1", "Ubxn1", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("LOC102923832", "Ctsq", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("Tgfbi", "Tgfbi", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("LOC102927525", "Cox7a2", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("Nfatc2", "Nfatc2", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("Lypd6b", "Lypd6b", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("Mettl24", "Mettl24", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("Got1", "Got1", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("Ppp4r2", "Ppp4r2", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("Bcat2", "Bcat2", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("Ucp2", "Ucp2", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("Pgd", "Pgd", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("Rab18", "Rab18", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("Slc16a3", "Slc16a3", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
 

# Mito genes

plot_gene_expression_faceted("Etfa", "Etfa", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("Etfb", "Etfb", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("Etfdh", "Etfdh", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("Acacb", "Acacb", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("Ehhadh", "Ehhadh", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("Abcd4", "Abcd4", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("LOC102903062", "CS", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("Cs", "Cs", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("Acly", "Acly", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("Ampd2", "Ampd2", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
plot_gene_expression_faceted("Gpd1", "Gpd1", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)


# For saving:
 save_faceted_plot("Etfa", "Etfa", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
 save_faceted_plot("Etfb", "Etfb", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
 save_faceted_plot("LOC102923285", "PKM-like", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
 save_faceted_plot("LOC102928417", "LDHA", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
 save_faceted_plot("Tgfbi", "Tgfbi", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)
 save_faceted_plot("Slc1a2", "Slc1a2", EP_vobjdream, LP_vobjDream, EP_Sample_Info, LP_Sample_Info)

