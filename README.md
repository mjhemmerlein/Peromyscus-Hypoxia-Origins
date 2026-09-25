# Read me associated with: Early metabolic reprogramming in the placenta shapes resilience to gestational hypoxia in high-elevation deer mice

# Pipeline Scripts
## Pipeline for Annotation Generation
### Input:
- HiFi fastq files available on NCBI: [link]

### Pipeline scripts:
- IsoQuant.sh

## Pipeline for RNA sequencing:
### Input: 
- Fastq files available on NCBI
- Early pregnancy files: [link]
- Late pregnancy files: [link]

### Pipeline scripts:
- `fastp_loop_EP.sh` / `fastp_loop_LP.sh` - Quality control filtering & trimming
- `hisat_build.sh` - Indexing Peromyscus genome
- `hisat_align_EP.sh` / `hisat_align_LP.sh` - align reads to genome
- `EP_featurecounts_ExtMFrac_Exon.sh` / `LP_featurecounts_ExtMFrac_Exon.sh` / `LP_JZ_featurecounts_ExtMFrac_Exon.sh` - counting genes that align to annotation

### Output:
- EP_Pman_ExtMMFrac_readcounts_Exon.xlsx
- LP_Pman_ExtMMFrac_readcounts_Exon.xlsx (LP = LP labyrinth zone)
- LP_JZ_Pman_ExtMMFrac_readcounts_Exon.xlsx

## Differential Expression Analysis for RNA sequencing
### Inputs found in `RNA_Seq_RawData`:
- EP_Pman_ExtMMFrac_readcounts_Exon.xlsx
- LP_Pman_ExtMMFrac_readcounts_Exon.xlsx (LP = LP labyrinth zone)
- LP_JZ_Pman_ExtMMFrac_readcounts_Exon.xlsx
- MetaData_EP.xlsx
- MetaData_LP.xlsx

### Scripts found in `RNA_Seq_RScripts`:
- `EP_Dream.R` - Early pregnancy DE counts from combined populations, lowland only, highland only
- `LP_Dream.R` - Late pregnancy labyrinth zone DE counts from combined populations, lowland only, highland only
- `LP_JZ_Dream.R` - Late pregnancy junctional zone DE counts from combined populations, lowland only, highland only

### Output:
- Differential expression counts tables found in `RNA_Seq_Output/Dream_RawFiles`
- Combined population files - Strain (Population), O2 (Hypoxia), IXN (interaction)
- Lowland only (BW) - Strain (Population), O2 (Hypoxia), IXN (interaction)
- Highland only (ME) - Strain (Population), O2 (Hypoxia), IXN (interaction)

## Summary files of differential expression analyses
### Input
- Differential expression counts tables found in Dream_RawFiles

### Scripts found in `RNA_Seq_RScripts`:
- `EP_LP_DreamOutput_Summarize.R`

### Output:
- EP_ISO_Ortho_Summary.xlsx
- LP_ISO_Ortho_Summary.xlsx
- LP_JZ_ISO_Ortho_Summary.xlsx
- Combination of these files is Dataset_S1.xlsx

# Figures
## Figure 1 - Fetal Vasculature
### Input:
- Placenta images available on Dryad [link]
- Images quantified in FIJI using macro (LamCyto_Analysis.ijm)
- Output of quantification: EP_BW_ME_Quantification.xlsx

### Scripts found in `Placental_Histology`
- `ProgenitorQuant_Plot.R` - Figure 1 plots
- `ProgenitorQuant.R` - Figure 1 statistics, Table S5

## Figure 2 - Hypoxia heatmap & respresentative gene plots
### Input:
- MetaData_EP.xlsx
- MetaData_LP.xlsx
- EP_Pman_ExtMMFrac_readcounts_Exon.xlsx
- LP_Pman_ExtMMFrac_readcounts_Exon.xlsx
- EP_ISO_Ortho_Summary.xlsx

### Scripts found in `RNA_Seq_RScripts`:
- EP_ZScore.R - Figure 2 heatmap plot
- EP_LP_Gene_Plots.R - Figure 2 representative gene plots

## Figure 3 - GSEA
### Input:
- EP_ISO_Ortho_Summary.xlsx

### Scripts found in `GSEA_RScripts`
- EP_GSEA_Plot.R - Figure 3 plot

## Figure 4 - Population differences across gestation
### Input:
- EP_ISO_Ortho_Summary.xlsx
- LP_ISO_Ortho_Summary.xlsx
- MetaData_EP.xlsx
- MetaData_LP.xlsx
- EP_Pman_ExtMMFrac_readcounts_Exon.xlsx
- LP_Pman_ExtMMFrac_readcounts_Exon.xlsx

### Scripts found in RNA_Seq_RScripts
- SharedStrain_Plot.R - Figure 4 plots
- EP_LP_Gene_Plots.R - Figure 4 representative gene plots

## Figure 4 - WGCNA_RScripts







	
