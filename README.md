### Read me associated with: 
# Early metabolic reprogramming in the placenta shapes resilience to gestational hypoxia in high-elevation deer mice

# Pipeline Scripts
## Pipeline for Annotation Generation
### Input:
- HiFi fastq files available on NCBI: [link]

### Pipeline scripts:
- `IsoQuant.sh`

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
- `EP_IXN_ReactionNorm.R` - graphs reaction norms for interaction genes, Table S6

### Output:
- Differential expression counts tables

## Summary files of differential expression analyses
### Input
- Differential expression counts tables

### Scripts found in `RNA_Seq_RScripts`:
- `EP_LP_DreamOutput_Summarize.R`

### Output:
- EP_ISO_Ortho_Summary.xlsx
- LP_ISO_Ortho_Summary.xlsx
- LP_JZ_ISO_Ortho_Summary.xlsx
- Combination of these files is `Dataset_S1.xlsx`, which is available on [FigShare](https://doi.org/10.6084/m9.figshare.32556435)

# Figures
## Figure 1 - Fetal Vasculature
### Input:
- Placenta images available on Dryad [link]
- Images quantified in FIJI using macro (LamCyto_Analysis.ijm)
- Output of quantification: EP_BW_ME_Quantification.xlsx

### Scripts found in `Placental_Histology`
- `ProgenitorQuant_Plot.R` - generate plot
- `ProgenitorQuant.R`- plot statistics, Table S5

## Figure 2 - Hypoxia heatmap & respresentative gene plots
### Input:
- MetaData_EP.xlsx
- MetaData_LP.xlsx
- EP_Pman_ExtMMFrac_readcounts_Exon.xlsx
- LP_Pman_ExtMMFrac_readcounts_Exon.xlsx
- EP_ISO_Ortho_Summary.xlsx

### Scripts found in `RNA_Seq_RScripts`:
- `EP_ZScore.R` - generate heatmap
- `EP_LP_Gene_Plots.R` - generate representative gene plots

## Figure 3 - GSEA
### Input:
- EP_ISO_Ortho_Summary.xlsx

### Scripts found in `GSEA_RScripts`
- `EP_GSEA_Plot.R` - generate plot, Table S7
- `EP_GSEA_PBS.R`, - determine leading edge genes under selection, Table S8

## Figure 4 - Population differences across gestation
### Input:
- EP_ISO_Ortho_Summary.xlsx
- LP_ISO_Ortho_Summary.xlsx
- MetaData_EP.xlsx
- MetaData_LP.xlsx
- EP_Pman_ExtMMFrac_readcounts_Exon.xlsx
- LP_Pman_ExtMMFrac_readcounts_Exon.xlsx

### Scripts found in RNA_Seq_RScripts
- `SharedStrain_Plot.R` - generate plot
- `EP_LP_Gene_Plots.R` - generate representative gene plots

## Figure 5 - WGCNA_RScripts
### Input:
- MetaData_EP.xlsx
- MetaData_LP.xlsx
- EP_Pman_ExtMMFrac_readcounts_Exon.xlsx
- LP_Pman_ExtMMFrac_readcounts_Exon.xlsx

### Scripts found in WGCNA_RScripts
- `EP_WGCNA_PowerTest.R` - testing multiple softPower values to achieve scale-free topology
- `EP_WGCNA_BWME.R` - generate network for lowlanders (BW) and highlanders (ME) separately, Table S9, S10
- `EP_WGCNAnetPres.R` - determine network preservation between populations in early pregnancy, generates plot
- `EP_WGCNA_GOEnrich.R` - run GO analysis on each lowland module, Table S12, S13, S14, S15
- `EP_Fishers_Nest.R` - Fishers exact test to determine enrichment/depletion of DE genes within each module, Table S11, generates plot





	
