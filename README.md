# ThermTol WGCNA analysis

This repository contains scripts to create gene co-expression networks from normalized
gene counts (processed with DESeq2 [1]) using WGCNA [2], pull enriched GO pathways in
modules correlated with a trait of interest, and permute null distributions for the
expected number of hub genes and genes with transgressive F1 expression in all modules. 
The dataset that this pipeline was developed for was tissue-level RNAseq of individuals
from two parent species and their F1 hybrids, exposed to a control and a high temperature
condition. Scripts were written to be run interactively in R. 

Used to perform thermal tolerance co-expression analysis in Payne et al 2021. ThermTol
expression data is available upon request.

[1] Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15:550. 10.1186/s13059-014-0550-8
[2] Langfelder P, Horvath S (2008) WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 2008, 9:559. https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559

## Dependencies

R-3.6.1

The following packages are required for each script:
ThermTol_WGCNA        ~ DESeq2, WGCNA, tximportData, tximport, GenomicFeatures, readr, rhdf5
GO_analysis_WGCNA     ~ GOstats, GSEABase, biomaRt
WGCNA-hub-permutation ~ dplyr

## Contact

For questions and/or comments, please email me at cypayne at stanford dot edu.
