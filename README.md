# blca_bcg_resistance
R code for BCG resistance in non-muscle invasive bladder cancer (NMIBC)

## Overview
This repository contains the R script for single-cell RNA sequencing data processing and integration described in our manuscript.
The workflow includes raw matrix extraction, doublet detection, quality control, normalization, HVG selection, PCA, Harmony batch correction, clustering, UMAP visualization and cell type annotation.

## Requirements
R (>=4.2.1)
Main packages: Seurat, Hmisc, survival, ggplot2, tidyverse, scDblFinder, org.Hs.eg.db

## Analysis workflow
1. Integrate public GSE269877 NMIBC scRNA-seq data and four in-house NMIBC samples
2. Doublet identification with scDblFinder and QC filtering
3. Normalization, HVG selection, PCA and Harmony batch correction
4. Graph-based clustering, UMAP visualization and cell annotation
5. Sub-clustering for refined cellular subtype classification

## Notes
- **No raw clinical data, raw count matrices or Seurat objects (.rds/.rda) are deposited here**, due to data copyright and patient privacy restrictions.
- This code reproduces the single-cell analysis pipeline in our manuscript.
