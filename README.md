# scRNA-seq Seurat Pipeline

A reproducible **single-cell RNA-seq analysis pipeline** built with **R and Seurat** for automated preprocessing, clustering, and marker discovery.

This project implements an end-to-end workflow for analyzing scRNA-seq datasets including quality control, normalization, dimensionality reduction, clustering, and identification of cluster-specific marker genes.

---

## Features

- Supports **Seurat `.rds` objects** or **10x Genomics datasets**
- Automated **quality control filtering**
- **SCTransform** or **LogNormalize** normalization
- **PCA → UMAP → graph-based clustering**
- **Cluster marker gene discovery**
- Automatic export of **figures, tables, and processed Seurat objects**

---

## Pipeline Steps

1. Load scRNA-seq dataset  
2. Compute QC metrics (`nFeature_RNA`, `nCount_RNA`, `percent.mt`)  
3. Filter low-quality cells  
4. Normalize gene expression data  
5. Perform PCA and construct neighbor graph  
6. Generate **UMAP embedding**  
7. Run **unsupervised clustering**  
8. Identify **cluster-specific marker genes**

---

## Output

Results are automatically saved to the following directories.

### Figures

outputs/figures

Examples include:

- QC violin plots  
- Feature QC plots  
- UMAP cluster visualizations  

### Tables

outputs/tables

Examples include:

- markers_per_cluster.csv

### Processed Data

data/processed

Examples include:

- seurat_qc.rds  
- seurat_clustered.rds  

---

## Requirements

Required R packages:

- Seurat  
- ggplot2  
- dplyr  
- readr  

---

## Run the Pipeline

Run the script in R:

source("scRNAseq_seurat_pipeline.R")

---

## Author

Bryan Gomez
Computer Science (Bioinformatics)