# Ovarian Cancer Spatial Transcriptomics with Seurat

## Overview

This project analyzes a public 10x Visium spatial transcriptomics sample from human ovarian cancer using Seurat.

The goal is to identify spatially distinct transcriptional regions in the tumor and examine epithelial, stromal, immune, and vascular gene-expression patterns across the tissue.

## Dataset

The project uses **GSM6506115** from **GSE211956**, a public ovarian tumor Visium dataset.

Sample:
- GSM6506115
- Ovarian tumour
- Human
- Processed with 10x Genomics Space Ranger
- Genome build: GRCh38

The processed count matrix and spatial image files are downloaded directly from NCBI GEO.

## Workflow

1. Download the processed Visium matrix and spatial files
2. Create a Seurat spatial object
3. Calculate mitochondrial percentage
4. Perform basic spot-level QC
5. Normalize expression
6. Identify highly variable genes
7. Scale the data
8. Run PCA
9. Build a nearest-neighbor graph
10. Perform graph-based clustering
11. Run UMAP
12. Visualize clusters in tissue space
13. Identify cluster marker genes
14. Visualize selected tumor, stromal, immune, and vascular markers
15. Calculate spatially variable genes using Moran's I

## Marker Panel

The script checks several marker groups when they are present in the dataset:

- Epithelial/tumor: EPCAM, KRT8, KRT18, KRT19
- Stromal/fibroblast: COL1A1, COL1A2, DCN, COL3A1
- Immune: PTPRC, CD3D, CD68, LST1
- Endothelial: PECAM1, VWF, KDR

These markers are used for visualization only. Final biological region labels should be assigned after reviewing cluster markers and spatial expression patterns.

## Project Structure

```text
ovarian_cancer_spatial_seurat/
├── scripts/
│   ├── 00_download_data.sh
│   └── 01_spatial_seurat_analysis.R
├── data/
│   └── raw/
├── results/
├── environment.yml
├── .gitignore
└── README.md
```

## How to Run

Create the Conda environment:

```bash
conda env create -f environment.yml
conda activate ovarian-spatial
```

Download the processed GEO data:

```bash
bash scripts/00_download_data.sh
```

Run the Seurat analysis:

```bash
Rscript scripts/01_spatial_seurat_analysis.R
```

## Expected Outputs

The analysis creates:

```text
results/
├── 01_qc_violin.pdf
├── 02_spatial_clusters.pdf
├── 03_umap_clusters.pdf
├── 04_spatial_marker_panel.pdf
├── 05_top_spatial_genes.pdf
├── cluster_markers.csv
├── top10_markers_per_cluster.csv
├── spatially_variable_genes.csv
├── spot_metadata.csv
├── ovarian_spatial_seurat.rds
└── sessionInfo.txt
```

