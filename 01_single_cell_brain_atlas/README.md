# Single-Cell Brain Atlas
**Author:** Srijitha Gangi  
**M.S. Bioinformatics**

## Overview

This project analyzes publicly available single-cell RNA-seq data from embryonic mouse brain using Seurat.

The goal is to perform quality control, preprocessing, dimensionality reduction, clustering, marker-gene analysis, and first-pass cell-type annotation.

## Dataset

The analysis uses a public 10x Genomics embryonic mouse brain dataset.

The filtered feature-barcode matrix is used as the input for the Seurat workflow.

## Workflow

1. Read the 10x gene expression matrix
2. Create a Seurat object
3. Calculate mitochondrial percentage
4. Perform quality-control filtering
5. Normalize the data
6. Identify highly variable genes
7. Scale the data
8. Run PCA
9. Build the neighbor graph
10. Cluster the cells
11. Run UMAP
12. Identify cluster marker genes
13. Review canonical brain markers
14. Assign first-pass cell-type labels based on cluster markers

## Analysis

The workflow uses:

- Seurat
- LogNormalize
- 2,000 highly variable genes
- PCA
- graph-based clustering
- UMAP
- FindAllMarkers for cluster marker identification

The final analysis produced 16 clusters.

## Cell-Type Annotation

Cell types were assigned manually based on cluster-specific marker genes.

The annotations include:

- excitatory neurons
- inhibitory interneurons
- neural progenitor populations
- Cajal-Retzius-like neurons
- oligodendrocyte precursor cells
- microglia

## Main Results

The analysis identified 16 transcriptionally distinct clusters.

Examples of marker-supported annotations include:

- SST interneurons based on markers such as `Sst` and `Lhx6`
- HTR3A-positive interneurons based on `Htr3a`
- intermediate neural progenitors based on `Eomes` and `Neurog2`
- oligodendrocyte precursor cells based on `Sox10`, `Pdgfra`, `Gpr17`, and `Bcas1`
- microglia based on `Adgre1`, `Trem2`, `Fcrls`, and `Siglech`

## Project Structure

```text
01_scRNA_brain_atlas/
├── scripts/
│   ├── 01_seurat_brain_atlas.R
│   └── 02_annotate_celltypes.R
│
├── results/
│   ├── 01_qc_violin.pdf
│   ├── 02_qc_scatter.pdf
│   ├── 03_pca.pdf
│   ├── 04_elbow.pdf
│   ├── 05_umap_clusters.pdf
│   ├── 06_marker_dotplot.pdf
│   ├── 07_umap_celltypes.pdf
│   ├── cluster_markers.csv
│   ├── top10_markers_per_cluster.csv
│   └── sessionInfo.txt
│
├── environment.yml
├── .gitignore
└── README.md
```

## How to Run

Create the conda environment:

```bash
conda env create -f environment.yml
conda activate scrna-brain
```

Place the filtered 10x feature-barcode matrix at:

```text
data/raw/filtered_feature_bc_matrix/
```

Run the main Seurat analysis:

```bash
Rscript scripts/01_seurat_brain_atlas.R
```

Run the cell-type annotation step:

```bash
Rscript scripts/02_annotate_celltypes.R
```

## Outputs

The main output files are:

- `01_qc_violin.pdf` - quality-control distributions
- `02_qc_scatter.pdf` - QC relationships
- `03_pca.pdf` - PCA visualization
- `04_elbow.pdf` - principal component selection
- `05_umap_clusters.pdf` - UMAP with cluster numbers
- `06_marker_dotplot.pdf` - expression of selected marker genes
- `07_umap_celltypes.pdf` - UMAP with cell-type labels
- `cluster_markers.csv` - marker genes for all clusters
- `top10_markers_per_cluster.csv` - top marker genes for each cluster
- `sessionInfo.txt` - R package and session information

## Tools Used

- R
- Seurat
- ggplot2
- dplyr
- patchwork
- Conda
