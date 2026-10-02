# Glioblastoma Bulk RNA-seq Analysis

## Overview

This project analyzes publicly available bulk RNA-seq data from glioblastoma and normal brain samples.

The goal was to identify genes that are differentially expressed in glioblastoma and understand the biological processes associated with those changes.

## Dataset

The analysis uses GEO dataset GSE147352.

For this project, I used:

- 85 glioblastoma samples
- 15 normal brain samples

The lower-grade glioma samples in the dataset were not included in this comparison.

## Workflow

The analysis includes the following steps:

1. Download gene count files and sample metadata from GEO
2. Build a gene expression count matrix
3. Filter genes with very low expression
4. Run differential expression analysis using DESeq2
5. Generate PCA, MA plot, volcano plot, and heatmap
6. Separate upregulated and downregulated genes
7. Perform Gene Ontology enrichment analysis for each group

## Differential Expression Analysis

The original count matrix contained 58,302 genes.

After filtering low-expression genes, 25,081 genes were included in the DESeq2 analysis.

Genes were considered significant when:

- adjusted p-value < 0.05
- absolute log2 fold change >= 1

A total of **10,238 differentially expressed genes** were identified.

## Main Results

The upregulated genes were mainly enriched in immune-related processes, including:

- leukocyte cell-cell adhesion
- leukocyte proliferation
- regulation of T cell activation
- leukocyte migration

The downregulated genes were mainly enriched in neuronal and synaptic processes, including:

- regulation of membrane potential
- vesicle-mediated transport in synapse
- neurotransmitter secretion
- synaptic vesicle cycle
- regulation of synaptic plasticity
- action potential

These results show strong differences between glioblastoma and normal brain tissue, including changes in immune activity and loss of normal neuronal expression patterns.

## Project Structure

```text
glioblastoma_bulk_rnaseq/
├── scripts/
│   ├── 00_download_geo_counts.sh
│   ├── 01_prepare_counts.R
│   ├── 02_deseq2_analysis.R
│   └── 03_go_enrichment.R
│
├── results/
│   ├── deseq2_all_genes.csv
│   ├── deseq2_significant_genes.csv
│   ├── deseq2_normalized_counts.csv
│   ├── pca_gbm_vs_normal.pdf
│   ├── ma_plot.pdf
│   ├── volcano_plot.pdf
│   ├── top40_deg_heatmap.pdf
│   ├── GO_upregulated.csv
│   ├── GO_downregulated.csv
│   ├── GO_upregulated_dotplot.pdf
│   └── GO_downregulated_dotplot.pdf
│
├── raw_pipeline/
├── environment.yml
└── README.md
```

## How to Run

Create the conda environment:

```bash
conda env create -f environment.yml
conda activate gbm-rnaseq
```

Download the GEO count files:

```bash
bash scripts/00_download_geo_counts.sh
```

Prepare the count matrix:

```bash
Rscript scripts/01_prepare_counts.R
```

Run differential expression analysis:

```bash
Rscript scripts/02_deseq2_analysis.R
```

Run GO enrichment analysis:

```bash
Rscript scripts/03_go_enrichment.R
```

## Tools Used

- R
- DESeq2
- GEOquery
- clusterProfiler
- org.Hs.eg.db
- ggplot2
- pheatmap
- Linux shell
- Conda

## Notes

This project was independently implemented using publicly available data to demonstrate a reproducible bulk RNA-seq analysis workflow.

Because this is bulk tissue RNA-seq, differences between glioblastoma and normal brain can reflect both gene expression changes and differences in the cell types present in the tissue.
