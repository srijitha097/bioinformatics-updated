# Schizophrenia snRNA-seq: Exploratory Results

## Dataset
NCBI GEO GSE254569 — human orbitofrontal cortex.

## Analysis
- 107,953 nuclei analyzed from 69 donors
- 36 schizophrenia donors and 33 controls
- Seven major brain cell populations
- Donor-level pseudobulk differential expression
- DESeq2 adjustment for age, PMI, RIN, sex and library batch

## Results
- 283 significant gene-cell-type associations (FDR < 0.05)
- Astrocytes: 205
- OPCs: 30
- Excitatory neurons: 15
- Inhibitory neurons: 7
- Microglia: 10
- Oligodendrocytes: 13
- Endothelial cells: 3

## Reactome enrichment
Three significant astrocyte-associated pathways:
1. Protein-protein interactions at synapses
2. Neurexins and neuroligins
3. Neuronal System

## Limitations
Exploratory subsampled secondary analysis, not the full cohort.
Neuronal RNA contamination in astrocytes remains unresolved.
Raw UMI provenance and full-dataset sensitivity require
additional verification. Findings are not independently replicated.
