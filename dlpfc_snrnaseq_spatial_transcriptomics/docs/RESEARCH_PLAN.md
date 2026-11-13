# Research protocol — integrating snRNA-seq and spatial transcriptomics

## Independently developed research question

Can broad neuronal and glial gene-expression signatures built from snRNA-seq data predict spatial cortical patterns in 10x Visium tissue sections? Are inferred cell-type patterns spatially autocorrelated within each section? Are model predictions consistent for genes withheld during fitting?

## Pre-analysis protocol

1. Access official LIBD preprocessed snRNA-seq (`spatialDLPFC_snRNAseq`) and 3-section Visium example subset (`spatialDLPFC_Visium_example_subset`).
2. Verify gene symbol harmonization and uniqueness; inspect cell-type annotations and sample IDs in the R object rather than assume filenames.
3. Construct a cell type reference from raw UMI counts with annotation-preserving downsampling; normalize to CPM.
4. Fit spot-level nonnegative least squares, normalizing fitted mixture coefficients as relative proportions.
5. Calculate kNN Moran I separately for each tissue section and obtain two-sided Monte Carlo permutation p-values; BH correct across tested section/cell-type combinations.
6. Withhold 20% of genes from fitting; evaluate reconstruction of held-out gene expression in each spot using Spearman correlation.

## Interpretation

- Spatial autocorrelation is evidence of a nonrandom geographic pattern, NOT proof of cell-cell interaction or neuronal circuit connectivity.
- Inferred proportions are affected by marker selection, batch, transcript capture efficiency, collinearity among glial/neuron signatures, and Visium resolution.
- Cortical layer patterns (when layer labels are available) should be described as exploratory unless independently replicated and validated histologically.
- Neurotypical DLPFC data are relevant to psychiatric disease **mechanisms and circuit organization**, but contain no case/control disease comparison.

## Extensions for a stronger staff-scientist submission

- Replace basic NNLS with SPOTlight/cell2location and compare performance on held-out genes.
- Add spatial domain segmentation or layered mixed-effects analysis using donor/section as hierarchy.
- Integrate synaptic gene programs and public psychiatric GWAS risk lists with provenance and gene-length/coverage bias controls.
- Validate selected spatial patterns with orthogonal RNA in situ hybridization if available.

## Checklist before GitHub publication

- [ ] Data and gene symbol intersection documented
- [ ] Real versus synthetic runs unambiguously separated
- [ ] Per-section results not treated as independent donors
- [ ] Reproducible fixed seeds and package versions recorded
- [ ] NNLS limitations explained
- [ ] Hold-out gene reconstruction and null/permutation results reported
- [ ] Figures visually reviewed with real tissue coordinates
- [ ] Original dataset generation attributed to LIBD
