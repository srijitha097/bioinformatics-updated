# Resume bullets — do not claim until real LIBD analysis is complete

**Human Prefrontal Cortex snRNA-seq and Spatial Transcriptomics Integration | R, Python, SpatialExperiment, 10x Visium**

- Built a reproducible cross-modality workflow deriving human DLPFC neuronal and glial transcriptional reference profiles from snRNA-seq and estimating cell-type composition in spatial transcriptomics sections.
- Quantified spatial cell-type organization with per-section Moran's I and permutation testing, and assessed deconvolution consistency using held-out gene-expression reconstruction.

**Use these only after `bash run_real.sh` succeeds and you have verified the actual output and biological interpretation.** Add only metrics generated and reviewed by your code, not numbers taken from the source study.

## Interview narrative

Harmonized gene symbols → reference expression estimates → nonnegative mixture model → spatial autocorrelation → held-out gene validation → nested section/donor structure → technical limitations → next experiment.
