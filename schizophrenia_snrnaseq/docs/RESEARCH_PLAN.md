# Research protocol — schizophrenia-associated cell-specific changes

## Scope and hypothesis

**H0:** After controlling for available donor covariates, no tested gene is associated with schizophrenia within a given broad cell class.

**H1:** One or more genes show adjusted disease-associated expression in at least one class. This is an association, NOT a demonstration of causal disease mechanisms.

## Pre-analysis decisions

- Restrict diagnosis to schizophrenia versus neurologically/psychiatrically healthy controls. Exclude other psychiatric diagnoses.
- Preserve donor IDs and select the same number of independent biological units across both groups where feasible. Donor ID must not be confused with per-nucleus barcode.
- Start with pilot (`--max-nuclei-per-donor 150`) to validate code, then scale after reviewing memory and group counts.
- Use only integer raw UMI counts for donor×celltype aggregation, not integrated or log-normalized RNA assays.
- Identify broad annotated cell classes and report authors' annotation provenance. Do not invent cell-type annotation labels.
- Covariates: age, sex, PMI where available and estimable; evaluate batch/RIN and possible diagnosis imbalance before final interpretation.

## Primary analysis

1. Inspect dataset metadata and raw count layers.
2. QC: track donors, genes, nuclei per donor, cell-type balance, distribution of libraries; review source quality-control policy.
3. Aggregate counts into independent donor × cell-type pseudobulk matrix.
4. For each class with 3+ independent donors in both groups, fit negative binomial DESeq2 gene tests with FDR BH.
5. Compare effect size distributions, rank reproducible gene programs, and test optional external pathway gene sets using the tested-gene universe.

## Sensitivity analysis suggestions (not yet implemented)

- Rerun with different nucleus caps / all nuclei to assess dependence on downsampling.
- Add RIN and batch where available; donor-level adjustment may be needed.
- Check effects of sample-level library sizes, PMI, sex balance, compositional differences, and influential donors.
- Evaluate multiple testing across all cell types as a secondary correction in addition to within-class FDR.
- External replication in a second schizophrenia brain cohort if access permits.

## Required results checklist

- [ ] GEO provenance and dataset version logged
- [ ] Donor/diagnosis/cell-type columns verified
- [ ] Raw UMI count source verified
- [ ] At least 3 donors in both groups for each analyzed cell class
- [ ] Pilot versus full dataset clearly labeled
- [ ] Regression formula recorded
- [ ] Number of genes tested and discoveries reported per class
- [ ] Volcano plot and summary table reviewed
- [ ] Relevant confounders and uncertainties documented
- [ ] Interpretation does not claim causation or patient utility
