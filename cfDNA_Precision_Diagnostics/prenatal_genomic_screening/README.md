# Prenatal Genomic Screening

Computational analysis of cell-free DNA sequencing data for fetal fraction estimation and chromosomal aneuploidy screening.

## Fetal Fraction Estimation

Informative biallelic SNPs were selected from Ensembl genomic variation annotations using reference allele consistency, local sequence uniqueness, GC content, and genomic spacing.

Fetal fraction was estimated from BAM-derived allele counts using maximum likelihood estimation and likelihood-based confidence intervals.

### Benchmark Results

| Metric | Result |
|---|---:|
| Informative SNPs | 30 |
| Sequencing samples | 6 |
| Original DNA fragments | 90,000 |
| Sequencing reads | 180,000 |
| Fetal fraction MAE | 0.1105 percentage points |
| Fetal fraction RMSE | 0.1262 percentage points |
| Confidence interval coverage | 5/5 evaluation samples |

## Chromosomal Aneuploidy Screening

Targeted sequencing windows were selected across chromosomes 1, 2, 3, 13, 18, and 21.

The workflow evaluates chromosomal dosage abnormalities using fragment quantification, GC-content normalization, euploid baseline calibration, and chromosome-specific Z-scores.

### Methodology

- Targeted GRCh38 reference preparation.
- Genomic window selection and sequence quality assessment.
- Simulation of euploid and trisomy sequencing samples.
- BWA-MEM alignment and SAMtools processing.
- Fragment-level counting across genomic windows.
- GC-bias correction using regression-based normalization.
- Euploid baseline calibration.
- Z-score classification and no-call assessment.

### Sequencing Benchmark

| Metric | Result |
|---|---:|
| Selected genomic windows | 180 |
| Euploid calibration samples | 20 |
| Evaluation samples | 60 |
| Original DNA fragments | 4,800,000 |
| Sequencing reads | 9,600,000 |
| Fragment recovery | 99.99983% |
| Alignment QC pass rate | 100% |

### Screening Performance

| Metric | Raw | GC Corrected |
|---|---:|---:|
| Accuracy among called samples | 50.00% | 66.67% |
| Trisomy detection among called affected samples | 33.33% | 55.56% |
| Specificity among called euploid samples | 100% | 100% |
| No-call rate | 20% | 20% |

### GC Correction

| Chromosome | Raw Baseline CV | Corrected Baseline CV |
|---|---:|---:|
| 13 | 1.8656% | 1.4093% |
| 18 | 2.5070% | 1.4572% |
| 21 | 2.4296% | 1.5463% |

GC normalization reduced chromosome-specific baseline variability by approximately 24-42%.

## Implementation

Analysis scripts are available in:

    ngs_pipeline/scripts/select_prenatal_snps.py
    ngs_pipeline/scripts/simulate_prenatal.py
    ngs_pipeline/scripts/estimate_prenatal_ff.py
    ngs_pipeline/scripts/select_aneuploidy_windows.py
    ngs_pipeline/scripts/simulate_aneuploidy.py
    ngs_pipeline/scripts/count_aneuploidy_fragments.py
    ngs_pipeline/scripts/analyze_aneuploidy.py

Benchmark reports and visualizations:

    results/ngs/fetal_fraction/
    results/ngs/aneuploidy/benchmark/

## Technologies

Python, SciPy, NumPy, Pandas, PySAM, BWA-MEM, SAMtools,
FastQC, MultiQC, Matplotlib, and Ensembl.

## Research Scope

Synthetic targeted sequencing experiments using simulated fetal fraction and chromosome dosage information. The results do not establish clinical screening performance.
