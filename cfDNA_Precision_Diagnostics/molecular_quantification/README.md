# cfDNA Molecular Quantification

Molecular quantification of cell-free DNA using unique molecular identifiers, directional error correction, and consensus-based variant allele frequency estimation.

## Methodology

- Simulation of DNA fragments, PCR amplification, sequencing errors, and UMI errors.
- Paired-end sequencing alignment against a targeted human GRCh38 reference.
- Extraction of molecular identifiers and fragment coordinates.
- Grouping of sequencing reads into molecular families.
- Directional UMI error correction using Hamming distance.
- Consensus allele identification using family size and agreement thresholds.
- Comparison of raw and consensus allele frequencies against simulated ground truth.

## Benchmark Results

| Metric | Result |
|---|---:|
| Original DNA molecules | 5,000 |
| Sequencing reads | 43,818 |
| Mapping rate | 100% |
| Properly paired reads | 100% |
| Mean mapping quality | 60 |
| True variant allele frequency | 1.8600% |
| Raw variant allele frequency | 1.9973% |
| Consensus variant allele frequency | 1.8509% |
| VAF error reduction | 93.36% |
| Accepted consensus molecular families | 3,836 |
| Exact molecular recovery | 76.72% |
| Allele accuracy among matched molecules | 100% |

## Molecular Family Threshold Analysis

| Minimum Family Size | Accepted Families | Molecular Recovery |
|---|---:|---:|
| 1 | 5,045 | 99.16% |
| 2 | 3,836 | 76.72% |
| 3 | 3,003 | 60.06% |
| 4 | 2,349 | 46.98% |
| 5 | 1,834 | 36.68% |

A minimum family size of two provided the lowest VAF estimation error among the tested thresholds for this benchmark.

Increasing the minimum family size improved consensus stringency but reduced molecular recovery.

## Implementation

Sequencing and molecular quantification scripts are available in:

    ngs_pipeline/scripts/simulate_cfdna.py
    ngs_pipeline/scripts/quantify_molecules.py

Benchmark results are available in:

    results/ngs/molecular/

## Technologies

Python, NumPy, Pandas, PySAM, BWA-MEM, SAMtools, FastQC, MultiQC, and Matplotlib.

## Research Scope

Synthetic sequencing experiments using a targeted human genomic reference. Performance estimates are specific to the simulated sequencing conditions.
