# cfDNA Precision Diagnostics

A bioinformatics framework for cell-free DNA molecular quantification, prenatal genomic screening, and single-gene genotype inference using targeted next-generation sequencing.

The framework integrates human GRCh38 reference sequences, ClinVar variant annotations, molecular consensus algorithms, statistical modeling, automated sequencing workflows, and laboratory quality monitoring.

## Scientific Workflows

### [cfDNA Molecular Quantification](molecular_quantification/)

UMI-based molecular quantification, PCR duplicate correction, consensus allele identification, and variant allele frequency estimation.

### [Prenatal Genomic Screening](prenatal_genomic_screening/)

SNP-based fetal fraction estimation, chromosome-level fragment quantification, GC-bias correction, and trisomy 13, 18, and 21 screening.

### [ClinVar-Integrated Single-Gene Screening](single_gene_screening/)

Pathogenic variant annotation, targeted molecular sequencing, and Bayesian fetal genotype inference across HBB, CFTR, HEXA, and GJB2.

### [Laboratory QC and Assay Monitoring](laboratory_qc/)

Interactive quality monitoring dashboard integrating sequencing metrics, molecular performance, automated quality alerts, and screening outcomes.

## Benchmark Summary

| Metric | Result |
|---|---:|
| Simulated sequencing samples | 222 |
| Original DNA molecules | 7,145,000 |
| Sequencing reads | 23,224,150 |
| Molecular quantification VAF error reduction | 93.36% |
| Fetal fraction estimation MAE | 0.1105 percentage points |
| Aneuploidy fragment recovery | 99.99983% |
| ClinVar records processed | 4,621,342 |
| Pathogenic and likely pathogenic SNVs retained | 1,357 |
| High-confidence pathogenic SNVs | 453 |
| Bayesian genotype evaluations | 1,080 |
| Accuracy among called genotypes | 97.37% |
| Genotype no-call rate | 61.30% |
| Sequencing samples monitored for alignment QC | 215 |
| Controlled QC failures detected | 9/9 |
| Assay monitoring events | 724 |
| Scientific unit tests passed | 34/34 |

## Workflow Architecture

    Paired-end FASTQ + GRCh38 Reference
                    |
          FastQC / BWA-MEM
                    |
         SAMtools Sort and Index
                    |
          Alignment QC / MultiQC
                    |
          Molecular Quantification
                    |
           Statistical Inference
                    |
        Laboratory QC Monitoring

Nextflow DSL2 automates reference indexing, sequencing quality assessment, alignment, alignment statistics, and MultiQC reporting. Downstream molecular and statistical analyses are implemented as separate Python workflows.

## Workflow Reproducibility

Native macOS and Docker-based Linux execution were compared using three synthetic sequencing samples.

| Validation | Result |
|---|---:|
| Native Nextflow tasks completed | 11/11 |
| Docker Nextflow tasks completed | 11/11 |
| Nextflow cache reuse | 100% |
| Sequencing reads compared | 712,916 |
| Matching alignment metric comparisons | 12/12 |
| Matching ordered alignment signatures | 3/3 |

Alignment signatures compare selected primary alignment fields across execution environments.

## Repository Structure

    molecular_quantification/
    prenatal_genomic_screening/
    single_gene_screening/
    laboratory_qc/
    ngs_pipeline/
    workflows/
    containers/
    tests/
    results/

## Technologies

Python, Bash, NumPy, Pandas, SciPy, PySAM, BWA-MEM, SAMtools,
BCFtools, BEDTools, FastQC, MultiQC, Nextflow, Docker,
Streamlit, Plotly, Matplotlib, and Pytest.

## Docker Environment

Build the bioinformatics container:

    docker build -f containers/Dockerfile -t cfdna-pipeline:1.0.0 .

Run scientific unit tests:

    docker run --rm -v "$PWD:/workspace" -w /workspace \
      cfdna-pipeline:1.0.0 python -m pytest tests/ -v

## Nextflow Execution

After preparing the reference and paired-end sequencing sample sheet:

    nextflow run workflows/main.nf \
      -c workflows/nextflow.config \
      -profile docker \
      --samplesheet workflows/samplesheet.csv \
      --outdir results/nextflow/docker_test

## Laboratory QC Dashboard

After generating upstream benchmark results:

    python laboratory_qc/scripts/aggregate_qc.py
    python laboratory_qc/scripts/assay_monitoring.py
    python -m streamlit run laboratory_qc/dashboard.py

## Data Sources

- [NCBI GRCh38](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/)
- [NCBI ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/)
- [Ensembl](https://www.ensembl.org/)
- Simulated paired-end cfDNA sequencing datasets

## Research Scope

All sequencing experiments use synthetic data generated from human reference sequences and publicly available variant annotations.

Reported performance reflects computational benchmarks rather than clinical diagnostic validation. Large sequencing datasets and downloaded annotation files are excluded from version control.
