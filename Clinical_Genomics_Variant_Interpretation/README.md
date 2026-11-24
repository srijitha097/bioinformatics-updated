# Clinical Genomics and Variant Interpretation

Three reproducible research workflows using public variant catalogs, published evidence, and reference sequencing datasets. Phenotype experiments include three synthetic scenarios and 30 published solved cases. The sequencing analyses use public reference or cell-line data, and the CNV caller exercise uses a separate simulated cohort. All reports are for research and education.

## Projects

| Project | Workflow | Main output |
|---|---|---|
| Rare disease | [Phenotype prioritization and ACMG/AMP evidence](rare_disease_variant_curation/README.md) | [GJB2 variant assessments](rare_disease_variant_curation/results/variant_assessments.md) |
| Somatic oncology | [CIViC assertions and AMP/ASCO/CAP evidence](oncology_somatic_variant_interpretation/README.md) | [Oncology assessments](oncology_somatic_variant_interpretation/results/oncology_assessments.md) |
| Sequencing | [Exome, gene-panel subset, and copy-ratio analysis](exome_gene_panel_cnv_analysis/README.md) | [Sequencing research report](exome_gene_panel_cnv_analysis/results/research_report.md) |

[Run summary](RUN_SUMMARY.md) · [Repository review](FINAL_GITHUB_REVIEW.md)

## Reproduction

Run commands from the repository root. Create the environment with `conda env create -f environment.yml` and activate it with `conda activate clinical-genomics-research`. The simulation and final target-level audit require NumPy. BAM analysis and the independent SNV cross-check use samtools and bcftools 1.25. Figure generation uses Matplotlib 3.11.2.

```bash
python3 -m unittest discover -s tests -v
bash reproduce.sh  # Full replay with archived inputs installed
```

Raw downloads belong in each project’s ignored `data/` directory. Results include source hashes where applicable. Run figure generation after reproducing all project outputs.

## Audited findings

Publication overlap materially affects rare-disease ranking. Clinical evidence breadth explains much of the oncology tier discrepancy. Exact exome SNV matching reproduced the baseline counts, and the CNV simulation showed a reproducible precision–recall trade-off across seeds.

[Final scientific audit](FINAL_AUDIT.md) · [Methods and reproduction commands](comparative_genomics/README.md) · [Data provenance](DATA_PROVENANCE.md) · [Limitations](LIMITATIONS.md)
