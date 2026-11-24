# Verified run summary

Run date 24 November 2026. Scope is research and education.

## 1. Rare disease and phenotype prioritization

| Metric | Result |
|---|---:|
| Eligible expert-reviewed or practice-guideline GRCh38 SNVs across the 20-gene set | 6,739 |
| Candidate–scenario rows after per-gene sampling | 200 |
| Unique candidate variants across three synthetic scenarios | 130 |
| GJB2 variant evidence reconstructions | 2 |

CFTR matched both CF-like HPO terms. The single-term cardiac and hearing-loss scenarios produced several full-score gene ties, including genes linked to occasional manifestations. This shows the limits of unweighted overlap in a restricted panel. The queue uses real ClinVar and HPO records, and the scenarios are synthetic. The evidence reports reconstruct published ClinGen classifications, so there is no independent agreement estimate. OMIM identifiers were sourced through HPO and GeneReviews because the direct OMIM page was inaccessible.

[Ranking table](rare_disease_variant_curation/results/phenotype_gene_scores.tsv) · [Research report](rare_disease_variant_curation/results/variant_assessments.md)

**Resume bullets**

- Built and executed a disease-scoped HPO prioritization workflow that screened 6,739 ClinVar SNVs and generated 200 ranked candidate–scenario pairs.
- Reconstructed two ClinGen ACMG/AMP GJB2 assessments and quantified ancestry-specific gnomAD exome and genome allele frequencies to document reduced-penetrance interpretation.

## 2. Somatic oncology and AMP tiering

| Metric | Result |
|---|---:|
| Accepted CIViC assertions with an AMP tier | 124 |
| Distinct disease labels | 85 |
| Tier I, Level A | 120 |
| Tier II, Level C | 4 |
| Literature reviews with verified CIViC evidence links | 6 |

The dated source was the 1 November 2026 accepted-only CIViC export. Four reviewed assertions were Tier I/Level A and two were Tier II/Level C. The audits verified linked CIViC evidence and PubMed IDs. The purposive sample and selection of already accepted assertions favor established evidence, and the source release excludes or underrepresents lower tiers. The six reviews are unblinded research assessments.

[Tier table](oncology_somatic_variant_interpretation/results/tier_by_assertion_type.tsv) · [Oncology report](oncology_somatic_variant_interpretation/results/oncology_assessments.md)

**Resume bullets**

- Adapted a CIViC parser to a verified accepted-only release and summarized 124 AMP-tiered assertions across 85 disease labels.
- Reviewed six tumor-specific predictive and diagnostic assertions and implemented evidence provenance checks linking CIViC evidence IDs to primary PubMed citations.

## 3. Exome, targeted analysis and CNV

| Run | Verified result |
|---|---|
| Real NA12878 exome, exact autosomal SNV genotype benchmark | 18,435 TP, 200 FP and 1,573 FN. Precision 98.93%, recall 92.14% |
| Real 10-gene WES-derived subset | 17 TP, 0 FP and 1 FN. Precision 100%, recall 94.44% across 18 truth SNVs |
| Real HCC1143 exome coverage | 18,986 target bases. Mean depth 25.23× normal and 28.52× tumor |
| Real HCC1143 exploratory copy ratios | 79 of 93 targets passed the two-sample depth floor |
| Synthetic CNV recovery | 36 TP, 30 FP and 10 FN. Target-level precision 54.55%, recall 78.26% |

Independent bcftools filtering and genotype intersection reproduced the Python SNV counts. The exome benchmark uses published calls from one genome. The panel is a restriction of that same run. The HCC1143 ratios are exploratory on a downsampled chromosome slice without independent CNV truth, and the CNV caller metrics come from a separate simulation. Both simulated two-target events were missed at their truth locations. No clinical diagnostic validation was performed.

[Sequencing report](exome_gene_panel_cnv_analysis/results/research_report.md) · [SNV cross-check](exome_gene_panel_cnv_analysis/results/benchmark_crosscheck.tsv)

**Resume bullets**

- Benchmarked public NA12878 exome calls against GIAB in 31.56 Mb of confident capture and a ten-gene subset, independently confirming exact genotype counts with bcftools.
- Quantified matched tumor–normal depth and relative copy ratios from 93 public exome targets, and measured target-level gain/loss recovery in a separate held-out simulation.

## Verification

All eight automated tests passed. Input hashes, evidence links, output tables and five PNG figures are retained in the project results folders. See [PROGRESS.md](PROGRESS.md) for the failures, fixes and reruns.

## Phase two · Working research results

The preceding sections report the preserved phase-one experiments. The follow-up used new cohorts and methods where documented.

| Project | Baseline | Phase-two result |
|---|---|---|
| Rare disease | Three synthetic scenarios and two evidence reconstructions | 30 published cases and 19 evidence reconstructions across five categories. Phenotype-only Top-1 was 95.0% in the selected candidate panel and 45.5% under the explicit distinct-publication restriction. |
| Oncology | Six unblinded literature reconstructions | 23 additional contexts across 22 molecular profiles. Initial literature review assigned 18 assessments and deferred five. Four assigned assessments matched the reference broad tier. Post-comparison regulatory evidence is documented separately. |
| Exome | SNV precision 98.93% and recall 92.14%, ten-gene subset | SNV counts reproduced. Exact indel precision 83.95%, recall 70.74%, F1 76.78%. Expanded 30-gene subset had 34 exact SNV matches and two misses. |
| CNV | Synthetic precision 54.55% and recall 78.26% | New simulation with GC/batch effects. On held-out data the fixed comparator had 11.03% precision and 63.04% recall, while the development-selected configuration had 85.29% precision and 63.04% recall. The original simulation and this simulation differ. Real array comparison covered 93 exome targets but lacked adequate truth-event diversity. |

[Phase-two experiments](comparative_genomics/README.md) · [Validation](VALIDATION_REPORT.md) · [Data provenance](DATA_PROVENANCE.md) · [Limitations](LIMITATIONS.md) · [Traceable resume bullets](RESUME_BULLETS.md)

## Final audit

All 15 tests passed and all eight plots reproduced. The citation audit found explicit case-PMID overlap in 21 of 30 cases. The original CNV holdout shared 113 true-positive targets between configurations, with 32 unique to each. [Final audit](FINAL_AUDIT.md)

The final saved-input replay and the figure-by-figure checks are recorded in the [repository review](FINAL_GITHUB_REVIEW.md).
