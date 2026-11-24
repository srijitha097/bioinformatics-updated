# Final scientific and reproducibility audit

All three computational research projects completed their audit. The saved-input rerun reproduced the phase-two output hashes, including all three figures. The five baseline plots were regenerated and passed the baseline checksum test. All 15 automated tests passed. The rerun covered case selection, cached gnomAD aggregation, scoring, citation restrictions, literature selection and evidence joins, reference normalization, benchmarking, BAM subsampling, copy-ratio intersection, simulation and plotting. Evidence judgments were replayed from encoded reviews. They have not received independent clinical adjudication.

## Rare-disease prioritization and curation

The original archive records confirm 30 distinct genes, alleles and case PMIDs, one solved causative interpretation per selected case, GRCh38 coordinates and at least two included HPO terms. All nine combinations of scoring method and perturbation used all 30 cases. The score assigned to a gene is identical for alleles of that gene. Ties contribute a fractional probability of being within k and the average reciprocal of each possible rank.

| Original phenotypes | Expected Top-1 | Top-5 | Top-10 | MRR |
|---|---:|---:|---:|---:|
| Phenotype | 95.00% | 100% | 100% | 0.9690 |
| Frequency | 3.33% | 16.67% | 33.33% | 0.1332 |
| Combined | 93.33% | 100% | 100% | 0.9667 |

The frequency-only null follows analytically from ranking the same 30 candidates exactly once per causal gene. Fourteen lack retrievable allele counts. The independently recomputed original and perturbation means agree with the [rank bounds](comparative_genomics/rare/ranks.tsv) and [summaries](comparative_genomics/rare/recovery.tsv). A full-score result can represent a tie.

In 21 cases the case PMID explicitly supports at least one observed term for a disease associated with the causal gene in the HPO annotation file. Under the different-PMID rule, seven cases retain zero causal-gene matches. Across the same 30 cases, restricted Top-1 is 45.48%, Top-5 is 60.10%, Top-10 is 70.75%, and MRR is 0.5407. This rule also discards annotations without identifiable PubMed support. It cannot isolate a causal effect of publication leakage, and distinct PMIDs need not represent independent patients. [Cases](comparative_genomics/rare/cases.tsv) · [Citation audit](comparative_genomics/rare/distinct_publication.tsv)

The 19 versioned SCVs match their archived variant identity, expert review status and submitted class. Counts are five pathogenic, four likely pathogenic, four uncertain, two likely benign and four benign. The table distinguishes gene-specific criteria, historical ENIGMA likelihood assessments and CFTR2 evidence. The interpretations reconstruct published decisions. No classifier was scored against those decisions. [Evidence table](comparative_genomics/rare/evidence_reviews.tsv)

## Oncology evidence audit

The analysis covered 23 contexts and 22 distinct profiles. Three exact tier/level matches were BRAF V600E and V600K melanoma combination therapy and KIT D816V imatinib resistance. FLT3 ITD in AML matched at Tier I with a B-versus-A level difference. Four broad matches among 18 assigned is 22.2%. Five deferred contexts are outside that denominator. All 14 assigned discrepancies were Tier II assessments against Tier I references.

| Assigned discrepancies | Interpretation |
|---|---|
| 4 therapeutic contexts, AIDs 8, 33, 34, 56 | Selected case, model or early clinical studies supported lower evidence levels in isolation. CIViC records disease-matched or fusion-class regulatory evidence. Omitting that evidence under-tiers a clinical synthesis. The FDA addendum documents the revised evidence scope for 8, 33 and 56. Assertion 34 also has regulatory and NSCLC guideline evidence in the release. |
| 9 diagnostic contexts, AIDs 21, 24, 76, 79, 80, 87, 97, 98, 99 | Discovery and diagnostic series provide limited standalone evidence. Reference narratives incorporate professional criteria, including WHO desirable diagnostic features for several fusions. A rare partner or a single assay can be less specifically supported than a broader defining rearrangement. The initial review omitted the guideline basis. |
| 1 prognostic context, AID 20 | The selected colorectal cohort supports a survival association. The CIViC assertion also cites colon cancer guidelines and six evidence items. Melanoma drug sensitivity supplies no evidence for that colorectal prognostic outcome. |

AIDs 5, 43, 51, 52 and 93 were deferred for missing, indirect or insufficiently profile-specific evidence. AID 43 has a Tier II-C reference. Some deferrals are appropriate for the narrow source and inadequate for a comprehensive tier assignment. No broad-tier arithmetic or join error was found. These were reference-selected contexts and earlier reference exposure prevents a blinded accuracy interpretation. WHO guideline content here was checked through CIViC’s source narrative rather than a full independent WHO review. [Original assessments and cited studies](comparative_genomics/oncology/report.md) · [Reference comparison](comparative_genomics/oncology/comparison.tsv) · [Regulatory evidence](comparative_genomics/oncology/regulatory_addendum.md)

## Exome and CNV

The query and truth are distinct published NA12878 callsets on GRCh37. The denominator is the intersection of autosomal capture and GIAB v3.3.2 high-confidence regions. Full REF spans must fit inside the mask. REF checking and left alignment of capture-overlapping PASS biallelic records against b37 yielded no mismatches or realignments. Independent BED boundary sweeping and normalized, unordered diploid-genotype set matching gave the following counts.

| Scope | TP | FP | FN | Precision | Recall | F1 |
|---|---:|---:|---:|---:|---:|---:|
| Exome SNV | 18,435 | 200 | 1,573 | 98.93% | 92.14% | 95.41% |
| Exome indel | 544 | 104 | 225 | 83.95% | 70.74% | 76.78% |
| Expanded panel SNV | 34 | 0 | 2 | 100% | 94.44% | 97.14% |
| Expanded panel indel | 0 | 1 | 0 | 0% | Undefined | 0% |

The expanded panel covers 96,565 confident bases. Exome genotype discrepancies contribute 77 SNV and 23 indel errors on each side. Exact matching leaves complex haplotype equivalence unresolved. BAM depth sensitivity used separate HCC1143 samples. [Counts](comparative_genomics/sequencing/normalized_crosscheck.tsv) · [Allele discrepancies](comparative_genomics/sequencing/discordances.tsv) · [Depth/GQ results](comparative_genomics/sequencing/metrics.tsv)

The CNV screen selects one of 48 settings using seed 504 and only the development F1. The reporting code evaluates the grid on seed 8841 but chooses the reported winner exclusively from development. Each cohort generates new target efficiencies, neutral controls and cases. Per-sample GC regression fits observed depth ratios without truth labels. The shared simulation assumptions limit external validity. On seed 8841 the false-positive count fell from 1,170 to 25. Both models had 145 true positives, of which **113 were the same targets**. Each model uniquely recovered another 32. Replaying the locked configuration on seeds 8842–8844 reduced false positives in every run but also reduced true positives. The locked recall across all four seeds ranged from 50.87% to 64.78%. [Seed results](comparative_genomics/sequencing/audit_seeds.tsv) · [Target identity cross-check](comparative_genomics/sequencing/audit_tp_overlap.tsv)

The CCLE comparison has 78 of 79 usable targets in one array segment and the remaining target in a second. Neither segment offers a neutral or loss comparison. Real-data CNV performance remains unmeasured. [Array intersection](comparative_genomics/sequencing/array_overlap_summary.tsv)

## References, figures and scope

The eight phase-two source hashes verified. A [reference link inventory](comparative_genomics/reference_audit.tsv) records HTTP checks. HTTP success checks availability only. Five historical FDA PDF paths returned 404. Current DailyMed entries supply alternate regulatory sources for [cobimetinib](https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=c387579e-cee0-4334-bd1e-73f93ac1bde6), [imatinib](https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=211ef2da-2868-4a77-8055-1cb2cd78e24b), [entrectinib](https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=c7c71b0c-2549-4495-86b6-c2807fa54908), [dabrafenib](https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=fee1e6b1-e1a5-4254-9f2e-a70e0f8dbdea) and [trastuzumab](https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=492dbdb2-077e-4064-bff3-372d6af0a7a2). GATK blocked automated HEAD requests. NCBI EFetch requires a parameterized request. The archive and [Phenopacket Store documentation](https://monarch-initiative.github.io/phenopacket-store/) establish dataset provenance. Full OMIM and independent full-text review of every underlying citation remain unavailable.

The five baseline figures are rendered directly from the phase-one queue, assertion types, SNV metrics, slice ratios and seeded CNV targets. The three expanded figures read recovery.tsv, assessments/comparison data and metrics.tsv. Byte-identical regeneration verifies rendering reproducibility. The original contributions are the bounded cross-source comparisons, perturbations, error audits and controlled simulations. The data resources, classification frameworks and established benchmarking approaches are credited in [data provenance](DATA_PROVENANCE.md) and [validation methods](VALIDATION_REPORT.md). These analyses provide research benchmarks only.

## Evidence for the six portfolio claims
### Rare disease

- Evaluated phenotype, allele-frequency and combined ranking on 30 published solved cases, including missing-term and noisy-term experiments with tie-aware Top-k and reciprocal rank. [Case manifest](comparative_genomics/rare/cases.tsv) · [Metrics](comparative_genomics/rare/recovery.tsv)
- Reconstructed 19 disease-specific germline expert assessments across five classification categories and audited phenotype rankings against independently cited disease annotations. [Evidence](comparative_genomics/rare/evidence_reviews.tsv) · [Citation audit](comparative_genomics/rare/distinct_publication.tsv)

### Somatic oncology

- Reviewed 23 molecular oncology contexts across 22 profiles, linking diagnostic, prognostic and therapeutic evidence to individual publications. [Contexts](comparative_genomics/oncology/contexts.tsv) · [Reports](comparative_genomics/oncology/report.md)
- Compared evidence-limited AMP tier assessments against accepted CIViC assertions and documented changes after targeted FDA-label review. [Initial comparison](comparative_genomics/oncology/comparison.tsv) · [Regulatory follow-up](comparative_genomics/oncology/regulatory_addendum.md)

### Exome and CNV

- Recalculated 18,435 exact SNV and 544 exact indel genotype matches in 31,561,021 confident capture bases and evaluated depth and genotype-quality sensitivity. [Independent counts](comparative_genomics/sequencing/normalized_crosscheck.tsv) · [Threshold results](comparative_genomics/sequencing/metrics.tsv)
- Compared 48 CNV screening configurations in simulated development data, evaluated a locked configuration on a separate simulated cohort, and intersected matched-cell-line exome ratios with orthogonal array segments. [Development](comparative_genomics/sequencing/cnv_development.tsv) · [Holdout](comparative_genomics/sequencing/cnv_holdout.tsv) · [Array comparison](comparative_genomics/sequencing/array_overlap_summary.tsv)

The later [repository quality review](FINAL_GITHUB_REVIEW.md) records the figure label updates and vector exports. All numerical results and baseline hashes survived the saved-input replay.
