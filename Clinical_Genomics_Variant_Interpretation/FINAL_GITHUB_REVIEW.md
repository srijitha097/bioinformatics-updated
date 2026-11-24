# Clinical genomics repository release review

The three research workflows reproduced their archived numerical results in an isolated working directory. All 15 automated tests passed. The release set contains the code, compact results and eight figures in PNG and PDF. The principal reproducibility constraint is access to the dated input snapshots.

## README quality

Each study README identifies the research question, specimens or case source, unit of analysis, method, findings and limits. The root README provides environment setup and the saved-input replay command. Figure captions distinguish the three synthetic phenotype scenarios from the published cases, the cell-line observations from simulated CNVs, and the baseline ten-gene subset from the expanded panel. The germline evidence records retain their reviewer attribution. The six quantitative portfolio claims have file-level support in the [scientific audit](FINAL_AUDIT.md).

## Code quality

Reviewed all 25 Python modules and three shell scripts, including test modules and shared parsing code. The repository contains no R source. Removed unused imports and standardized formatting. Project-relative paths and tools on PATH replace workstation-specific paths. The replay builds the two autosomal BED masks from deposited intervals. BED coordinates are zero-based and half-open. A VCF allele at POS occupies the interval beginning at POS minus one and extending for the length of REF.

Fixed an interval-sweep edge case in which overlapping intervals from one BED could count twice toward a multi-file intersection. A regression fixture now includes an overlap that is present in only one input. The published genotype counts were unchanged. gnomAD aggregation distinguishes an absent variant from an API error. The environment records the Python package and alignment-tool versions used. [Interval cross-check](comparative_genomics/check_benchmark.py) · [Tests](tests/test_comparative_genomics.py) · [Environment](environment.yml)

## Figure-by-figure verification

Inspected all eight rendered PNGs, including axes, legends, labels, clipping and shaded truth regions. Captured the plotted Matplotlib data and compared the bars, lines or scatter coordinates with the indicated source records. The two CIViC panels use different horizontal scales, explicitly identified on the figure. Each PNG has an adjacent vector PDF. Local Markdown rendering checks covered 18 documents, 14 tables and 195 relative link or image references. All eight image embeds resolve to release files.

| Figure | Numerical check |
|---|---|
| [Synthetic phenotype overlap](rare_disease_variant_curation/results/phenotype_gene_scores.png) | Eleven gene–scenario bars match the grouped scores for the three scenarios in the [candidate queue](rare_disease_variant_curation/results/candidate_queue.tsv). |
| [Published case ranking](comparative_genomics/rare/ranking.png) | All nine method-by-condition MRR points match the [recovery table](comparative_genomics/rare/recovery.tsv), with 30 cases per point. |
| [Published CIViC tiers](oncology_somatic_variant_interpretation/results/tier_by_assertion_type.png) | Independent counts from the [assertions](oncology_somatic_variant_interpretation/results/accepted_assertions.tsv) give diagnostic 90/2, predictive 29/1 and prognostic 1/1 for I-A/II-C. |
| [Initial oncology assessments](comparative_genomics/oncology/tiers.png) | The [context table](comparative_genomics/oncology/comparison.tsv) gives 3 I-A, 1 I-B, 10 II-C, 4 II-D and 5 deferred. |
| [Exact SNV benchmark](exome_gene_panel_cnv_analysis/results/benchmark_metrics.png) | Four bars reproduce precision and recall from the [confusion counts](exome_gene_panel_cnv_analysis/results/benchmark_metrics.tsv). The baseline truth sets contain 20,008 exome and 18 panel SNVs. |
| [Query depth sensitivity](comparative_genomics/sequencing/depth.png) | At DP cutoffs 0, 10, 20 and 30, each point equals TP divided by the fixed truth denominator of 20,008 SNVs or 769 indels, with no GQ cutoff. [Metrics](comparative_genomics/sequencing/metrics.tsv) |
| [Cell-line depth ratios](exome_gene_panel_cnv_analysis/results/real_copy_ratios.png) | All 79 coordinates equal the interval midpoint in Mb and the centered log₂ ratio among 93 capture targets. [Target table](exome_gene_panel_cnv_analysis/results/real_copy_ratios.tsv) |
| [Synthetic CNV profiles](exome_gene_panel_cnv_analysis/results/synthetic/cnv_examples.png) | All 100 ratio values per displayed sample, four truth spans and the −0.50 and +0.35 thresholds match the [target data](exome_gene_panel_cnv_analysis/results/synthetic/synthetic_targets.tsv) and caller settings. |

## Scientific accuracy

- **Rare disease.** In 30 solved cases, phenotype-only Top-1/5/10 were 95%/100%/100%, with MRR 0.9690. The case PMID supported at least one causal-gene term in 21 cases. Requiring explicit different-PMID support gave 45.48%/60.10%/70.75% and MRR 0.5407. Seven cases retained no causal-gene term match. This intervention changes both overlap and annotation coverage. The 19 germline evidence reviews reconstruct expert decisions. [Ranks](comparative_genomics/rare/ranks.tsv) · [Citation audit](comparative_genomics/rare/distinct_publication.tsv) · [Germline evidence](comparative_genomics/rare/evidence_reviews.tsv)
- **Oncology.** Four of 18 assigned contexts matched the broad CIViC tier, three exactly matched tier and level, and five further contexts were deferred. The 14 assigned discrepancies involved four therapeutic, nine diagnostic and one prognostic interpretation. Restricted evidence selection omitted disease-specific labels, diagnostic criteria or guidelines used by the comparator. The subsequent regulatory review was informed by the discrepancies. [Case explanations](FINAL_AUDIT.md) · [Regulatory follow-up](comparative_genomics/oncology/regulatory_addendum.md)
- **Exome.** In 31,561,021 GRCh37 confident capture bases, normalized exact diploid matching yielded SNV TP/FP/FN of 18,435/200/1,573 and indel counts of 544/104/225. Precision/recall/F1 were 98.93%/92.14%/95.41% for SNVs and 83.95%/70.74%/76.78% for indels. The expanded panel had zero truth indels and undefined indel recall. [Independent counts](comparative_genomics/sequencing/normalized_crosscheck.tsv)
- **CNV.** The initial simulation yielded 36 TP, 30 FP and 10 FN. On the separate comparative holdout, the development-selected configuration reduced FP from 1,170 to 25. Both configurations recovered 145 targets, but only 113 were shared. On three further seeds the locked model lost sensitivity each time. Two positive CCLE segments across 79 usable targets cannot establish real-data CNV accuracy. [Initial simulation](exome_gene_panel_cnv_analysis/results/synthetic/cnv_summary.json) · [Overlap](comparative_genomics/sequencing/audit_tp_overlap.tsv) · [Seeds](comparative_genomics/sequencing/audit_seeds.tsv) · [Real-data summary](comparative_genomics/sequencing/array_overlap_summary.tsv)

## Automated tests

All 15 tests passed both in the saved-input replay and in a separate release-files-only copy with no raw data. They cover parser and evidence constraints, interval boundaries, exact genotypes, coverage calculations, ranking ties, frequency-only expectations, context joins, benchmark denominators, development-only selection and the baseline checksum manifest. Ruff checks for undefined names, unused imports and import order passed. The coverage command also reproduced the archived normal-sample coverage JSON independently of its caller.

## Reproducibility

The isolated replay used the proposed release files and private copies of the archived inputs. It executed the baseline workflows, comparative analyses, normalization, physical BAM subsampling, simulation, evidence joins and plotting. All 38 baseline artifact hashes and eight comparative input hashes verified. The resulting numerical tables and 16 figure exports matched the maintained copies. Manually encoded evidence judgments replay deterministically. Historical reports and the earlier demonstration retained in the baseline manifest were checked by hash.

An initial run resolved samtools and bcftools 1.21 ahead of 1.25. Six baseline checksums differed. The normal-sample mean depth was 25.2374× with 1.21 and 25.2333× with 1.25. The SNV confusion counts agreed in both runs. The replay now checks tool versions at startup. The successful run used the installed pinned tools with existing Python libraries. An independent Conda environment solve and a fresh acquisition of all mutable inputs remain unverified. [Replay](reproduce.sh) · [Input provenance](DATA_PROVENANCE.md)

## Files proposed for tracking

The [tracking manifest](TRACKING_MANIFEST.txt) names each release file. Git ignore rules omit raw downloads, intermediate per-base depths, local credentials and caches. The row-level legacy CCLE join stays local pending file-specific redistribution terms. The tracked aggregate records how many targets and segments contributed. The two historical demo files and original curation templates are distinct checksum fixtures. Dataset attribution and license details accompany the provenance record. The source tree has no Git history or configured publication remote.

## Repository size

122 proposed files total 1,388,215 bytes (1.32 MiB). This is the sum of regular-file lengths in the manifest, before Git compression. The local workspace occupies approximately 4.5 GiB on disk including raw inputs and reference indexes.

## Remaining limitations

All studies are research-only. The case panel is selected, the reference-informed oncology assessments lack blinded adjudication, and the sequencing benchmark uses one genome and an older GIAB truth release. Model tuning and holdout simulation share a generative model. The eight figures passed local rendering and value checks. GitHub-hosted rendering can only be confirmed after publication. Of 100 external URLs checked with HTTP HEAD, 95 returned 2xx, two GATK pages returned 403, two archived FDA PDFs returned 404, and the bare EFetch endpoint returned 405. Availability alone does not verify scientific content. Archived baseline reports retain expired FDA links, with replacement regulatory sources in the scientific audit. [Limitations](LIMITATIONS.md) · [External link checks](comparative_genomics/release_link_checks.tsv)

## Recommended commit message

`Add reproducible clinical genomics benchmarks, evidence reviews and audited figures`
