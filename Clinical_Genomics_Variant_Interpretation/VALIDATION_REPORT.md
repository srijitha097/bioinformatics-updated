# Research validation · Phase two

## Rare disease

The baseline’s three synthetic scenarios yielded a queue without diagnostic truth. This phase used published causal assignments from the Phenopacket Store, withholding gene labels and variant classifications from scoring. One literature case per gene gives a different denominator from the baseline. Baseline-versus-revised ranking accuracy cannot be calculated. Top-k and reciprocal rank were recalculated from rank intervals, and the frequency-only means matched their analytical null values. Cases are independently documented, but their publications can overlap the literature underlying HPO annotations.

Danis and colleagues provide the underlying real case representations. Exomiser integrates sequence-level and phenotype evidence in genome-scale analysis. This experiment asks how missing terms and rarity affect a fixed, small gene panel. It does not benchmark Exomiser. [Phenopacket Store publication](https://pubmed.ncbi.nlm.nih.gov/39394689/) · [Exomiser study](https://pubmed.ncbi.nlm.nih.gov/39695184/) · [Recovery table](comparative_genomics/rare/recovery.tsv)

## Molecular oncology

The initial six assessments primarily reconstructed known contexts. Phase two expanded to 23 additional contexts. Reference tiers were excluded from the scoring table and joined afterward. Study selection still used CIViC evidence links, and some reference tiers had appeared during the baseline run. The exercise evaluates evidence sufficiency in a restricted review. Disagreement with an expert summary is not evidence that the expert tier is wrong. The follow-up regulatory review was informed by discrepancies and is recorded separately.

Li et al. define clinical significance by the strength and clinical setting of the evidence. Danos et al. describe CIViC assertions as syntheses of several evidence items and guidelines. A single selected study will often have less evidence than the assertion, especially for rare diagnostic fusions. [AMP/ASCO/CAP](https://pubmed.ncbi.nlm.nih.gov/27993330/) · [CIViC curation methodology](https://pmc.ncbi.nlm.nih.gov/articles/PMC6883603/) · [Context-by-context comparison](comparative_genomics/oncology/comparison.tsv)

## Exome and CNV

The baseline denominator was 31,561,021 confident autosomal capture bases. The second implementation reconstructed the same denominator by sweeping capture and confidence BED boundaries. Reference validation and left normalization used public GRCh37 sequence. bcftools norm accepted all 20,964 selected truth records and 24,004 selected query records without realignment or REF mismatch. Those counts precede the confidence intersection. Depth-filtered results use the query FORMAT/DP field. Missing DP/GQ fails an active threshold. Under exact matching, a wrong genotype contributes both FP and FN.

| Scope | TP | FP | FN | Precision | Recall | F1 |
|---|---:|---:|---:|---:|---:|---:|
| Exome SNV | 18,435 | 200 | 1,573 | 98.93% | 92.14% | 95.41% |
| Exome indel | 544 | 104 | 225 | 83.95% | 70.74% | 76.78% |
| 30-gene SNV | 34 | 0 | 2 | 100% | 94.44% | 97.14% |
| 30-gene indel | 0 | 1 | 0 | 0% | Undefined | 0% |

Of 1,573 SNV false negatives, 77 had a different genotype, 120 had a filtered query allele and 1,376 lacked an eligible query allele. Of 225 indel false negatives, 23 had a different genotype, one had a filtered query allele and 201 lacked an eligible query allele. Published variant-only records cannot distinguish uncovered bases from caller failures for absent alleles. At minimum query depth 30, SNV recall fell to 83.80% and indel recall to 63.72%. Haplotype-aware matching can rescue alternate representations, particularly for indels. [Krusche et al.](https://pubmed.ncbi.nlm.nih.gov/30858580/) · [Independent counts](comparative_genomics/sequencing/normalized_crosscheck.tsv)

The initial CNV simulation had 28 false gain targets and two false loss targets. Both two-target events were missed under five-target smoothing. In the new GC-confounded simulation, the selected development configuration used library normalization, pooled controls, GC regression, a three-target minimum run and a +0.5 gain threshold. Its held-out precision was 85.29%, recall 63.04%, and F1 72.50%. The baseline simulator used a different data-generating process. Its 54.55% precision and 78.26% recall are historical context only. Actual HCC1143 array comparisons have two overlapping segment levels and no suitable held-out truth events. [Held-out comparison](comparative_genomics/sequencing/cnv_holdout.tsv) · [GATK coverage modeling](https://gatk.broadinstitute.org/hc/en-us/articles/27007964420635-GermlineCNVCaller)

## Automated checks

All 15 tests passed, including seven new phase-two checks. They cover interval boundaries, rank ties, the frequency-only null, assessment joins, independent normalized counts, the frozen CNV selection, and baseline-file integrity. [Tests](tests/test_comparative_genomics.py)

## Read-level depth sensitivity

Seeded samtools subsampling at 25%, 50%, 75% and full input measured physical depth sensitivity separately from the VCF DP filter. Across three 50% replicates, the fraction of matched ICE target bases at ≥20× fell from 50.86% to 19.73% for the normal and from 57.18% to 27.14% for the tumor. The full-input endpoints reproduced phase-one coverage. These measurements concern the HCC1143 slices and cannot be used to explain individual NA12878 missed alleles. [Alignment subsampling results](comparative_genomics/sequencing/alignment_subsampling.tsv)

## Final independent checks

The saved-input replay reproduced the archived outputs. The audit checked packet eligibility, rank denominators, SCV versions, variant identities, reference classes and independently recomputed the overlap of correctly called CNV targets. [Audit code](comparative_genomics/final_checks.py) · [Results and discrepancy analysis](FINAL_AUDIT.md)

The [release review](FINAL_GITHUB_REVIEW.md) records the clean-directory replay, the version mismatch encountered during verification and the plotted-value checks.
