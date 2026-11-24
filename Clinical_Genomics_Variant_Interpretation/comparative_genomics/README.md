# Comparative Analyses of Phenotype Ranking, Somatic Evidence and Genotype Recovery

Run the commands below from the repository root with the archived inputs described in the data provenance record.

## Can phenotype or population rarity identify the reported gene within a fixed candidate panel?

Thirty published cases have a solved interpretation, one causative variant record on GRCh38, and at least two observed HPO terms. Selection retains the first eligible case per gene in archive order and samples genes with seed 181, giving preference to requested genes when available. DMD met these criteria. CFTR, HBB, PAH and BRCA examples are represented in the separate evidence review. The same 30 alleles form every case’s candidate list. The causal gene is used only for panel assembly and scoring the result. Neither ClinVar class nor phenopacket pathogenicity class enters a ranking score.

Phenotype score is the maximum fraction of observed terms found within one HPO-annotated disorder for the candidate gene. Frequency score uses the higher available overall gnomAD exome/genome frequency, transformed as capped negative log frequency divided by eight. Fourteen alleles lacked retrievable counts and received a neutral score of 0.5. The combined score adds one quarter of the frequency score to phenotype overlap. Tied ranks are averaged over all possible tie orders. Perturbations remove half the terms or add three randomly selected other panel terms with a fixed seed.

With original terms, phenotype-only and combined Top-1 recovery were 95.0% and 93.3%. With half the terms retained they were 76.9% and 80.0%. Frequency-only Top-1 is 1/30 by construction because the same list is used for each case and each gene is the answer exactly once. This comparison measures phenotype discrimination within a selected panel. [All ranks and tie bounds](rare/ranks.tsv) · [Top-k and MRR](rare/recovery.tsv) · [Figure](rare/ranking.png)

Nineteen variant–disease evidence reconstructions cover all five germline classes. Each row identifies the SCV, applicable panel codes or historical framework, clinical context and rationale. The reconstructions reuse published evidence. They supply no independent classification accuracy estimate. [Evidence reviews](rare/evidence_reviews.tsv)

## How much of a reference AMP assessment can be supported by a bounded literature review?

The study selected 23 assertion contexts spanning 22 molecular profiles, taking distinct profiles in export order with an additional BRAF V600E prognostic context. The selector uses accepted, unflagged assertions and one unflagged linked PubMed evidence item of matching clinical type per context. It excludes the six initial reviews. The first assessment pass used the publication abstract or title and selected FDA labeling. The code records these judgments separately before joining CIViC reference tiers. Earlier inspection of the release and reference-based study selection prevent a blinded interpretation.

Eighteen contexts received initial tier assessments and five were deferred. Four of the 18 assigned contexts matched the broad reference tier and three matched tier and level. The large difference reflects restricted paper selection and uneven label coverage, especially for fusion diagnostics. Subsequent targeted regulatory review identified qualifying higher-tier indications for six contexts. Those follow-up judgments were recorded separately. The reports include drug sensitivity, KIT D816V resistance, diagnostic markers and a colorectal prognostic cohort. [Reports](oncology/report.md) · [Comparison](oncology/comparison.tsv) · [Regulatory follow-up](oncology/regulatory_addendum.md) · [Figure](oncology/tiers.png)

## How do depth, variant class and analysis scope affect error counts?

The autosomal ARUP/GIAB v3.3.2 comparison requires a PASS biallelic diploid nonreference genotype and the entire REF span within confident capture. Python matches alleles and unordered genotypes. A second implementation reads bcftools-normalized records and uses a separately constructed interval mask. Both recover the baseline SNV denominator and counts. Indels have lower precision and recall. Raising minimum query depth trades recall for SNV precision. Filtering a published VCF cannot simulate resequencing at lower depth. The expanded panel intersects 30 gene spans with the same exome capture design. [Metrics](sequencing/metrics.tsv) · [Independent counts](sequencing/normalized_crosscheck.tsv) · [Discordances](sequencing/discordances.tsv) · [Depth figure](sequencing/depth.png)

In a separate simulation, 48 configurations varied library normalization, selection of batch-matched controls, linear GC correction, minimum segment length and gain threshold. A five-target running median preceded thresholding. Development and held-out seeds generated distinct target efficiencies, control panels, noise and event locations. Each cohort contained 24 controls and 60 cases. The configuration with maximum development F1 was applied to the held-out cohort alongside a fixed comparator. Both recovered 145 altered targets. The selected configuration produced 25 false positives versus 1,170 for the comparator in this deliberately GC-confounded simulation. [Development grid](sequencing/cnv_development.tsv) · [Held-out results](sequencing/cnv_holdout.tsv)

An orthogonal CCLE segmentation file covered all 93 HCC1143 capture intervals. Among the 79 with usable sequencing ratios, 78 map to the same segment and one maps to another. Both segments have positive array means. These data cannot support a real gain/loss accuracy estimate. [Array overlap counts](sequencing/array_overlap_summary.tsv)

## Reproduce from the saved inputs

Install the environment specified in the root README. Populate the three study input directories using their README instructions and restore the comparative cache using the filenames in the provenance table. The replay rebuilds autosomal BED masks, indexes source VCFs as needed, normalizes alleles against indexed b37, and runs the tests and input checksums.

```bash
bash reproduce.sh
```

Source files and versions are documented in [data provenance](../DATA_PROVENANCE.md). The original project results and phase-two outputs can be checked separately in the [validation report](../VALIDATION_REPORT.md).

The [citation-restricted audit](rare/distinct_publication.tsv) retained only annotations with an explicit PMID different from the case PMID. Top-1 fell to 45.5% and MRR to 0.541. This changes both publication overlap and annotation coverage.

[Read-level depth sensitivity](sequencing/alignment_subsampling.tsv) uses seeded BAM subsampling, independently of the NA12878 VCF filtering experiment.
