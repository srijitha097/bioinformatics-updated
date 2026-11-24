# Sequencing research report

## Published NA12878 exome calls against GIAB v3.3.2

| GRCh37 scope | Target bases | Confident target bases | TP | FP | FN | Precision | Recall | F1 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| ARUP SeqCap EZ autosomal capture | 34,479,984 | 31,561,021 | 18,435 | 200 | 1,573 | 98.93% | 92.14% | 95.41% |
| Ten hereditary cancer genes within capture | 44,260 | 41,713 | 17 | 0 | 1 | 100.00% | 94.44% | 97.14% |

The exome comparison found 77 shared alleles with different genotypes. Each contributes one false positive and one false negative. Exact genotype counts agree with the bcftools cross-check. The ten-gene subset contains only 18 truth SNVs. Nominal Wilson 95% intervals are 81.57–100.00% for precision and 74.24–99.01% for recall. These intervals are descriptive and do not capture genetic or regional dependence. This is a comparison of published WES calls within confident autosomal capture regions. It excludes indels and complex haplotypes and evaluates one reference genome. The panel reuses the WES capture and variant calls.

[Benchmark metrics](benchmark_metrics.tsv) · [Cross-check](benchmark_crosscheck.tsv) · [Figure](benchmark_metrics.png) · [Data release](https://doi.org/10.5281/zenodo.3636192)

## HCC1143 chromosome 20 exome slices

| Public cell line | ICE target bases | Mean depth | Bases ≥20× | Bases ≥30× |
|---|---:|---:|---:|---:|
| HCC1143_BL blood normal, SM-74NEG | 18,986 | 25.23× | 50.86% | 31.83% |
| HCC1143 tumor, SM-74P4M | 18,986 | 28.52× | 57.18% | 39.86% |

Seventy-nine of 93 targets reached a mean of at least 10× in both slices. The raw median log₂ tumor/normal ratio was 0.233, and centered retained ratios ranged from −0.736 to +0.749. Centering removes an overall shift within the slice. Coverage and local ratio variation can reflect capture efficiency, downsampling, or copy-number variation. There was no independent segment truth, purity adjustment, or absolute copy-number inference for these slices.

[Per-target results](real_copy_ratios.tsv) · [Figure](real_copy_ratios.png) · [GATK tumor/normal provenance](https://gatk.broadinstitute.org/hc/en-us/articles/360035890011--How-to-part-II-Sensitively-detect-copy-ratio-alterations-and-allelic-segments)

## Synthetic CNV benchmark

| Unit | TP | FP | FN | Precision | Recall |
|---|---:|---:|---:|---:|---:|
| Target and direction | 36 | 30 | 10 | 54.55% | 78.26% |

The simulation uses 16 copy-neutral reference samples and 12 held-out test samples across 100 targets. Eight test samples contain one implanted event and four are neutral. The screen missed both two-target events at their truth locations. Twenty-eight of 30 false-positive targets were gains. The five-target median attenuates short events, and the fixed gain threshold permits noise excursions. These metrics describe the specified seed and simulation assumptions, with no real-data accuracy estimate.

[Synthetic summary](synthetic/cnv_summary.json) · [Target table](synthetic/synthetic_targets.tsv) · [Segment table](synthetic/synthetic_segments.tsv) · [Figure](synthetic/cnv_examples.png)
